"""Extract core CGNS fields with the container's matching 32-bit CGNS C API."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root
lib=ct.CDLL('/home/mdolabuser/packages/CGNS-4.3.0/opt-gcc/lib/libcgns.so')
lib.cg_get_error.restype=ct.c_char_p
size_t=ct.c_int32

def call(name,*args):
    code=getattr(lib,name)(*args)
    if code:raise RuntimeError(name+': '+lib.cg_get_error().decode())

def extract(volume):
    fn=ct.c_int();call('cg_open',str(volume).encode(),0,ct.byref(fn))
    data={};stats={}
    try:
        nz=ct.c_int();call('cg_nzones',fn,1,ct.byref(nz))
        for z in range(1,nz.value+1):
            name=ct.create_string_buffer(33);sizes=(size_t*9)()
            call('cg_zone_read',fn,1,z,name,sizes)
            tag='near' if b'airfoil' in name.value else 'background'
            vdims=tuple(sizes[:3]);cdims=tuple(sizes[3:6])
            project=(lambda v:v[:,0,:]) if tag=='near' else (lambda v:v[:,:,0])
            def read(field,coord=False):
                dims=vdims if coord else cdims
                lo=(size_t*3)(1,1,1);hi=(size_t*3)(*dims)
                buf=np.empty(int(np.prod(dims)),dtype=np.float64)
                params=(fn,1,z,field.encode(),4,lo,hi,ct.c_void_p(buf.ctypes.data))
                if coord:call('cg_coord_read',*params)
                else:call('cg_field_read',fn,1,z,1,*params[3:])
                return project(buf.reshape(dims,order='F'))
            for name in ['CoordinateX','CoordinateY']:
                data[tag+'_'+name]=read(name,True)
            for name in ['Pressure','Density','VelocityX','VelocityY','Mach','CoefPressure','Iblank']:
                data[tag+'_'+name]=read(name)
            xyzx=data[tag+'_CoordinateX'];xyzy=data[tag+'_CoordinateY']
            # Signed areas of all four corner orientations, excluding z span.
            corners=[np.stack([xx,yy],axis=-1) for xx,yy in [
                (xyzx[:-1,:-1],xyzy[:-1,:-1]),(xyzx[1:,:-1],xyzy[1:,:-1]),
                (xyzx[1:,1:],xyzy[1:,1:]),(xyzx[:-1,1:],xyzy[:-1,1:])]]
            corner_cross=[]
            for k in range(4):
                v=corners[(k+1)%4]-corners[k];w=corners[(k-1)%4]-corners[k]
                corner_cross.append(v[...,0]*w[...,1]-v[...,1]*w[...,0])
            area=np.stack(corner_cross)
            # Native near-grid index order is i,z,j, so planar i,j is reversed.
            signed=-area if tag=='near' else area
            stats[tag]={'grid_shape':vdims,'cell_shape':cdims,
                'minimum_corner_area_m2':float(signed.min()),
                'nonpositive_corner_areas':int((signed<=0).sum()),
                'minimum_grid_y_m':float(xyzy.min()),
                'iblank_counts':{str(int(k)):int(v) for k,v in zip(*np.unique(data[tag+'_Iblank'],return_counts=True))}}
    finally:call('cg_close',fn)
    return data,stats

settings=json.loads((root/'settings.json').read_text())
mach=settings['flow_derived']['mach'];gamma=1.4;pref=settings['pressure_Pa']
q=.5*gamma*pref*mach**2;p0=pref*(1+.2*mach**2)**3.5
for label in ['baseline','final','best_feasible']:
    folder=root/label
    if not (folder/'flow_000_vol.cgns').exists():continue
    fields,stats=extract(folder/'flow_000_vol.cgns')
    for tag in ['near','background']:
        m=fields[tag+'_Mach'];pressure=fields[tag+'_Pressure']
        fields[tag+'_TotalPressure_Pa']=pref*pressure*(1+.2*m**2)**3.5
        ib=fields[tag+'_Iblank'];valid=np.isin(ib,[1,-1])
        cpcalc=(pressure-1)*pref/q
        stats[tag]['cp_reconstruction_max_error']=float(np.max(np.abs(cpcalc[valid]-fields[tag+'_CoefPressure'][valid])))
        assert stats[tag]['cp_reconstruction_max_error']<1e-4,stats[tag]
        assert stats[tag]['nonpositive_corner_areas']==0,stats[tag]
    fields.update({'gamma':gamma,'pressure_reference_Pa':pref,'dynamic_pressure_inf_Pa':q,
                   'total_pressure_inf_Pa':p0,'mach_inf':mach})
    np.savez_compressed(folder/'fields.npz',**fields)
    rows=[]
    for line in (folder/'flow_000_slices.dat').read_text().splitlines():
        v=line.split()
        if len(v)==9:
            try:rows.append([float(x) for x in v])
            except ValueError:pass
    surface=np.asarray(rows)
    assert surface.shape[0]>200
    shape=np.loadtxt(folder/'surface_geometry.csv',delimiter=',',skiprows=1)
    le=shape[136,:2];te=(shape[0,:2]+shape[272,:2])/2;ch=te-le
    unit=ch/np.linalg.norm(ch);normal=np.array([-unit[1],unit[0]])
    surface[:,3]=(surface[:,:2]-le)@unit/.5
    surface[:,4]=(surface[:,:2]-le)@normal/.5
    np.savetxt(folder/'surface_cp.csv',surface,delimiter=',',
       header='x_m,y_m,z_m,chord_x_over_c0,chord_y_over_c0,z_over_c0,Cp,Cf,yplus',comments='')
    result=json.loads((folder/'result.json').read_text())
    active_mach=np.concatenate([fields[tag+'_Mach'][np.isin(fields[tag+'_Iblank'],[1,-1])] for tag in ['near','background']])
    active_p0=np.concatenate([fields[tag+'_TotalPressure_Pa'][np.isin(fields[tag+'_Iblank'],[1,-1])]/p0 for tag in ['near','background']])
    result.update({'surface_cp_min':float(surface[:,6].min()),'surface_cp_max':float(surface[:,6].max()),
        'maximum_active_mach':float(active_mach.max()),
        'minimum_total_pressure_ratio':float(active_p0.min()),'maximum_total_pressure_ratio':float(active_p0.max()),
        'surface_yplus_max':float(surface[:,8].max()),
        'surface_yplus_median':float(np.median(surface[:,8])), 'grid_validation':stats})
    (folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(label,json.dumps(stats),flush=True)
