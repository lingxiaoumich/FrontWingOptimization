"""Read surface y-plus with the installed CGNS C API (no ADF ABI assumptions)."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('surface',type=Path)
p.add_argument('--library',default='/home/mdolabuser/packages/CGNS-4.3.0/opt-gcc/lib/libcgns.so')
p.add_argument('--cgsize-bits',type=int,choices=[32,64],default=32)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
lib=ct.CDLL(args.library)
lib.cg_get_error.restype=ct.c_char_p
size_t=ct.c_int32 if args.cgsize_bits==32 else ct.c_int64
def call(name,*items):
    code=getattr(lib,name)(*items)
    if code:
        raise RuntimeError(name+': '+lib.cg_get_error().decode())
fn=ct.c_int()
call('cg_open',str(args.surface).encode(),0,ct.byref(fn))
stats={}
try:
    nb=ct.c_int();call('cg_nbases',fn,ct.byref(nb))
    for b in range(1,nb.value+1):
        nz=ct.c_int();call('cg_nzones',fn,b,ct.byref(nz))
        for z in range(1,nz.value+1):
            name=ct.create_string_buffer(33);sizes=(size_t*9)()
            call('cg_zone_read',fn,b,z,name,sizes)
            ns=ct.c_int();call('cg_nsols',fn,b,z,ct.byref(ns))
            for sol in range(1,ns.value+1):
                dim=ct.c_int();dims=(size_t*12)()
                call('cg_sol_size',fn,b,z,sol,ct.byref(dim),dims)
                shape=tuple(dims[:dim.value]);count=int(np.prod(shape))
                nf=ct.c_int();call('cg_nfields',fn,b,z,sol,ct.byref(nf))
                for field in range(1,nf.value+1):
                    fname=ct.create_string_buffer(33);kind=ct.c_int()
                    call('cg_field_info',fn,b,z,sol,field,ct.byref(kind),fname)
                    if fname.value.decode().lower() not in {'yplus','y_plus'}:continue
                    lo=(size_t*dim.value)(*[1]*dim.value)
                    hi=(size_t*dim.value)(*shape)
                    values=np.empty(count,dtype=np.float64)
                    call('cg_field_read',fn,b,z,sol,fname.value,4,lo,hi,ct.c_void_p(values.ctypes.data))
                    stats[name.value.decode()]={'yplus_max':float(values.max()),
                        'yplus_median':float(np.median(values)),'points':count}
finally:
    call('cg_close',fn)
assert stats,'Surface file contains no y-plus fields'
args.output.write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
