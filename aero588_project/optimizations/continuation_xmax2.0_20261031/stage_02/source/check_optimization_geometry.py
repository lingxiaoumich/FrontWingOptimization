"""Check the custom actual-surface box/thickness Jacobians at a selected design."""
from pathlib import Path
setup=Path(__file__).with_name('run_slip_optimization.py').read_text().split('options={')[0]
exec(compile(setup,'shared_geometry_setup','exec'))
selected=root/'source/design_to_check.json'
x=({k:np.array(v) for k,v in json.loads(selected.read_text()).items()}
   if selected.exists() else copy.deepcopy(xinitial))
geo.setDesignVars(x);pts=geo.update('surface_check');jac=geometry_sens(pts)
# Use a nonzero design to avoid the baseline's coincident upper/lower knots,
# at which an exact piecewise-linear maximum can have multiple active branches.
if not selected.exists():
    for j,name in enumerate(shape_groups):
        x[name]=.002*np.sin(np.arange(len(x[name]))*.71+.37*(j+1))
    x['incidence']=np.array([-3.02])
geo.setDesignVars(x);pts=geo.update('surface_check');jac=geometry_sens(pts)
tests=[('incidence',np.ones(1),1e-4,'incidence')]
for name in shape_groups:
    tests.extend([(name,np.sin(np.arange(len(x[name]))+.37),1e-7,'varying_'+name),
                  (name,np.tile([-1.,1.],nstation),1e-7,'upper_lower_'+name)])
errors={}
for key,direction,h,label in tests:
    direction=direction/np.linalg.norm(direction)
    xp=copy.deepcopy(x);xm=copy.deepcopy(x)
    xp[key]+=h*direction;xm[key]-=h*direction
    geo.setDesignVars(xp);vp=geometry_values(geo.update('surface_check'))
    geo.setDesignVars(xm);vm=geometry_values(geo.update('surface_check'))
    for name in jac:
        fd=(vp[name]-vm[name])/(2*h)
        exact=np.asarray(jac[name][key])@direction
        errors[label+'_'+name]=float(np.max(np.abs(fd-exact)))
assert max(errors.values())<1e-6,errors
dump(root/'geometry_derivative_check_extended.json',{'maximum_absolute_errors':errors,'design':realdict(x)})
if comm.rank==0:print(json.dumps(errors,indent=2))
