import ast,json
from pathlib import Path
from dataclasses import asdict,is_dataclass
import numpy as np
from pyoptsparse.pyOpt_solution import SolutionInform
source=Path('/project/run_slip_optimization.py').read_text()
fn=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='json_default')
exec(compile(ast.Module(body=[fn],type_ignores=[]),'json_default_under_test','exec'))
encoded=json.dumps({'inform':SolutionInform(11,'Maximum number of iterations exceeded'),
                    'array':np.array([1.,2.]),'scalar':np.int64(4)},default=json_default)
r=json.loads(encoded)
assert r['inform']['value']==11 and r['array']==[1.,2.] and r['scalar']==4
Path('/project/optimization/source/export_serialization_check.json').write_text(json.dumps({'passed':True,'round_trip':r},indent=2)+'\n')
print('Serialization round-trip passed for actual pyOptSparse SolutionInform, NumPy array and scalar.')
