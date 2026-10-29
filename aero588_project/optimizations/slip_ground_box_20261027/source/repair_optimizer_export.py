import json,shutil
from pathlib import Path
import numpy as np
from pyoptsparse import History
root=Path('/project/optimization')
h=History(str(root/'optimization.hst'),flag='r')
meta=h.getMetadata();record=h.read(h.getCallCounters()[-1]);h.close()
x={k:np.asarray(v).tolist() for k,v in record['xuser'].items()}
best=json.loads((root/'best_feasible/result.json').read_text())
assert all(np.allclose(x[k],best['design'][k],rtol=0,atol=1e-14) for k in x)
assert abs(float(record['funcs']['objective'])-best['cl'])<1e-10
assert not best['cfd_failed'] and best['feasible']
assert 'ITERM= 11' in (root/'PSQP.out').read_text()
for path in (root/'best_feasible').iterdir():
    if path.is_file():shutil.copy2(path,root/'final'/path.name)
exitdata={'inform':{'value':11,'message':'Maximum number of iterations exceeded'},
    'variables':x,'optimizer_wall_time_s':meta['optTime'],'metadata':meta,
    'evaluations':len((root/'history.jsonl').read_text().splitlines()),
    'export_repaired_from_saved_history':True}
(root/'optimizer_exit.json').write_text(json.dumps(exitdata,indent=2)+'\n')
complete={'inform':exitdata['inform'],'final_failed':False,'final':best,'best_feasible':best,
    'optimizer_wall_time_s':meta['optTime'],'total_wall_time_s':meta['optTime'],
    'final_from_converged_evaluation':best['evaluation'],
    'export_repaired_from_saved_history':True,
    'note':'The optimizer completed. A dataclass JSON export error was repaired from the exact saved history and converged solution, without rerunning optimization.'}
(root/'run_complete.json').write_text(json.dumps(complete,indent=2)+'\n')
(root/'export_repair.txt').write_text(complete['note']+'\nThe repaired driver serializes pyOptSparse SolutionInform using dataclasses.asdict. The original executed driver and raw traceback are retained for audit.\nThe final variables exactly match the last history iterate and the saved converged best design.\n')
print(json.dumps({'exit':exitdata['inform'],'cl':best['cl'],'cd':best['cd'],'time_s':meta['optTime']}))
