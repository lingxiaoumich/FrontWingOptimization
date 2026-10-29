"""Mirror stable optimization artifacts from local Docker staging to the project."""
import argparse
import json
from pathlib import Path
import re
import shutil
import time

p=argparse.ArgumentParser()
p.add_argument('--stage',type=Path,required=True)
p.add_argument('--destination',type=Path,required=True)
p.add_argument('--watch',action='store_true')
a=p.parse_args()
a.destination.mkdir(parents=True,exist_ok=True)
last_message=None
while True:
    for path in a.stage.rglob('*'):
        if not path.is_file() or path.name.endswith(('.hst','.lock')):
            continue
        stat=path.stat()
        if time.time()-stat.st_mtime < 4:
            continue
        target=a.destination/path.relative_to(a.stage)
        if target.exists() and target.stat().st_size==stat.st_size and target.stat().st_mtime==stat.st_mtime:
            continue
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target)
    records=[]
    history=a.stage/'history.jsonl'
    if history.exists():
        for line in history.read_text().splitlines():
            try:records.append(json.loads(line))
            except json.JSONDecodeError:pass
    iteration=None
    log=a.stage/'PSQP.out'
    if log.exists():
        values=re.findall(r'NIT=\s*(\d+)',log.read_text(errors='replace'))
        if values:iteration=int(values[-1])
    if records:
        latest=records[-1]
        message={k:latest[k] for k in ['evaluation','cl','cd','feasible','cfd_failed','elapsed_wall_time_s']}
        message.update(major_iteration=iteration,
                       downforce_change_percent=100*(latest['cl']/records[0]['cl']-1),
                       failed_cfd_evaluations=sum(r['cfd_failed'] for r in records))
        if message!=last_message:
            print(json.dumps(message),flush=True)
            last_message=message
            (a.destination/'progress.json').write_text(json.dumps(message,indent=2)+'\n')
    complete=(a.stage/'run_complete.json').exists()
    if complete:
        print('OPTIMIZER_RETURNED',flush=True)
        break
    if not a.watch:break
    time.sleep(15)
