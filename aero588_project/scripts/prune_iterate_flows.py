"""Remove the saved flow files of intermediate PSQP iterates; keep what recreates them.

  python scripts/prune_iterate_flows.py optimizations            # list what would go
  python scripts/prune_iterate_flows.py optimizations --delete   # remove it

In every run or chain stage under the given folders (any folder with an iterates/ subfolder)
the first, the best and the last iterate keep their flow files: the best is the design to
report, the last is what a remesh or restart starts from. Every other iterate loses
flow_000_vol.cgns (14 MB) and, where they exist, flow_000_surf.cgns and flow_000_surf.plt.
Its design.json, result.json, surface_geometry.csv and flow_000_slices.dat stay, so
regenerate_iterate_flow.py can solve the design again. A folder without design.json is an
unfinished write (a crash); its flow file goes too. Runs whose solver is running are skipped.
Each pruned iterates/ folder gets a PRUNED.txt saying what was kept.
"""
import argparse
import json
import subprocess
import time
from pathlib import Path

FLOW = ('flow_000_vol.cgns', 'flow_000_surf.cgns', 'flow_000_surf.plt')
p = argparse.ArgumentParser()
p.add_argument('folders', type=Path, nargs='+')
p.add_argument('--delete', action='store_true', help='remove the files; without it only the list is printed')
a = p.parse_args()

number = lambda d: int(d.name.split('_')[1])
size = lambda d: sum((d/f).stat().st_size for f in FLOW if (d/f).exists())
total = kept_total = 0
for top in a.folders:
    for its in sorted(top.resolve().rglob('iterates')):
        run = its.parent
        if subprocess.run(['pgrep', '-f', f'[p]ython -u {run}/source/'], capture_output=True).stdout.strip():
            print(f'{run.name}: its solver is running; skipped'); continue
        dirs = sorted((d for d in its.glob('eval_*') if d.is_dir()), key=number)
        whole = [d for d in dirs if (d/'design.json').exists() and (d/'result.json').exists()]
        broken = [d for d in dirs if d not in whole and size(d)]
        flows = [d for d in whole if (d/'flow_000_vol.cgns').exists()]
        if not flows and not broken: continue
        res = {d: json.loads((d/'result.json').read_text()) for d in flows}
        pool = [d for d in flows if res[d]['feasible']] or flows  # saved iterates are converged solves
        keep = {flows[0], flows[-1], min(pool, key=lambda d: res[d]['cl'])} if flows else set()
        go = [d for d in flows if d not in keep]+broken
        freed = sum(size(d) for d in go); total += freed; kept_total += sum(size(d) for d in keep)
        if not go: continue
        print(f"{str(run.relative_to(top.resolve().parent)):78s} {len(flows):3d} flows, keep {', '.join(d.name[5:] for d in sorted(keep, key=number)):14s} "
              f"remove {len(go):3d} ({freed/1e9:.2f} GB){' incl. unfinished '+', '.join(d.name for d in broken) if broken else ''}")
        if a.delete:
            for d in go:
                for f in FLOW: (d/f).unlink(missing_ok=True)
                if d in broken and not any(d.iterdir()): d.rmdir()
            with (its/'PRUNED.txt').open('a') as f:
                f.write(f"{time.strftime('%c')}: prune_iterate_flows.py removed the flow files ({', '.join(FLOW)}) of {len(go)} iterates "
                        f"({freed/1e9:.2f} GB). Kept with their flow: {', '.join(d.name for d in sorted(keep, key=number))} (first, best, last). "
                        'Every iterate keeps design.json, result.json and surface_geometry.csv; '
                        'scripts/regenerate_iterate_flow.py <iterate folder> solves a design again and puts its flow back.\n')
print(f"{'removed' if a.delete else 'would remove'} {total/1e9:.2f} GB; flow files kept in the first, best and last iterates: {kept_total/1e9:.2f} GB")
