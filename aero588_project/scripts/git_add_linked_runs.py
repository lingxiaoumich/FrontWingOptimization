"""Stage the small files of run folders that are symlinks (runs on the 4T HD, linked into optimizations/).

git does not look past a symlink, so `git add` cannot see these runs. This writes each kept file as a blob
(git hash-object -w) and adds it to the index at its repository path, marked skip-worktree so the local
symlink does not show the files as deleted. The symlinks themselves go to .git/info/exclude.
Kept: everything except the per-run .gitignore patterns (flow solutions *.cgns/*.plt, optimizer databases
*.hst, evaluations/, cache_*/, written into each linked run as in runs 003-015) and files over --max-mb.
  python scripts/git_add_linked_runs.py [--dry] [--max-mb 10] [run names...]
"""
import argparse
import fnmatch
import os
import subprocess
from pathlib import Path

REPO = Path('/home/lingxiao77/AER588')
OPT = REPO/'aero588_project/optimizations'
RULES = '# Flow solutions, optimizer databases and render caches stay local\n*.cgns\n*.plt\n*.hst\nevaluations/\ncache_*/\n'

def ignored(rel, patterns):
    parts = rel.split('/')
    for p in patterns:
        if p.endswith('/'):
            if any(fnmatch.fnmatch(d, p[:-1]) for d in parts[:-1]): return True
        elif fnmatch.fnmatch(parts[-1], p): return True
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='*')
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--max-mb', type=float, default=10.)
    a = ap.parse_args()
    links = sorted(p for p in OPT.iterdir() if p.is_symlink() and p.is_dir() and (not a.runs or p.name in a.runs))
    entries, skipped, total = [], [], 0
    for link in links:
        real = link.resolve()
        gi = real/'.gitignore'
        if not gi.exists() and not a.dry: gi.write_text(RULES)
        patterns = [l.strip() for l in (gi.read_text() if gi.exists() else RULES).splitlines() if l.strip() and not l.startswith('#')]
        for root, dirs, files in os.walk(real):
            dirs[:] = [d for d in dirs if not ignored(os.path.relpath(os.path.join(root, d), real)+'/x', patterns)]
            for f in files:
                path = Path(root)/f; rel = os.path.relpath(path, real)
                if path.is_symlink() or ignored(rel, patterns) or f.endswith(('.pyc',)) or '__pycache__' in rel: continue
                size = path.stat().st_size
                if size > a.max_mb*1e6: skipped.append((size, f'{link.name}/{rel}')); continue
                entries.append((path, f'{link.relative_to(REPO)}/{rel}')); total += size
    print(f'{len(links)} linked runs, {len(entries)} files, {total/1e6:.1f} MB kept; {len(skipped)} over {a.max_mb:g} MB left out')
    for s, p in sorted(skipped, reverse=True)[:15]: print(f'  left out {s/1e6:6.1f} MB  {p}')
    if a.dry: return
    # write blobs in bulk, then the index entries
    paths = [str(p) for p, _ in entries]
    shas = subprocess.run(['git', '-C', str(REPO), 'hash-object', '-w', '--stdin-paths'], input='\n'.join(paths)+'\n',
                          capture_output=True, text=True, check=True).stdout.split()
    assert len(shas) == len(entries)
    info = ''.join(f"{'100755' if os.access(p, os.X_OK) else '100644'} {sha}\t{rel}\n" for (p, rel), sha in zip(entries, shas))
    subprocess.run(['git', '-C', str(REPO), 'update-index', '--add', '--index-info'], input=info, text=True, check=True)
    subprocess.run(['git', '-C', str(REPO), 'update-index', '--skip-worktree', '--stdin'], input='\n'.join(r for _, r in entries)+'\n',
                   text=True, check=True)
    ex = REPO/'.git/info/exclude'; have = ex.read_text() if ex.exists() else ''
    new = [f'/{l.relative_to(REPO)}' for l in links if f'/{l.relative_to(REPO)}\n' not in have]
    if new: ex.write_text(have+('' if have.endswith('\n') or not have else '\n')+'# runs on the 4T HD, symlinked; their small files are staged by scripts/git_add_linked_runs.py\n'+'\n'.join(new)+'\n')
    print(f'staged {len(entries)} files (skip-worktree); {len(new)} symlinks added to .git/info/exclude')

if __name__ == '__main__':
    main()
