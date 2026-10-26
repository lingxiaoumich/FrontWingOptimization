#!/bin/bash
# From Git Bash / PowerShell on Windows: run a script file inside the fwo container.
#   bash env/run_in.sh env/scratch/x.sh      (path relative to aero588_project)
exec wsl.exe -d Ubuntu-24.04 -- bash -c "cd /mnt/v/WORK/Research/AER588/Project_FWO/aero588_project/env && sed -i 's/\r\$//' ../$1 && ./fwo.sh bash $1"
