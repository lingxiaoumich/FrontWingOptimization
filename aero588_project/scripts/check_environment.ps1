param([int]$Ranks = 2)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$localWork = Join-Path $env:USERPROFILE 'Documents\Codex\aero588_solver_work'
New-Item -ItemType Directory -Path (Join-Path $localWork 'scripts'),(Join-Path $localWork 'reports') -Force | Out-Null
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'verify_environment.py') -Destination (Join-Path $localWork 'scripts\verify_environment.py') -Force
$linuxWork = ((& wsl -d Ubuntu --exec wslpath -a $localWork) -join '').Trim()
if ($LASTEXITCODE -ne 0) { throw 'Ubuntu startup or path conversion failed.' }
& wsl -d Ubuntu -u root --exec /usr/sbin/service docker start
if ($LASTEXITCODE -ne 0) { throw 'Docker startup failed.' }
$mount = "type=bind,source=$linuxWork,target=/project"
& wsl -d Ubuntu -u root --exec docker run --rm --mount $mount mdolab/public:u24-gcc-ompi-stable /bin/bash -c "source /home/mdolabuser/.bashrc_mdolab && mpirun -np $Ranks python /project/scripts/verify_environment.py"
if ($LASTEXITCODE -ne 0) { throw 'MACH-Aero environment verification failed.' }
Copy-Item -LiteralPath (Join-Path $localWork 'reports\environment_check.json') -Destination (Join-Path $projectRoot 'reports\environment_check.json') -Force
Write-Host 'Environment check passed; report saved in the project reports folder. No CFD run.'

