#!/bin/bash
# Launch / attach / run commands in the MACH-Aero container for this repo (run inside WSL).
#
#   ./fwo.sh            -> interactive shell in the container (creates it on first use)
#   ./fwo.sh <cmd...>   -> run a command inside the container, in /project
#   ./fwo.sh --rm       -> remove the container (image is kept)
#
# Mounts V:\WORK\Research\AER588\Project_FWO\aero588_project at /project, the root the
# repo scripts expect. Same image as the course container; ~/.bashrc_mdolab is sourced
# explicitly because non-interactive shells skip it.
set -euo pipefail

IMAGE=mdolab/public:u22-gcc-ompi-stable
NAME=fwo
HOST_DIR=/mnt/v/WORK/Research/AER588/Project_FWO/aero588_project
MOUNT_DIR=/project

if [[ "${1:-}" == "--rm" ]]; then
  docker rm -f "$NAME" >/dev/null 2>&1 && echo "removed container $NAME" || echo "no container $NAME"
  exit 0
fi

if ! docker container inspect "$NAME" >/dev/null 2>&1; then
  docker create -it --name "$NAME" \
    --mount "type=bind,src=$HOST_DIR,target=$MOUNT_DIR" \
    -w "$MOUNT_DIR" --shm-size=2g "$IMAGE" /bin/bash >/dev/null
  echo "created container $NAME from $IMAGE"
fi
if [[ "$(docker container inspect -f '{{.State.Running}}' "$NAME")" != "true" ]]; then
  docker start "$NAME" >/dev/null
fi

# OpenMPI sees fewer physical cores than the VM's vCPUs; allow -n 8.
MPIENV=(-e OMPI_MCA_rmaps_base_oversubscribe=true)
if [[ $# -eq 0 ]]; then
  exec docker exec -it "${MPIENV[@]}" -w "$MOUNT_DIR" "$NAME" /bin/bash
else
  exec docker exec "${MPIENV[@]}" -w "$MOUNT_DIR" "$NAME" /bin/bash -c "source ~/.bashrc_mdolab; $*"
fi
