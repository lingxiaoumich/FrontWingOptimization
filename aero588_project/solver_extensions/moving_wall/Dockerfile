FROM mdolab/public:u24-gcc-ompi-stable
COPY patch_adflow.py prescribedWallMotion.F90 /tmp/moving-wall/
RUN /bin/bash -c 'source /home/mdolabuser/.bashrc_mdolab && python /tmp/moving-wall/patch_adflow.py && cd /home/mdolabuser/repos/adflow && make'
RUN /bin/bash -c 'source /home/mdolabuser/.bashrc_mdolab && pip install --no-deps --no-build-isolation /home/mdolabuser/repos/adflow'
