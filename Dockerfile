FROM ubuntu:20.04
RUN yes | apt-get update
RUN yes | apt-get upgrade
RUN yes | apt-get install --no-install-recommends git wget
RUN cd home; env GIT_SSL_NO_VERIFY=true git clone https://gitfront.io/r/mmcgibbon/cff5aca4804137fbf88d4c7357b32f0cd1c20ad1/SCORCH.git
RUN TZ=Europe/London
RUN cd /home/SCORCH; chmod +x setup.sh
RUN cd /home/SCORCH; ./setup.sh
CMD source /root/.bashrc && conda clean -a
CMD conda activate scorch; cd /home/SCORCH; python scorch.py -receptor examples/predocked_1a0q/1a0q_receptor.pdbqt -ligand examples/predocked_1a0q/ligands/1a0q_docked_ligand.pdbqt
