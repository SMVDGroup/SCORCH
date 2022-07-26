FROM ubuntu:20.04
RUN yes | apt-get update
RUN yes | apt-get upgrade
RUN yes | apt-get install --no-install-recommends git wget
RUN git config --global http.sslverify false; cd home; git clone https://github.com/SMVDGroup/SCORCH.git
RUN TZ=Europe/London
RUN cd /home/SCORCH; chmod +x setup.sh
RUN cd /home/SCORCH; ./setup.sh
CMD source /root/.bashrc && conda clean -a
CMD conda activate scorch; cd /home/SCORCH; python scorch.py -receptor examples/predocked_1a0q/1a0q_receptor.pdbqt -ligand examples/predocked_1a0q/ligands/1a0q_docked_ligand.pdbqt
