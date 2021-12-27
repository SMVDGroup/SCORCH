=======================================================================
GWOVina 1.0
Michael K. M. Wong & Shirley W. I. SIU                                   
                                                                       
Computational Biology and Bioinformatics Lab (CBBio)                        
University of Macau                                                    
                                                                       
http://cbbio.cis.umac.mo/software/gwovina
=======================================================================

0. WHAT IS IT?
   A fast protein-ligand docking tool based on the implementation of 
   AutoDock Vina and grey wolf optimization (GWO) method.

1. INSTALLATION
   For successful compilation, please install Boost (version 1.59.0) 
   from http://www.boost.org. For preparing molecules for docking, 
   please install AutoDockTools (ADT) from http://mgltools.scripps.edu.

   The installation basically follows the installation of AutoDock Vina. 
   The steps are simple:

     a. unpack the files
     b. cd gwovina-1.0/build/<your-platform>/release
     c. modify Makefile to suit your system setting
     d. type "make" to compile
    
    The binary gsovina will be generated at the current directory. You can 
    copy this binary to a directory in your PATH e.g. /usr/local/bin, or add
    the path of the current directory to your PATH.

2. RUNNING GWOVINA

   You can run gwovina as the way you run vina but an additional parameter
   to specify the population size can be given:

   % <path-to-gwovina>/gwovina

   GWO parameters (optional):
       --num_wolves arg (=12)      Number of wolves

   For example, docking Kifunensine in the Mannosidase enzyme (PDBID 1ps3 from
   the PDBbind v2012 dataset) using GWOVina:

   % <path-to-AutoDockTools>/prepare_ligand4.py -l 1ps3_ligand.mol2 \
     -o 1ps3_ligand.pdbqt -A 'hydrogens' -U 'nphs_lps_waters'

   % <path-to-AutoDockTools>/prepare_receptor4.py -r 1ps3_protein.pdb \
     -o 1ps3_protein.pdbqt -A 'hydrogens' -U 'nphs_lps_waters' 

   % <path-to-gwovina>/gwovina \
     --receptor 1ps3_protein.pdbqt --ligand 1ps3_ligand.pdbqt \
     --center_x  31.951 --center_y 65.5053 --center_z 7.63888 \
     --size_x    33.452 --size_y   27.612  --size_z   35.136  \ 
     --cpu 8  
   

3. DEVELOP GWOVINA

   If you are interested in the source code of GWOVina for any academic
   purposes, please note that the following files were newly developed   
   in our work or modified based on GWOVina:
	src/main/main.cpp
        src/lib/gwo.h
        src/lib/gwo.cpp
        src/lib/parallel_gwo.h
        src/lib/parallel_gwo.cpp
        src/lib/gwo_mutate.h
        src/lib/gwo_mutate.cpp

4. CITATION
   Please cite our paper if you have used any version of GWOVina for your work. 

   Kin Meng Wong, Hio Kuan Tai, and Shirley W. I. Siu*. Evaluation of Grey Wolf Optimization Algorithm on Rigid and Flexible Receptor Docking. (Submitted)
   
   Please check out our homepage for the updated citation.

5. CONTACT US
    
   Developers: 
     Michael Wong (emailtoming15@hotmail.com)

   Project P.I.: 
     Shirley W. I. Siu (shirleysiu@umac.mo)

   Computational Biology and Bioinformatics Lab, University of Macau
   http://cbbio.cis.umac.mo
   http://www.cis.umac.mo/~shirleysiu
