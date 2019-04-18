-=-=-=-=-=-=-=-=-=-=-=-=-
Movement analysis scripts
-=-=-=-=-=-=-=-=-=-=-=-=-

Analyze the motions of rigid bodies in MD trajectories or a series of structures thought to represent points along a motion.  Returns this for each body:

movement 1	distance (A)    direction
          	6.96 (-2.71, 0.43, 6.4)
rotation angles x	y	z
            	-1.6	-1.25	-4.82

movement 2	distance (A)    direction
          	2.16 (-2.12, 0.14, 0.43)
rotation angles	x	y	z
            	-9.65	-4.09	-1.18

movement 3	distance (A)    direction
          	3.61 (-1.82, 1.23, 2.86)
rotation angles	x	y	z
            	0.26	-1.2	-4.2

movement 4	distance (A)    direction
          	9.13 (-7.28, -3.22, -4.47)
rotation angles	x	y	z
            	4.14	-2.43	-8.57

Along with a .bild file for visulization.

Also reports the deviation of all Ca atoms between each frame and the last to analyze deformations of the body as it moves.

('Ca deviations between', ('P2_trajectory_frame_1.pdb', 'P2_trajectory_frame_2.pdb'))
Ca deviation (RMSD,min,max)   0.43	0.06	0.96

and generates a command to runin chimera to map the Ca deviations to the structure.

-=-=-= 
To use
-=-=-=

1) The models must have identical sequences if they don't, use the script seq_norm_pdbs.py to normalize them all to the common sequence

USAGE seq_norm_pdbs.py <pdb 1> <pdb 2> ... <pdb n>

2) Next align the models on the part of the sequence you want to hold constant using align_pdbs_on_body.py.

USAGE: align_pdbs_on_body.py <pdb list file> <start residue> <end residue> <chain>

List file is one pdb per line, 1st one is used as the reference

Start and end residues define the part of the model that is held static during the alignments

3) Finally run the script.  It requires a copy of UCSF chimera, edit the script to point to your installation.

USAGE: movement_analysis.py <body definition file> <models file>

body definition file = four columns text file
body_name       start_AA        end_AA      chain

models file -  text file list of models one per line -  models must be in order of the motion    

