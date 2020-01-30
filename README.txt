-=-=-=-=-=-=-=-=-=-=-=-=-
Movement analysis scripts
-=-=-=-=-=-=-=-=-=-=-=-=-

Analyze the motions of rigid bodies in MD trajectories or a series of structures thought to represent points along a motion.

Returns a .bild file for visulization and Results/movement_data.txt which contains all transformation matrices for all the motions.

Also reports the deviation of all Ca atoms between each frame and the last to analyze deformations of the body as it moves and generates a command to run in chimera to map the Ca deviations to the structure.
The comands are found in Results/Deviation_analysis/


Centre of mass (COM) analysis compares different structures based on the locations of the rigid bodies and returns a correlation matrix for all structures.

-=-=-= 
To use
-=-=-=

1) The models must have identical sequences if they don't, use the script seq_norm_pdbs.py to normalize them all to the common sequence

USAGE seq_norm_pdbs.py <pdb 1> <pdb 2> ... <pdb n>

2) Next align the models on the part of the sequence you want to hold constant using align_pdbs_on_body.py.

USAGE: align_pdbs_on_body.py <pdb list file> <start residue> <end residue> <chain>

List file is one pdb per line, 1st one is used as the reference
The order of the pdbs in the list file is assumed to bethe order of the motions/steps in the trajectory
Start and end residues define the part of the model that is held static during the alignments

3) Finally run the script.  It requires a copy of UCSF chimera, edit the script to point to your installation.

USAGE: movement_analysis.py <body definition file> <models file>

body definition file = four columns text file
body_name       start_AA        end_AA      chain

models file -  text file list of models one per line -  models must be in order of the motion    

4) To compare structures with COM analysis:

USAGE: COM_analysis.py <body definition file> <models file> <max value for colour scaling - optional>

The body definition and pdb list files are as above
Max value for color scaling is the max used for the colors in the correlatoin matrix.  This can be set if you want to make matrices for several different datasets but want them all on the same absolute color scale. If this is left blank the maximum value for the set you are running on is used.

5) If you need to use the resulting bildfiles from either program in figures:
Use the script bildfile_figure.py to change the sizes of the spheres and replace the vectors connecting them with thicker cylinders, which will look much better.

USAGE: bildfile_figure.py <bild file> <cylinder 1 thickness> <cylinder 2 thickness> <sphere size>
cylinder 1 is the connections between each centre of mass
cylinder 2 is the three cylinders that represent the x,y,z axis on each sphere.  Set this to 0 to leave them as vectors.
