#!/usr/bin/env python
## align multiple pdb files on a specified rigid body

import sys
import subprocess
import os

if len(sys.argv) < 5:
    sys.exit('USAGE: align_pdbs_on_body.py <pdb list file> <start residue> <end residue> <chain>\nlist file is one pdb per line, 1st one is used as the reference\nstart and end residues define the part of the model that is held static during the alignments')
chimerapath='/fbs/emsoftware2/LINUX/fbsmi/Chimera-1.11.2-linux/bin/chimera'
files = open(sys.argv[1],'r').readlines()
path = os.getcwd()

chicom = ['open #0 {0}/{1};wait;write #0 {0}/aligned_{1};wait;'.format(path,files[0].replace('\n',''))]
for i in files[1:]:
    chicom.append('open #1 {0}/{1};wait;'.format(path,i.replace('\n','')))
    chicom.append('mmaker #0:{0}-{1}.{2} #1:{0}-{1}.{2};wait;'.format(sys.argv[2],sys.argv[3],sys.argv[4]))
    chicom.append('write relative #0 #1 {0}/aligned_{1};wait;'.format(path,i.replace('\n','')))
    chicom.append('close #1;wait;')

chimeraout = open('chimeracommand.cmd','w')
chimeraout.write(''.join(chicom))
chimeraout.close()

run_chimera = subprocess.Popen('{0} --nogui {1}/chimeracommand.cmd'.format(chimerapath,path), shell=True, stdout=subprocess.PIPE)
screenbarf = run_chimera.stdout.read()
print (screenbarf)
