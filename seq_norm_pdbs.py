#!/usr/bin/env python

# normalize the sequence of multiple pdb files...
# make sure they have the exact same atoms and no more

import sys

if len(sys.argv) < 2:
    sys.exit('\nUSAGE seq_norm_pdbs.py <pdb 1> <pdb 2> ... <pdb n>')

pdbs = sys.argv[1:]
    

def get_atoms(pdbdata):
    atomlines_dic = {}
    for line in pdbdata:
        if line[0:4] =='ATOM':
            linename = '{0}_{1}_{2}_{3}'.format(line[23:26].replace(' ',''),line[21],line[17:20].replace(' ',''),line[12:16].replace(' ',''))
            if 'HSD' in line:
                line = line.replace('HSD','HIS')
            atomlines_dic[linename] = line.replace('\n','')
    return(atomlines_dic)

def get_atom_ids(lines_dic):
    idlist = {}
    for i in lines_dic:
        idlist[int(lines_dic[i][4:11].replace(' ',''))] = [i]
    return(idlist)


keylist =[]
for i in pdbs:
    pdblines = open(i,'r').readlines()
    the_data = get_atoms(pdblines)
    keylist.append(set(the_data.keys()))

results_intersect = set.intersection(*keylist)

for i in pdbs:
    i=os.path.abspath(i.replace('\n',''))
    pdblines = open(i,'r').readlines()
    filename,path = i.split('/')[-1],'/'.join(i.split('/')[:-1])
    output = open('SN_{1}'.format(filename),'w')
    good_ids = []
    the_data = get_atoms(pdblines)
    for linekey in the_data:
        print (linekey)
        if linekey in results_intersect:
            good_ids.append(int(the_data[linekey][4:11].replace(' ','')))
    id_dic = get_atom_ids(the_data)
    good_ids.sort()

    for i in good_ids:
        print (i)
        print (id_dic[i])
        output.write('{0}\n'.format(the_data[id_dic[i][0]]))
    output.close()
