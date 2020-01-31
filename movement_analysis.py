#!/usr/bin/env python

#### update the path to chimera here ############
chimerapath='/Applications/Chimera.app/Contents/MacOS/chimera'
#################################################

import os
import sys
import math
import subprocess
import numpy as np

vers = 0.5
if os.path.isfile(chimerapath) == False:
	sys.exit("\nERROR: can't find UCSF Chimera at {0}\nupdate the path by editing the script".format(chimerapath))
alignedpath = '{0}/aligned'.format(os.getcwd())
resultspath = '{0}/Results'.format(os.getcwd())
tmppath = '{0}/TMP'.format(os.getcwd())
devapath = '{0}/Results/deviation_analysis'.format(os.getcwd())

errormsg = '''
USAGE: movement_analysis.py <body definition file> <models file>
---
body definition file = four columns text file
body_name       start_AA        end_AA      chain

models file -  text file list of models one per line -  models must be in order of the motion    

---'''

if len(sys.argv) < 3:
    sys.exit(errormsg)
if os.path.isfile(sys.argv[1]) == False or os.path.isfile(sys.argv[2]) == False:
    sys.exit(errormsg)

def make_directories():
    if os.path.isdir(alignedpath) == False:
        subprocess.call(['mkdir',alignedpath])
    if os.path.isdir(resultspath) == False:
        subprocess.call(['mkdir',resultspath])
    if os.path.isdir(tmppath) == False:
        subprocess.call(['mkdir',tmppath])
    if os.path.isdir(devapath) == False:
        subprocess.call(['mkdir',devapath])

def get_files():
    print('-------------------\nfiles to operate on\n-------------------')
    pdbfiles = open(sys.argv[2],'r').readlines()
    pdblist = []
    for i in pdbfiles:
        pdblist.append(i.replace('\n',''))
        print(i.replace('\n',''))
    return(pdblist)

def get_bodies():
    allbodies = []
    bodydeffile = open(sys.argv[1],'r').readlines()
    print('\n----------------\nbody definitions\n----------------\nname    start   end     chain')
    for i in bodydeffile:
        line = i.split()
        print('{0}\t{1}\t{2}\t{3}'.format(line[0],line[1],line[2],line[3]))
        allbodies.append((line[0],range(int(line[1]),int(line[2])),line[3]))
    return(allbodies)

def slice_n_save(pdbfile,potrarange):
    '''save each POTRA from a pdb as an individual pdb file'''
    pdbdata = open(pdbfile,'r').readlines()
    goodlines = []
    
    # get the lines
    for line in pdbdata:
        try:
                if int(line[23:26]) in potrarange[1] and line[21] == potrarange[2] and line[0:4] =='ATOM' and line[13:16] == 'CA ':
                    goodlines.append(line)
        except:
            pass

    #write the output
    filename ='{0}_{1}'.format(potrarange[0],pdbfile)
    output = open('{0}/{1}'.format(alignedpath,filename),'w')
    output.write('{0}'.format(goodlines[0]))
    for i in goodlines[1:]:
        output.write('\n{0}'.format(i))
    output.close()
    print('{0}'.format(filename))
    return(filename)

def write_chimera_script(path,expPDB,refPDB):
    output.write('open #0 {0}/{1};'.format(path,expPDB))
    output.write('wait;')
    output.write('open #1 {0}/{1};'.format(path,refPDB))
    output.write('wait;')
    output.write('match #0 #1 showMatrix true;')
    output.write('wait;')
    output.write('write relative #1 #0 {0}/{1}_ft_{2};'.format(alignedpath,expPDB.split('.')[0],refPDB))
    output.write('wait;')    
    output.write('close #0 #1;')
    output.write('wait;')

def make_sub_pdbs(pdblist,allbodies):
    outlist = []        # list of sliced pdb files made
    print('''\n-------------------------\nmaking sub-body pdb files\n-------------------------''')
    for file in pdblist:
        for potra in allbodies:
            outlist.append((slice_n_save(file,potra),potra[0]))
    return(outlist)

def order_bodies(subpdblist):
    bodies = {}
    for i in subpdblist:
        if i[1] not in bodies:
            bodies[i[1]] = [i[0]]
        else:
            bodies[i[1]].append(i[0])
        bkeys = bodies.keys()
    bkeys.sort()
    
    print('''\n----------------------\nbody models (in order)\n----------------------''')
    for i in bkeys:
        print (i,bodies[i])
    return(bodies,bkeys)

def make_body_pairs(bodkeys,bodydic):
    bodypairs = []
    for i in bodkeys:
        for j in range(len(bodydic[i])):
            try:
                bodypairs.append((bodydic[i].pop(0),bodydic[i][0]))
            except:
                pass
            
    print('''\n------------------------------------------------------------------\nsets for comparison (reference pdb / mobile pdb)\n------------------------------------------------------------------''')
    for i in bodypairs:
        print ('{0} / {1}'.format(i[0],i[1]))
    return(bodypairs)

def parse_chimera_out(chimera_data,files_in_order):
    '''does exactly as the name says...'''
    #get files shifts and axes data
    files,matrices,rmsds = [],[],[]
    lc = 0
    for i in chimera_data:
        if 'Matrix rotation and translation' in i:
            matrixdata = chimera_data[lc+1:lc+4]
            tmatrix = np.array([matrixdata[0].split(),matrixdata[1].split(),matrixdata[2].split(),[0.0,0.0,0.0,1.0]],dtype=float)
            matrices.append(tmatrix)
        if 'RMSD between' in i:
            rmsds.append(float(i.split()[6]))
        lc+=1    
    
    alldata = zip(files_in_order,matrices,rmsds)
    return(alldata)
    
def calculateCOM(pdbfile):
    '''calculate the calpha center of mass for a set of atoms'''
    x,y,z = [],[],[]
    for line in pdbfile:
        if line[13:16] == 'CA ':
            coords = line[31:56].split()
            x.append(float(coords[0]))
            y.append(float(coords[1]))
            z.append(float(coords[2]))    
    centerofmass = np.array([[np.mean(x)],[np.mean(y)],[np.mean(z)],[1]],dtype=float)
    return(centerofmass)

def divide_into_bodies(clean_chimera_data):
    '''separate the stuff to draw data into individual bodies'''
    bodies = {}
    for i in clean_chimera_data:
        bodname = i[0][0].split('_')[0]
        if bodname not in bodies:
            bodies[bodname] = [i]
        else:
            bodies[bodname].append(i)
    return(bodies)

def draw_globes(chimeraline,ovec1,ovec2,ovec3):
    '''draws the globes representing the motions (rot and trans) of the bodies'''
    vectors = [[],[]]
    first = True
    colors = ('red','blue','yellow','purple','pink','green')
    cc = 0
    ## calculate the points for the transformed COM and xyz axes
    count = 1
    for i in chimeraline:
        if cc == 6:
            cc = 0
        centerofmass =  calculateCOM(open('{0}/{1}'.format(alignedpath,i[0][0]),'r').readlines())
        vectors[0].append(centerofmass)
        ovpoint1 = np.array([[centerofmass.item(0,0)+ovec1.item(0,0)],[centerofmass.item(1,0)+ovec1.item(1,0)],[centerofmass.item(2,0)+ovec1.item(2,0)],[1]],dtype=float)
        ovpoint2 = np.array([[centerofmass.item(0,0)+ovec2.item(0,0)],[centerofmass.item(1,0)+ovec2.item(1,0)],[centerofmass.item(2,0)+ovec2.item(2,0)],[1]],dtype=float)
        ovpoint3 = np.array([[centerofmass.item(0,0)+ovec3.item(0,0)],[centerofmass.item(1,0)+ovec3.item(1,0)],[centerofmass.item(2,0)+ovec3.item(2,0)],[1]],dtype=float)
  
        xformed_COM = np.dot(i[1],centerofmass)
        xovpoint1 = np.dot(i[1],ovpoint1)
        xovpoint2 = np.dot(i[1],ovpoint2)
        xovpoint3 = np.dot(i[1],ovpoint3)
        vectors[1].append(xformed_COM)

        ## draw globes here
        if first == True:
            output.write('.comment movement {0} to {1}\n'.format(chimeraline[0][0][0],chimeraline[0][0][1]))
            output.write('.color {0}\n'.format(colors[cc]))
            output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1 \n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
            output.write('.sphere {0} {1} {2} 0.5\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0)))
            output.write('.color white\n')
            output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint1.item(0,0),ovpoint1.item(1,0),ovpoint1.item(2,0)))
            output.write('.color orange\n')
            output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint2.item(0,0),ovpoint2.item(1,0),ovpoint2.item(2,0)))
            output.write('.color green\n')
            output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint3.item(0,0),ovpoint3.item(1,0),ovpoint3.item(2,0)))
            cc+=1
        if first == False:
            output.write('{0} {1} {2} 0.1\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
        first = False
        output.write('.comment movement {0} to {1}\n'.format(chimeraline[0][0][0],chimeraline[0][0][1]))
        output.write('.color {0}\n'.format(colors[cc]))
        output.write('.sphere {0} {1} {2} 0.5\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
        output.write('.color white\n')
        output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint1.item(0,0),xovpoint1.item(1,0),xovpoint1.item(2,0)))
        output.write('.color orange\n')
        output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint2.item(0,0),xovpoint2.item(1,0),xovpoint2.item(2,0)))
        output.write('.color green\n')
        output.write('.v {0} {1} {2} {3} {4} {5}\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint3.item(0,0),xovpoint3.item(1,0),xovpoint3.item(2,0)))
        try:
            test = chimeraline[count]
            output.write('.color {0}\n'.format(colors[cc]))
            output.write('.cylinder {0} {1} {2} '.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
            count +=1
        except:
            pass
        cc+=1
        
    return(i[1])

def getpoints(pdbfile):
    aas = []
    coordlist= []
    for line in pdbfile:
        if line[13:16] == 'CA ':
            coords = line[31:56].split()
            coordlist.append((float(coords[0]),float(coords[1]),float(coords[2])))
            aas.append(line[23:26])
    return(coordlist,aas)

def calc_dist(xyz1,xyz2):
    xdif = ((xyz1[0]-xyz2[0])**2)
    ydif = ((xyz1[1]-xyz2[1])**2)
    zdif = ((xyz1[2]-xyz2[2])**2)
    distance = math.sqrt(xdif+ydif+zdif)
    return distance,xdif+ydif+zdif

def do_deviation_analysis(modpdb,refpdb):
    devaout = open('{0}/deva_{1}_ft_{2}.txt'.format(devapath,modpdb.replace('.pdb',''),refpdb.replace('.pdb','')),'w')
    devaout.write('aa,dev(Angstrom)')
    coords1,aanos1 =  getpoints(open('{0}/{1}'.format(alignedpath,modpdb),'r').readlines())
    coords2,aanos2 =  getpoints(open('{0}/{1}'.format(alignedpath,refpdb),'r').readlines())
    rmsdrunning = []
    distances = {}
    alldists = []
    for i in zip(coords1,coords2,aanos1):
        distances[i[2]],sdval=(calc_dist(i[0],i[1]))
        alldists.append(float(distances[i[2]]))
        rmsdrunning.append(sdval)
        devaout.write('\n{0}\t{1}'.format(i[2],float(distances[i[2]])))
    dmax = max(alldists)
    dmin = min(alldists)
    rmsd = math.sqrt(sum(rmsdrunning)/float(len(rmsdrunning)))
    chimeraout = ['open #0 {0}/{1};background solid white;'.format(alignedpath,refpdb)]
    for i in distances:
        try:
            color = 255 - int(255*((distances[i]-dmin)/(dmax-dmin)))
        except:
            color = 0
        chimeraout.append('color #FF{0}{0} :{1};'.format("%0.2X" % color,i))
    chimeraout.append('colorkey  0.75,0.65  0.85,0.85 labelcolor black {0} white  {1} red'.format(dmin,dmax))
    
    print('{0}\t{1}\t{2}\t{3}\t{4} / {5}'.format(round(rmsd,2),round(dmin,2),round(dmax,2),refpdb.split('_')[0],'_'.join(refpdb.split('_')[1:]).replace('.pdb',''),'_'.join(modpdb.split('_')[1:]).replace('_ft_{0}'.format(refpdb),'')))
    outfile = open('{0}/deva_{1}.cmd'.format(devapath,modpdb.replace('.pbd','')),'w')
    outfile.write(''.join(chimeraout))
    devaout.close()

def deviation_analysis(bodydic):
    deva_outs = {}
    bdkeys = list(bodydic)
    bdkeys.sort()
    for i in bodydic:
        for j in bodydic[i]:
            try:
                deva_outs[j[0][0]].append(do_deviation_analysis('{0}_ft_{1}'.format(j[0][0].replace('.pdb',''),j[0][1]),'{0}'.format(j[0][1])))
            except:
                deva_outs[j[0][0]]=[do_deviation_analysis('{0}_ft_{1}'.format(j[0][0].replace('.pdb',''),j[0][1]),'{0}'.format(j[0][1]))]


########------------------------------------------------------------------------

make_directories()

pdb_files = get_files()                             # [filename1, ... filenameN]

all_bodies = get_bodies()                           # {body name: [reidues],chain}

sub_pdbs = make_sub_pdbs(pdb_files,all_bodies)      # [(pdb file, body), ... (pdb file, body) ]

(body_dic,bkeys) = order_bodies(sub_pdbs)                   # {body name: [subpdb files in order ]}, [keys to that list in order]

body_pairs_inorder = make_body_pairs(bkeys,body_dic)     # [pair of bodies to analyze in order]

## write the chimera script
output = open('{0}/chimera_script.cmd'.format(tmppath),'w')         
for i in body_pairs_inorder:
    write_chimera_script(alignedpath,i[0],i[1]) 
output.close()   

# run the chimera script - capture the output:
runchimera = subprocess.Popen('{0} --nogui {1}/chimera_script.cmd 2>/dev/null'.format(chimerapath,tmppath), shell=True, stdout=subprocess.PIPE)
chimeraout = runchimera.stdout.read()
chimera_outfile = open('{0}/movements_raw.txt'.format(resultspath),'w')
chimera_data = []
for i in chimeraout.split('\n'):
    chimera_outfile.write('{0}\n'.format(i))
    chimera_data.append(i)
all_data = parse_chimera_out(chimera_data,body_pairs_inorder)

# divide the data into the individual bodies
boddic = divide_into_bodies(all_data)

# do the movement analysis and write the bild files
output = open('{0}/movements.bild'.format(resultspath),'w')
outdata = open('{0}/movements_data.bild'.format(resultspath),'w')
ortho1 = np.array([[10],[0],[0],[1]])
ortho2 = np.array([[0],[10],[0],[1]])
ortho3 = np.array([[0],[0],[10],[1]])
analysisdic = {}
for bod in boddic:
    (draw_globes(boddic[bod],ortho1,ortho2,ortho3))

## write the output file:
bodykeys = list(boddic)
bodykeys.sort()
for i in bodykeys:
    outdata.write('\nBody {0}\n'.format(i))
    for j in boddic[i]:
        outdata.write('\nmovement {0}\n'.format(' to '.join([x.replace('.pdb','') for x in j[0]])))
        outdata.write('{0}\n'.format(str(j[1][0:3]).replace('[',' ').replace(']',' ')))


#do the Ca deviation analysis
print('''\n---------------------\nCa deviation analysis\n---------------------\nRMSD\tmin\tmax\tbody\tpdbs''')
deviation_analysis(boddic)


### cleanup
print('\n----------------------\nCleaning up temp files\n----------------------')
subprocess.call(['rm','-r','TMP'])
print('FINISHED')
