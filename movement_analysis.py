#!/usr/bin/env python

# quantify the movements of the POTRAs using Chimera mmaker and fitmap

#
#TO DO:
# add center of mass normalization
# add sanity check.. calculate Ca RMSD of fit pdb vs actual - compare to mm RMSD - should be close... flag if too high

import os
import sys
import math
import subprocess
import numpy as np

vers = 0.3
chimerapath='/fbs/emsoftware2/LINUX/fbsmi/Chimera-1.11.2-linux/bin/chimera'
alignedpath = '{0}/aligned'.format(os.getcwd())
resultspath = '{0}/Results'.format(os.getcwd())
tmppath = '{0}/TMP'.format(os.getcwd())

if len(sys.argv) < 3:
    sys.exit('''
USAGE: calculate_POTRA_movements.py <body definition file> <models file>
---
body definition file = four columns text file
body_name       start_AA        end_AA      chain

models file -  text file list of models one per line -  models must be in order of the motion    

---''')


def make_directories():
    if os.path.isdir(alignedpath) == False:
        subprocess.call(['mkdir',alignedpath])
    if os.path.isdir(resultspath) == False:
        subprocess.call(['mkdir',resultspath])
    if os.path.isdir(tmppath) == False:
        subprocess.call(['mkdir',tmppath])

def get_files():
    print('----- files to operate on -----')
    pdbfiles = open(sys.argv[2],'r').readlines()
    pdblist = []
    for i in pdbfiles:
        pdblist.append(i.replace('\n',''))
        print(i.replace('\n',''))
    return(pdblist)

def get_bodies():
    allbodies = []
    bodydeffile = open(sys.argv[1],'r').readlines()
    print('----- body definitions ------\nname    start   end     chain')
    for i in bodydeffile:
        line = i.split()
        print('{0}\t{1}\t{2}\t{3}'.format(line[0],line[1],line[2],line[3]))
        allbodies.append((line[0],range(int(line[1]),int(line[2])),line[3]))
    print('-----------------------------')
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
    print('making sub body pdb: {0}'.format(filename))
    return(filename)

def write_chimera_script(path,expPDB,refPDB):
    print('writing chimera script in: TMP/chimera_script.cmd')
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
    
    print('''----------------------\nbody models (in order)\n----------------------''')
    for i in bkeys:
        print i,bodies[i]
    return(bodies,bkeys)

def make_body_pairs(bodkeys,bodydic):
    bodypairs = []
    for i in bodkeys:
        for j in range(len(bodydic[i])):
            try:
                bodypairs.append((bodydic[i].pop(0),bodydic[i][0]))
            except:
                pass
            
    print('''------------------------------------------------------------------\nsets for comparison (mobile pdb, reference pdb, reference density)\n------------------------------------------------------------------''')
    for i in bodypairs:
        print i
    return(bodypairs)

def parse_chimera_out(chimera_data,files_in_order):
    '''does exactly as the name says...'''
    #get files shifts and axes data
    files,matrices,rmsds = [],[],[]
    lc = 0
    for i in chimera_data:
        #print i
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

def RMtoEuler(R) : 
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
 
    return([x, y, z])

def draw_globes(chimeraline,ovec1,ovec2,ovec3):
    '''draws the globes representing the motions (rot and trans) of the bodies'''
    vectors = [[],[]]
    eulerlist = []
    xyzrotlist = []
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

        #  calculate the euler angles
        eulerlist.append(RMtoEuler(i[1]))


        ## draw globes here
        if first == True:
            output.write('.color {0}\n'.format(colors[cc]))
            output.write('.v {0} {1} {2} {3} {4} {5} \n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
            output.write('.sphere {0} {1} {2} 0.15\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0)))
            output.write('.color black\n')
            output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint1.item(0,0),ovpoint1.item(1,0),ovpoint1.item(2,0)))
            output.write('.color orange\n')
            output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint2.item(0,0),ovpoint2.item(1,0),ovpoint2.item(2,0)))
            output.write('.color green\n')
            output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(centerofmass.item(0,0),centerofmass.item(1,0),centerofmass.item(2,0),ovpoint3.item(0,0),ovpoint3.item(1,0),ovpoint3.item(2,0)))
            cc+=1
        if first == False:
            output.write('{} {} {}\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
        first = False
        output.write('.color {0}\n'.format(colors[cc]))
        output.write('.sphere {0} {1} {2} 0.3\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
        output.write('.color white\n')
        output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint1.item(0,0),xovpoint1.item(1,0),xovpoint1.item(2,0)))
        output.write('.color orange\n')
        output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint2.item(0,0),xovpoint2.item(1,0),xovpoint2.item(2,0)))
        output.write('.color green\n')
        output.write('.cylinder {0} {1} {2} {3} {4} {5} 0.1\n'.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0),xovpoint3.item(0,0),xovpoint3.item(1,0),xovpoint3.item(2,0)))
        try:
            test = chimeraline[count]
            output.write('.color {0}\n'.format(colors[cc]))
            output.write('.v {0} {1} {2} '.format(xformed_COM.item(0,0),xformed_COM.item(1,0),xformed_COM.item(2,0)))
            count +=1
        except:
            pass
        cc+=1
        
    return((vectors,eulerlist))

def write_movements_output(adkeys,analysisdic):
    movementsout = open('{0}/movements_output.txt'.format(resultspath),'w')
    movementsout.write('body, movement, distance, directionx, directiony, directionz, rotx, roty, rotz')
    for i in adkeys:
        print('''----------\n{0}\n----------'''.format(i))
        mcount = 1
        for group in analysisdic[i]:
            for pair in zip(group[0][0],group[0][1],group[1]):
                dvec= pair[0]-pair[1]
                print ('movement {0}\tdistance (A)    direction'.format(mcount))
                x2 = (pair[0].item(0,0)-pair[1].item(0,0))**2
                y2 = (pair[0].item(1,0)-pair[1].item(1,0))**2
                z2 = (pair[0].item(2,0)-pair[1].item(2,0))**2
                print('          \t{0} ({1}, {2}, {3})'.format(np.round(math.sqrt(x2+y2+z2),2),np.round(dvec.item(0,0),2),np.round(dvec.item(1,0),2),np.round(dvec.item(2,0),2)))
                print 'rotation angles\tx\ty\tz'
                print('            \t{0}\t{1}\t{2}'.format(round(np.degrees(pair[2][0]),2),round(np.degrees(pair[2][1]),2),round(np.degrees(pair[2][2]),2)))
                movementsout.write('\n{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t'.format(i,mcount,np.round(math.sqrt(x2+y2+z2),2),np.round(dvec.item(0,0),2),np.round(dvec.item(1,0),2),np.round(dvec.item(2,0),2),round(np.degrees(pair[2][0]),2),round(np.degrees(pair[2][1]),2),round(np.degrees(pair[2][2]),2)))
                mcount+=1
    
    movementsout.close()

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
    devaout = open('{0}/deva_{1}_ft_{2}.txt'.format(resultspath,modpdb.replace('.pdb',''),refpdb.replace('.pdb','')),'w')
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
    
    print('Ca deviation (RMSD,min,max)   {0}\t{1}\t{2}'.format(round(rmsd,2),round(dmin,2),round(dmax,2)))
    #print(chimeraout)
    outfile = open('{0}/deva_{1}.cmd'.format(resultspath,modpdb.replace('.pbd','')),'w')
    outfile.write(''.join(chimeraout))
    devaout.close()
def deviation_analysis(bodydic):
    deva_outs = {}
    for i in bodydic:
        for j in bodydic[i]:
            print('Ca deviations between',j[0])
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
print('running chimera')
runchimera = subprocess.Popen('{0} --nogui {1}/chimera_script.cmd'.format(chimerapath,tmppath), shell=True, stdout=subprocess.PIPE)
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
ortho1 = np.array([[10],[0],[0],[1]])
ortho2 = np.array([[0],[10],[0],[1]])
ortho3 = np.array([[0],[0],[10],[1]])
analysisdic = {}
for bod in boddic:
    try:
        analysisdic[bod].append(draw_globes(boddic[bod],ortho1,ortho2,ortho3))
    except:
        analysisdic[bod] = [draw_globes(boddic[bod],ortho1,ortho2,ortho3)]
adkeys = analysisdic.keys()
adkeys.sort()
write_movements_output(adkeys,analysisdic)

#do the Ca deviation analysis
print('''---------------------\nCa deviation analysis\n---------------------\n''')
deviation_analysis(boddic)


### cleanup
print('clean up temp files')
subprocess.call(['rm','-r','TMP'])
print('FINISHED')

#########TESTING

#print pdb_files
#print all_bodies
#print sub_pdbs
#print body_dic
#print bkeys
#print body_pairs_inorder
#print(chimeraoutput)
#print all_data
