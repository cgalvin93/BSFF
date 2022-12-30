#time python enumerate_hybrid_motifs.py clusterresults.json n_polar_clust_each_frag lig_resname nmotifres ligparams
#time python ~/desktop/bsff_scripts/enumerate_hybrid_motifs.py ~/desktop/polmotifdev/a8s_polar_contacts_clustered.json 5 a8s
#time python ~/bsff_scripts/enumerate_hybrid_motifs.py a8s_polar_contacts_clustered.json 20 a8s 3 lig.parans

import json
import sys
import os
from pyrosetta import *
init('-ignore_unrecognized_res -load_PDB_components False')
import math
import numpy as  np
import random


#first load the clustering results
clusterjson=sys.argv[1]
df=open(clusterjson,'r')
clusterdf=json.load(df)
df.close()

n_clust_each_frag=int(sys.argv[2])

ligand_resname=sys.argv[3]

nmotifres=int(sys.argv[4])

ligparams = [sys.argv[5]]
'''
#dev
df=open('a8s_polar_contacts_clustered.json','r')
clusterdf=json.load(df)
df.close()
n_clust_each_frag=5
ligand_resname='a8s'
ligparams=['/wynton/home/kortemme/cgalvin/drug_targets/a8s/Inputs/Rosetta_Inputs/a8s.params']
'''
'''
here developing the new stuff

i simply need to access the nonpolar clusters, and then convert them to cst blocks in
precisely the same way as i did in the original script,
but to somehow sabe these bloks somewhere that i can draw from when building the cst files
hoping i can leave the original loop that i do this in but add something to skip ones that contain
multiple args and to add an additional block for one of the nonpolar constraints

okay firstly which fragments do i wanna use for the nonpolar contacts and how di i want to specify them?

how about i get the top N most populated nonpolar residue clusters from across all fragments and then
convert them to cst blocks that i append to a list and then for each hbond motif i make a version
containing each of these n constraints, so if m is the number of polar clusters i am enumerating through,
i will get (n^3)*m-(z multi arg) motifs.
    i dont want to match with more than, say, 20k motifs, at least for this test run,
    and i must recall that ill need to be even more selective when matching with lucs designs,
    so that means if m is 5, which seems like it should be bare minimum, n needs to be the third root of 4k, which is 15
        given that i can expect a good number of multi arg motifs to be discarded, it feels appropriate
        to stick with using 20 clusters to enumerate (worst case 40k if none discarded)

mmmm, alternatively i can just randomly append a single np contact to each motif so i get same number,
and mauybe i also DO want to exclude polar fragments that i use for hbonds, because perhaps they
will be more likely to clash with hbond contacts
'''
#dont look at clusters of polar frags that i already reclustered...
excludepolarfrags=[]
for key in clusterdf.keys():
    excludepolarfrags.append(key)
# excludepolarfrags=[]
# for key in clusterdf.keys():
#     if key!='Fragment_2':
#         excludepolarfrags.append(key)
#okay so first let me get to the fragments and np clusters
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY']
npclusts=[] #this should hold all np cluster paths and frequencies for all fragments
for frag in os.listdir(os.path.join(ligand_resname,'Transformed_Aligned_PDBs')):
    if frag not in excludepolarfrags:
        nall=0
        clusters_path=os.path.join(ligand_resname,'Transformed_Aligned_PDBs',frag,'clustered_contacts')
        clusts=[i for i in os.listdir(clusters_path) if i[-3:]=='pdb']
        for clust in clusts:
            s=clust.split('.')[0]
            ss=s.split('_')
            restype=ss[1]
            if restype not in polar_residues:
                ninclust=int(ss[3])
                nall+=ninclust
        for clust in clusts:
           s=clust.split('.')[0]
           ss=s.split('_')
           restype=ss[1]
           if restype not in polar_residues:
               ninclust=int(ss[3])
               clustfreq=float(ninclust)/nall
               npclusts.append((os.path.join(clusters_path,clust),clustfreq))

nps=sd=sorted(npclusts, key=lambda x: x[1],reverse=True)
'''
okay yeah this appears to have gotten me the list that I need, gotta convert top N of these to
cst blocks now and store them in a list ill later draw from
'''
#functions to get distance between 2 points
def displace(p1,p2):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    return (x,y,z)
def norm(x):
    return math.sqrt(sum(i**2 for i in x))
def dist(p1,p2):
    v = displace(p1,p2)
    return norm(v)
#function to get angle between 3 points
def getAngle(p1, p2, p3):
    a=np.array(p1)
    b=np.array(p2)
    c=np.array(p3)
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)
#function to get dihedral from 4 points
def dihedral(i1,i2,i3,i4):
    p0=np.array(i1)
    p1=np.array(i2)
    p2=np.array(i3)
    p3=np.array(i4)
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)
    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
    y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)
    return np.degrees(np.arctan2(y, x))

#this function i will use just to add np cst at end, polar cst blocks are built manually
def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.35, angle_A_tolerance=15, angle_B_tolerance=15,
                                          torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                          torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                          distance_constraint_sample_number=0):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue3'


    constraint_block = [
        '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(ligand_atoms)),
        '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
        '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(residue_atoms)),
        '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, residue_resname),
        '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
            distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
        '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionA), torsion_A_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionB), torsion_B_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionAB), torsion_AB_tolerance, 100,
            torsion_constraint_sample_number)]
    return constraint_block


'''
hmmmm okay now how do i wanna define np csts for top clusts
    can score and use lowest scoring
    can simply take the first residue
i think the np contacts are sufficiently similar within each cluster that taking the first res
should be fine, so ill try that
'''
np_cst_blocks=[]
for clust,freq in nps[:25]: #####can change this from 20 to anything
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    ligand_resnum=p.total_residue()
    contact_pose=Pose()
    ligand_pose=p.residue(ligand_resnum).clone()
    res_pose=p.residue(1).clone()
    contact_pose.append_residue_by_jump(res_pose, 1)
    contact_pose.append_residue_by_jump(ligand_pose, 1)
    res_resname=contact_pose.residue(1).name()[0:3]
    lig_resname=contact_pose.residue(2).name()[0:3]
    #get x,y,z coordinates for every atom in residue and ligand
    ligand_atoms_xyz={}#atomindex=(x,y,z,index)
    residue_atoms_xyz={}
    n_residue_atoms=contact_pose.residue(1).natoms()
    n_ligand_atoms=contact_pose.residue(2).natoms()
    for k in range(1,n_ligand_atoms):
        x,y,z=contact_pose.residue(2).atom(k).xyz()
        ligand_atoms_xyz[(contact_pose.residue(2).atom_name(k)).strip()]=(x,y,z,k)
    for j in range(1,n_residue_atoms):
        x,y,z=contact_pose.residue(1).atom(j).xyz()
        residue_atoms_xyz[(contact_pose.residue(1).atom_name(j)).strip()]=(x,y,z,j)
    #find 2 atoms with shortest distance, will define atom1 for each res in constraint block
    distances=[]
    for key in ligand_atoms_xyz.keys():
        p1=ligand_atoms_xyz[key][:3]
        index1=ligand_atoms_xyz[key][3]
        for key2 in residue_atoms_xyz.keys():
            p2=residue_atoms_xyz[key2][:3]
            index2=residue_atoms_xyz[key2][3]
            d=dist(p1,p2)
            distances.append((key,key2,d,index1,index2))
    sd=sorted(distances, key=lambda x: x[2])
    #remove hydrogens
    for r in range(10):#not clear to me why but need to repeat to be thorough
        for i in sd:
            if 'H' in i[0] or 'H' in i[1]:
                sd.remove(i)
    ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[0]
    #now find base atoms for res and lig, will be atoms 2 and 3 for each
    ligatomsequence=[]
    bondedto1=[]
    for ia in range(1,n_ligand_atoms):
        atom1 = AtomID(indexlig, 2)
        if ia!=indexlig:
            atom2 = AtomID(ia, 2)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto1.append(ia)
    for x in bondedto1:
        for ia2 in range(1,n_ligand_atoms):
            atom1 = AtomID(x, 2)
            if ia2!=x and ia2!=indexlig:
                atom2 = AtomID(ia2, 2)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    if 'H' in (contact_pose.residue(2).atom_name(ia2)).strip():
                        continue
                    else:
                        ligatomsequence.append((indexlig,x,ia2))
    resatomsequence=[]
    bondedto1r=[]
    for ia in range(1,n_residue_atoms):
        atom1 = AtomID(indexres, 1)
        if ia!=indexres:
            atom2 = AtomID(ia, 1)
            if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                bondedto1r.append(ia)
    for x in bondedto1r:
        for ia2 in range(1,n_residue_atoms):
            atom1 = AtomID(x, 1)
            if ia2!=x and ia2!=indexres:
                atom2 = AtomID(ia2, 1)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    if 'H' in (contact_pose.residue(1).atom_name(ia2)).strip():
                        continue
                    else:
                        resatomsequence.append((indexres,x,ia2))
    #
    if len(ligatomsequence)>0:
        ligbase1=ligatomsequence[0][1]
        ligbase2=ligatomsequence[0][2]
        ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
        ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
    else:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(resatomsequence)>0:
        resbase1=resatomsequence[0][1]
        resbase2=resatomsequence[0][2]
        residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
        residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
    else:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    #fix oxt to O in res if present
    if residue_atom1=='OXT':
        residue_atom1='O'
    if residue_atom2=='OXT':
        residue_atom2='O'
    if residue_atom3=='OXT':
        residue_atom3='O'
    #save res and ligand atoms
    lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
    res_atoms=[residue_atom1,residue_atom2,residue_atom3]
    if len(set(res_atoms))<3:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT RES ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(set(lig_atoms))<3:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    res_atom_coords=[]
    lig_atom_coords=[]
    x,y,z=contact_pose.residue(1).atom(indexres).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(1).atom(resbase1).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(1).atom(resbase2).xyz()
    res_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(indexlig).xyz()
    lig_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(ligbase1).xyz()
    lig_atom_coords.append((x,y,z))
    x,y,z=contact_pose.residue(2).atom(ligbase2).xyz()
    lig_atom_coords.append((x,y,z))
    #okay getting angles
    #RES1 IN CONSTRAINT FILE IS GONNA BE LIGAND
    angle_A=getAngle(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
    angle_B=getAngle(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
    #finally, getting dihedrals
    torsion_A=dihedral(lig_atom_coords[2],lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
    torsion_AB=dihedral(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
    torsion_B=dihedral(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1],res_atom_coords[2])
    #get constraint block
    cstblock=generate_single_constraint_block_base(distance_AB,angle_A,angle_B,torsion_A,torsion_AB,torsion_B,
                                                   res_resname,lig_resname,lig_atoms,res_atoms,
                                                   distance_tolerance=0.35, angle_A_tolerance=15, angle_B_tolerance=15,
                                                   torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                   torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                   distance_constraint_sample_number=0)
    for wp in range(int(1000*freq)):
        np_cst_blocks.append(cstblock)


'''
now doing the old stuff but with addition of a randomly selected cst block from the np_cst_blocks list
c hcance weighted by cluster frequency->number of entries cst block has in list =freq*1000
'''

#load the list of bin values for each parameter
#
dist_ind=[i/10 for i in range(14,36,1)]
ang_ind=[i for i in range(0,190,10)]
tor_ind=[i for i in range(-180,190,10)]
#
dist_ind_pairs=[]
for i,e in enumerate(dist_ind[:-1]):
    upperlim=dist_ind[i+1]
    dist_ind_pairs.append((e,upperlim))
ang_ind_pairs=[]
for i,e in enumerate(ang_ind[:-1]):
    upperlim=ang_ind[i+1]
    ang_ind_pairs.append((e,upperlim))
tor_ind_pairs=[]
for i,e in enumerate(tor_ind[:-1]):
    upperlim=tor_ind[i+1]
    tor_ind_pairs.append((e,upperlim))




motif_clust_indices=[]
#now how can I properly enumerate the indices of the motifs that I want
if nmotifres==3:
    for x in range(n_clust_each_frag):
        for y in range(n_clust_each_frag):
            for z in range(n_clust_each_frag):
                motif_clust_indices.append((x,y,z))
elif nmotifres==4:
    for x in range(n_clust_each_frag):
        for y in range(n_clust_each_frag):
            for z in range(n_clust_each_frag):
                for j in range(n_clust_each_frag):
                    motif_clust_indices.append((x,y,z,j))

#lists of middle values of bins for each parameter
#will use this as input cst value with tolerance aligning with bin range
#ie ang bins = 10 degrees, ex val = 175, tol = +/-5
dist_ind=[i/100 for i in range(145,350,10)]
ang_ind=[i for i in range(5,180,10)]
tor_ind=[i for i in range(-175,180,10)]
residue_tag = 'residue3'
torsion_constraint_sample_number=8  #IDK WHAT IDEAL VALUE IS FOR THIS, TRY 8 AND SEE WHAT HAPPENS BUT DONT FORGET TO PLAY AROUND WITH THIS
#okay now writing motif csts
count=1
for indset in motif_clust_indices:
    cstname='hybridcst_'+str(count)+'.cst'
    of=open(cstname,'w')
    motifresnames=[]
    for frag_ind in range(len(indset)):
        contact=clusterdf[list(clusterdf.keys())[frag_ind]][indset[frag_ind]][1]
        contact_data=contact.split('_')
        res_atoms=[contact_data[0],contact_data[1],contact_data[2]]
        if 'H' in res_atoms[0]:
            resdonor=True
        else:
            resdonor=False
        lig_atoms=[contact_data[3],contact_data[4],contact_data[5]]
        res_resname=contact_data[6]
        motifresnames.append(res_resname)
        dbin=int(contact_data[7])-1
        phibin=int(contact_data[8])-1
        psibin=int(contact_data[9])-1
        try:
            chibin=int(contact_data[10])-1
        except:
            pass
        distanceAB=dist_ind[dbin]
        distance_tolerance=0.05
        distance_constraint_sample_number=0
        if resdonor==True:
            angle_A=ang_ind[psibin]
            angle_B=ang_ind[phibin]
            try:
                torsion_A=tor_ind[int(contact_data[10])-1]
                torsion_A_tolerance=5
                torsion_A_constraint_sample_number=1
            except:
                torsion_A=0
                torsion_A_tolerance=180
                torsion_A_constraint_sample_number=torsion_constraint_sample_number
            torsion_AB=0
            torsion_B=0
            angle_A_tolerance=5
            angle_B_tolerance=5
            torsion_AB_tolerance=180
            torsion_B_tolerance=180
            cstblock = [
                '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(lig_atoms)),
                '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
                '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(res_atoms)),
                '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, res_resname),
                '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                    distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
                '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_A), angle_A_tolerance, 100, 1),
                '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_B), angle_B_tolerance, 100, 1),
                '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_A), torsion_A_tolerance, 100,
                    torsion_A_constraint_sample_number),
                '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_B), torsion_B_tolerance, 100,
                    torsion_constraint_sample_number),
                '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_AB), torsion_AB_tolerance, 100,
                    torsion_constraint_sample_number)]
        elif resdonor==False:
            angle_A=ang_ind[phibin]
            angle_B=ang_ind[psibin]
            try:
                torsion_B=tor_ind[int(contact_data[10])-1]
                torsion_B_tolerance=5
                torsion_B_constraint_sample_number=1
            except:
                torsion_B=0
                torsion_B_tolerance=180
                torsion_B_constraint_sample_number=torsion_constraint_sample_number
            torsion_AB=0
            torsion_A=0
            angle_A_tolerance=5
            angle_B_tolerance=5
            torsion_AB_tolerance=180
            torsion_A_tolerance=180
            cstblock = [
                '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(lig_atoms)),
                '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
                '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(res_atoms)),
                '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, res_resname),
                '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                    distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
                '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_A), angle_A_tolerance, 100, 1),
                '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_B), angle_B_tolerance, 100, 1),
                '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_A), torsion_A_tolerance, 100,
                    torsion_constraint_sample_number),
                '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_B), torsion_B_tolerance, 100,
                    torsion_B_constraint_sample_number),
                '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_AB), torsion_AB_tolerance, 100,
                    torsion_constraint_sample_number)]
        of.write('CST::BEGIN\n')
        of.write('\n'.join(cstblock))
        of.write('\n')
        of.write('CST::END\n')
        of.write('\n')
    narg=0
    for tn in motifresnames:
            if tn=='ARG':
                narg+=1
    if narg>1:
        of.close()
        os.system('rm '+cstname)
        continue
    elif len(set(motifresnames))<2:
        of.close()
        os.system('rm '+cstname)
        continue
    else:
        rnnp=random.randrange(0,len(np_cst_blocks))
        npcstnlock=np_cst_blocks[rnnp]
        of.write('CST::BEGIN\n')
        of.write('\n'.join(npcstnlock))
        of.write('\n')
        of.write('CST::END\n')
        of.write('\n')
        of.close()
        count+=1
