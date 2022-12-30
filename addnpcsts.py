#
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random
##################################################
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs/a8s_Fragment_4_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params']
# polar_cst_paths='/wynton/home/kortemme/cgalvin/r6matching/a8s_4p'
motifsoutdir='npmotifpdbs'
# polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/a8s/newpm_test_polar_csts_rand/three_res'
polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/a8s/new_t3_polar_csts_rand/three_res'
##################################################
##################################################
# fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_4_residue_contacts','/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_5_residue_contacts','/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_6_residue_contacts']
# ligparams=['/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params']
# # polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e_3res_polar_polar_csts'
# # polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e4p_polar_csts'
# # polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e_polar_csts_rand/three_res'
# polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e_polar_csts_rand/four_res'
# motifsoutdir='npmotifpdbs'
##################################################
os.makedirs(motifsoutdir,exist_ok=True)

clust_paths=[]
for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))

clust_paths_pops=[]
for clustpath in clust_paths:
    pop=int(clustpath.split('/')[-1].split('.')[0].split('_')[3])
    clust_paths_pops.append((clustpath,pop))

nps=sorted(clust_paths_pops, reverse=True, key=lambda nmem: nmem[1])


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
                                          distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
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


sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)

# min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
# mm4060 = pyrosetta.rosetta.core.kinematics.MoveMap()
# mm4060.set_bb(False)
# mm4060.set_chi(False)
# mm4060.set_jump(1,True)
# min_mover.movemap(mm4060)
# min_mover.score_function(sf)
# min_mover.max_iter(1)


np_cst_blocks=[]
accepted_clusts=[]
ac=0
pdboutpaths=[]
clustnames_blocks=[]
for clust,freq in nps[:25]:
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    ligand_resnum=p.total_residue()
    clust_cont_scores=[]
    for indd in range(1,p.total_residue()):
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(indd).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        res_resname=contact_pose.residue(1).name()[0:3]
        lig_resname=contact_pose.residue(2).name()[0:3]
        if res_resname==lig_resname:
            print('uh oh')
            continue
        else:
            pass
        # min_mover.apply(contact_pose)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((indd,score_this_one))
    sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[1])
    if sccs[0][1]<-2.5:
        accepted_clusts.append((clust,freq,sccs[0][1]))
        ac+=1
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(sccs[0][0]).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        pdboutpath=os.path.join(motifsoutdir,'npmotif_'+str(ac)+'.pdb')
        contact_pose.dump_pdb(pdboutpath)
        pdboutpaths.append(pdboutpath)
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
        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
            try:
                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[1]
                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                    try:
                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[2]
                        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                            try:
                                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[3]
                                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                                    try:
                                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[4]
                                    except:
                                        pass
                            except:
                                pass
                    except:
                        pass
            except:
                pass
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
            continue
        if len(resatomsequence)>0:
            resbase1=resatomsequence[0][1]
            resbase2=resatomsequence[0][2]
            residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
            residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
        else:
            continue
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
            continue
        if len(set(lig_atoms))<3:
            continue
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
                                                       distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                                                       torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                       torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                       distance_constraint_sample_number=0)
        for wp in range(freq):
            np_cst_blocks.append(cstblock)
        clustnames_blocks.append((pdboutpath,cstblock))



# accepted_clusts


os.mkdir('1np')
os.chdir('1np')
#to add 1 np cst to csts
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(1):          #this is how many unique hybrid cst i want to make per polar cst
            rblock=random.choice(np_cst_blocks)
            if rblock not in bs:
                bs.append(rblock)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblock))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                for line in lines:
                    of.write(line)
                of.close()
                print(c)
                c+=1
            else:
                rblock=random.choice(np_cst_blocks)
                if rblock not in bs:
                    bs.append(rblock)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    for line in lines:
                        of.write(line)
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass


os.chdir('..')

import os
l=[i for i in os.listdir() if i[-3:]=='cst']
len(l)
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)<40:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
len(l)

#to make 2 residue motifs
clust_cont_scores=[]
for cind,clust in enumerate(pdboutpaths):
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    for clust2 in pdboutpaths[cind+1:]:
        p2 = Pose()
        generate_nonstandard_residue_set(p2,ligparams)
        pose_from_file(p2, clust2)
        #
        contact_pose=Pose()
        res_pose=p.residue(1).clone()
        res_pose2=p2.residue(1).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(res_pose2, 1)
        ligand_pose=p.residue(2).clone()
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((clust,clust2,score_this_one))


sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[2])
cstblockss=[]
for x in sccs:
    if x[2]<-3.0:
        cstblocks=[]
        for y in clustnames_blocks:
            if x[0]==y[0]:
                cstblocks.append(y[1])
            elif x[1]==y[0]:
                cstblocks.append(y[1])
        cstblockss.append((cstblocks[0],cstblocks[1]))



os.mkdir('2np')
os.chdir('2np')
#should do this in directory where you want them output
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(1):
            rblocks=random.choice(cstblockss)
            if rblocks not in bs:
                bs.append(rblocks)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                for line in lines:
                    of.write(line)
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[0]))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[1]))
                of.write('\n')
                of.write('CST::END\n')
                of.close()
                print(c)
                c+=1
            else:
                rblocks=random.choice(cstblockss)
                if rblocks not in bs:
                    bs.append(rblocks)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    for line in lines:
                        of.write(line)
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[0]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[1]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass
os.chdir('..')




'''

































now to do 4 rewidue polar

for a8s, wanna copy fragment 3 to a fragment 5
'''
import json
#first load the clustering results
clusterjson='a8s_polar_contacts_clustered_liberal.json'
df=open(clusterjson,'r')
clusterdf=json.load(df)
df.close()

d2={}
for key in clusterdf.keys():
    d2[key]=clusterdf[key]
#
d2['Fragment_5']=clusterdf['Fragment_3']
#
out_filejson = open('a8s_polar_clusters_4res.json','w')
json.dump(d2, out_filejson, indent = 6)
out_filejson.close()

'''
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_clusters_4res.json 20 a8s 4 a8s4p 250 5
compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_clusters_4res.json 40 a8s 4 a8s4p 250 5
qstat -j "$JOB_ID"
'>epm4_$compound_dir.sh
chmod ugo+x epm4_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output epm4_$compound_dir.sh


'''






















'''
38E
i need the nonpolar clustering to finish in order to generate the hybrids
and i need to do the 4 res polar motifs (and then their hybrids)
    also need to do a8s 4 res hybrids
'''













#
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random
##################################################
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_4_residue_contacts','/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_5_residue_contacts','/wynton/home/kortemme/cgalvin/round6/38e/38efragfbs/38e_Fragment_6_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params']
# polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e_3res_polar_polar_csts'
polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/38e/38e4p_polar_csts'

motifsoutdir='npmotifpdbs'
##################################################
os.makedirs(motifsoutdir,exist_ok=True)

clust_paths=[]
for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))

clust_paths_pops=[]
for clustpath in clust_paths:
    pop=int(clustpath.split('/')[-1].split('.')[0].split('_')[3])
    clust_paths_pops.append((clustpath,pop))

nps=sorted(clust_paths_pops, reverse=True, key=lambda nmem: nmem[1])


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


sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)

# min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
# mm4060 = pyrosetta.rosetta.core.kinematics.MoveMap()
# mm4060.set_bb(False)
# mm4060.set_chi(False)
# mm4060.set_jump(1,True)
# min_mover.movemap(mm4060)
# min_mover.score_function(sf)
# min_mover.max_iter(1)


np_cst_blocks=[]
accepted_clusts=[]
ac=0
pdboutpaths=[]
clustnames_blocks=[]
for clust,freq in nps[:25]:
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    ligand_resnum=p.total_residue()
    clust_cont_scores=[]
    for indd in range(1,p.total_residue()):
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(indd).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        res_resname=contact_pose.residue(1).name()[0:3]
        lig_resname=contact_pose.residue(2).name()[0:3]
        if res_resname==lig_resname:
            print('uh oh')
            continue
        else:
            pass
        # min_mover.apply(contact_pose)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((indd,score_this_one))
    sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[1])
    if sccs[0][1]<-1.5:
        accepted_clusts.append((clust,freq,sccs[0][1]))
        ac+=1
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(sccs[0][0]).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        pdboutpath=os.path.join(motifsoutdir,'npmotif_'+str(ac)+'.pdb')
        contact_pose.dump_pdb(pdboutpath)
        pdboutpaths.append(pdboutpath)
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
        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
            try:
                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[1]
                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                    try:
                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[2]
                        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                            try:
                                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[3]
                                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                                    try:
                                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[4]
                                    except:
                                        pass
                            except:
                                pass
                    except:
                        pass
            except:
                pass
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
            continue
        if len(resatomsequence)>0:
            resbase1=resatomsequence[0][1]
            resbase2=resatomsequence[0][2]
            residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
            residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
        else:
            continue
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
            continue
        if len(set(lig_atoms))<3:
            continue
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
        for wp in range(freq):
            np_cst_blocks.append(cstblock)
        clustnames_blocks.append((pdboutpath,cstblock))





os.mkdir('1np')
os.chdir('1np')
#to add 1 np cst to csts
#should do this in directory where you want them output
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(8):
            rblock=random.choice(np_cst_blocks)
            if rblock not in bs:
                bs.append(rblock)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                for line in lines:
                    of.write(line)
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblock))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                of.close()
                print(c)
                c+=1
            else:
                rblock=random.choice(np_cst_blocks)
                if rblock not in bs:
                    bs.append(rblock)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    for line in lines:
                        of.write(line)
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass


os.chdir('..')


#to make 2 residue motifs
clust_cont_scores=[]
for cind,clust in enumerate(pdboutpaths):
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    for clust2 in pdboutpaths[cind+1:]:
        p2 = Pose()
        generate_nonstandard_residue_set(p2,ligparams)
        pose_from_file(p2, clust2)
        #
        contact_pose=Pose()
        res_pose=p.residue(1).clone()
        res_pose2=p2.residue(1).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(res_pose2, 1)
        ligand_pose=p.residue(2).clone()
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((clust,clust2,score_this_one))


sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[2])
cstblockss=[]
for x in sccs:
    if x[2]<-3.0:
        cstblocks=[]
        for y in clustnames_blocks:
            if x[0]==y[0]:
                cstblocks.append(y[1])
            elif x[1]==y[0]:
                cstblocks.append(y[1])
        cstblockss.append((cstblocks[0],cstblocks[1]))



os.mkdir('2np')
os.chdir('2np')
#should do this in directory where you want them output
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(8):
            rblocks=random.choice(cstblockss)
            if rblocks not in bs:
                bs.append(rblocks)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                for line in lines:
                    of.write(line)
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[0]))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[1]))
                of.write('\n')
                of.write('CST::END\n')
                of.close()
                print(c)
                c+=1
            else:
                rblocks=random.choice(cstblockss)
                if rblocks not in bs:
                    bs.append(rblocks)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    for line in lines:
                        of.write(line)
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[0]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[1]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass

'''




now to do 4 rewidue polar

for 38e, wanna copy fragment 1 as fragment 7
'''
import json
#first load the clustering results
clusterjson='38e_polar_contacts_clustered_liberal.json'
df=open(clusterjson,'r')
clusterdf=json.load(df)
df.close()

d2={}
for key in clusterdf.keys():
    d2[key]=clusterdf[key]
#
d2['Fragment_7']=clusterdf['Fragment_1']
#
out_filejson = open('38e_polar_clusters_4res.json','w')
json.dump(d2, out_filejson, indent = 6)
out_filejson.close()

'''
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_clusters_4res.json 30 38e 4 38e4p 250 5
compound_dir='38e'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_clusters_4res.json 30 38e 4 38e4p 250 5
qstat -j "$JOB_ID"
'>epm4_$compound_dir.sh
chmod ugo+x epm4_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output epm4_$compound_dir.sh


'''
