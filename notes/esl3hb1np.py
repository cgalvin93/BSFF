
######
import random
import os
#input ligand name and number of unique polar csts you want
n_desired_motifs=120
lig_resname='esl'
allowed_problematic_motifs=20
#new directory where motifs will be output
polar_cst_paths='1np3hb'
# polar_cst_paths='hb5'


#for each polar fragment, define whether acceptor is sp2 hybridized or not
#and what are the 3 constrained atoms
relevant_atoms={'Fragment_1':['n','O3','C18','C17'],
                'Fragment_2':['n','O2','C8','C7'],
                'Fragment_3':['n','O1','C5','C1']}

# #to put 2 hb on carbonyl
# relevant_atoms={'Fragment_1':['yes','O5','C23','C22'],
#                 'Fragment_2':['n','O1','C1','C3'],
#                 'Fragment_3':['n','O2','C6','C8'],
#                 'Fragment_4':['n','O3','C19','C18'],
#                 'Fragment_5':['yes','O5','C23','C22']}

#define wheteher or not vcharged residues are allowed to be used with frag
charged_allowed={'Fragment_1':'yes',
                'Fragment_2':'yes',
                'Fragment_3':'yes'}

#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
q_polar_residues=['LYS','ARG','SER','THR','GLN','ASN','TYR','TRP','GLY'] #no his
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS'] #no his
problematic_residues=['LYS','ARG','ASP','GLU','GLN','ASN','HIS']#tyr,ser,thr
#3 to 1 letter code dict
one_lcd={'LYS':'K',
         'ARG':'R',
         'SER':'S',
         'THR':'T',
         'GLN':'Q',
         'ASN':'N',
         'HIS':'H',
         'TYR':'Y',
         'TRP':'W',
         'GLY':'G'}







#pass as residue type either Hpol for all sc or HNbb for bb gly
def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.15, angle_A_tolerance=15, angle_B_tolerance=15,
                                          torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                          torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                          distance_constraint_sample_number=0):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue1'


    constraint_block = [
        '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(ligand_atoms)),
        '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
        '  TEMPLATE::   ATOM_MAP: 2 atom_type: {}'.format(residue_atoms),
        '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, one_lcd[residue_resname]),
        '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
            distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
        '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  180.00  {3:3}'.format(
            float(torsionA), torsion_A_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionB), torsion_B_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionAB), torsion_AB_tolerance, 100,
            torsion_constraint_sample_number)]
    return constraint_block




#for each fragment, create a cst block
motifs=[]
nproblematic=0
ntries=0
while len(motifs)<n_desired_motifs:
    problem_res=0
    ngly=0
    current_motif=[]
    residues_used=[]
    for frag in relevant_atoms.keys():
        sp2_acc=relevant_atoms[frag][0]
        lig_atoms=relevant_atoms[frag][1:]
        if charged_allowed[frag]=='yes':
            res_resname=random.choice(q_polar_residues)
        else:
            res_resname=random.choice(polar_residues)
        residues_used.append(res_resname)
        if res_resname=='GLY':
            res_atoms='HNbb,'
            ngly+=1
        else:
            res_atoms='Hpol,'
        if res_resname in problematic_residues:
            problem_res+=1
        if sp2_acc=='yes':
            cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                           residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                           distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=20, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                           torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                           distance_constraint_sample_number=3)
        else:
            cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=180,torsionA=180,torsionAB=0,torsionB=0,
                                                           residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                           distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                           torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                           distance_constraint_sample_number=5)
        current_motif.append(cstblock)
    if problem_res==1:
        if nproblematic<=allowed_problematic_motifs:
            if ngly<=1:
                if len(list(set(residues_used)))>=len(list(relevant_atoms.keys()))-1:
                    if current_motif not in motifs:
                        nproblematic+=1
                        motifs.append(current_motif)
    elif problem_res==0:
        if ngly<=1:
            if len(list(set(residues_used)))>=len(list(relevant_atoms.keys()))-1:
                if current_motif not in motifs:
                    motifs.append(current_motif)
    ntries+=1
    if ntries>=10000:
        break



#output the csts
os.makedirs(polar_cst_paths,exist_ok=True)
c=1
for cstf in motifs:
    cstfn='cst_'+str(c)+'.cst'
    of=open(os.path.join(polar_cst_paths,cstfn),'w')
    for block in cstf:
        of.write('CST::BEGIN\n')
        of.write('\n'.join(block))
        of.write('\n')
        of.write('CST::END\n')
    of.close()
    c+=1




















'''
ADDING NP CSTS
'''


#
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random
#####
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep,fa_elec#, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)
min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
mm4060 = pyrosetta.rosetta.core.kinematics.MoveMap()
mm4060.set_bb(False)
mm4060.set_chi(1,True)
# mm4060.set_chi(1,True)
# mm4060.set_jump(1,True)
min_mover.movemap(mm4060)
min_mover.score_function(sf)
# min_mover.max_iter(1)
min_mover.tolerance(1e-6)
##################################################
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_1_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params']
motifsoutdir='npmotifpdbs'
polar_cst_paths='/wynton/home/kortemme/cgalvin/esl/1np3hb'
############
single_contact_threshold=-2.0
double_contact_score_threshold=-5.0
triple_contact_score_threshold=-9.0
##################################################
##################################################
##################################################
#we make sure we are in the right working directory
if not os.getcwd()==polar_cst_paths:
    os.chdir(polar_cst_paths)
#make the nop motif pdb output directory
os.makedirs(motifsoutdir,exist_ok=True)




#check the polar csts to make sure theres no funky ones
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)<40:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))


#get the paths of the residue fuzzball pdbs to open
clust_paths=[]
for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))
#get the populations
clust_paths_pops=[]
for clustpath in clust_paths:
    pop=int(clustpath.split('/')[-1].split('.')[0].split('_')[3])
    clust_paths_pops.append((clustpath,pop))

#sort by population across all fragments
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
                                          distance_constraint_sample_number=2):
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






np_cst_blocks=[]
accepted_clusts=[]
ac=0
pdboutpaths=[]
clustnames_blocks=[]
for clust,freq in nps[:20]:         ####################################################################
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
    if sccs[0][1]<=single_contact_threshold:
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
                                                       distance_constraint_sample_number=2)
        for wp in range(freq):
            np_cst_blocks.append(cstblock)
        clustnames_blocks.append((pdboutpath,cstblock))



# accepted_clusts
'''
 minimize_motif = rosetta.protocols.minimization_packing.MinMover()
minimize_motif.movemap(motif_movemap)
minimize_motif.score_function(sfxn)
minimize_motif.min_type('lbfgs_armijo')
minimize_motif.tolerance(1e-6)
minimize_motif.apply(match_pose_clone)
'''

if os.path.exists('1np'):
    os.system('rm -r 1np')
    os.mkdir('1np')
else:
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
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock))
                    of.write('\n')
                    of.write('CST::END\n')
                    for line in lines:
                        of.write(line)
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass


l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)<40:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))

os.chdir('..')


'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/npmotifpdbs ~/desktop/npmotifs
'''



































'''
                            MATCHING
'''
###########################################################################
###########################################################################
###########################################################################
###########################################################################
import os
lig_name='esl'
paramspath='/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params'
allmotifsdir=os.getcwd()
#
shellfile_suffix='_esl1n3h.sh'
#
scaffold_path='/wynton/home/kortemme/cgalvin/ntf2_round3_best'
###########################################################################
############################################################################
############################################################################
############################################################################
#check that all scaffold pdb/pos paths exist
scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='pdb']
for i in scaffolds:
    p=os.path.join(scaffold_path,i.strip('.pdb')+'.pos')
    p2=os.path.join(scaffold_path,i)
    if not os.path.exists(p):
        print(i)
        scaffolds.remove(i)
    if not os.path.exists(p2):
        print(i)
        scaffolds.remove(i)
#
motifs=[i for i in os.listdir(allmotifsdir) if i[-4:]=='.cst']
paths=[os.path.join(allmotifsdir,i) for i in motifs]
motif_paths={}
motif_paths[lig_name]=paths
#
os.makedirs('cluster_out',exist_ok=True)
for key in list(motif_paths.keys()):
    lig_name=str(key)
    count=1
    idx=[j for j in range(0,len(motif_paths[key]),50)]  #how many csts per job
    if len(motif_paths[key]) not in idx:
        idx.append(len(motif_paths[key]))
    for ei, ix in enumerate(idx[:-1]):
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        of.write('\n')
        of.write('tasks=(0\n')
        for scaff in scaffolds[:-1]:
            of.write('       '+os.path.join(scaffold_path,scaff.strip('.pdb'))+'\n')
        of.write('       '+os.path.join(scaffold_path,scaffolds[-1].strip('.pdb'))+')')
        of.write('\n')
        for cst_path in motif_paths[key][ix:idx[ei+1]]:
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', '${tasks[$SGE_TASK_ID]}.pdb',  ###########
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', '${tasks[$SGE_TASK_ID]}.pos',   #########
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match::lig_name', lig_name,
                           '-load_PDB_components', 'false']
            cmd=' '.join(matcher_cmd)
            of.write(cmd+'\n')
        print(count)
        count+=1
        of.write('\nqstat -j "$JOB_ID"')
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)
print('n motifs')
print(len(paths))
print('n scaffolds')
print(str(len(scaffolds)))
'''
n motifs
230
n scaffolds
515

~/main/source/bin/match.linuxgccrelease -s /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_278534_design_2_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001.pdb -extra_res_fa /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params -match:geometric_constraint_file /wynton/home/kortemme/cgalvin/esl/hb3/2np/hybrid_25.cst -match::scaffold_active_site_residues /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_278534_design_2_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001.pos -match:output_format PDB -match:match_grouper SameSequenceGrouper -match:consolidate_matches -match:output_matches_per_group 1 -use_input_sc -ex1 -ex2 -extrachi_cutoff 0 -enumerate_ligand_rotamers false -match::lig_name esl -load_PDB_components false

File: src/protocols/match/MatcherTask.cc:1188
[ ERROR ] UtilityExitException
ERROR: No allowed residue types seen for downstream residue for constraint 3 block 1

FORGOT TO ADD 3 ATOMS FOR LIGAND IN HBOND CSTS

~/main/source/bin/match.linuxgccrelease -s /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_209757_design_9_unrelaxed_model_4_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001.pdb -extra_res_fa /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params -match:geometric_constraint_file /wynton/home/kortemme/cgalvin/esl/hb3/2np/hybrid_25.cst -match::scaffold_active_site_residues /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_209757_design_9_unrelaxed_model_4_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001.pos -match:output_format PDB -match:match_grouper SameSequenceGrouper -match:consolidate_matches -match:output_matches_per_group 1 -use_input_sc -ex1 -ex2 -extrachi_cutoff 0 -enumerate_ligand_rotamers false -match::lig_name esl -load_PDB_components false


/wynton/home/kortemme/cgalvin/esl/hb3/2np
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('_1_a8s_3p1np.sh')
#
for j in jobss:
    cmd='qsub -cwd -t 1-515 -l mem_free=20G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''
MATCHING ONLY TAKES ABOUT 1-3 HRS WITH 25 CSTS PER MATCH HERE
they also are taking less than 1g mem so i can definitely specify way less lmao

6-x hrs with occ now

15 hrs 50 ea job 1np3hb, maxvmem still no higher than 3 gb

RAN OUT OF MEMORY DURING MATCHING, I THINK PROBABLY BECAUSE OF CORE DUMPS
beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

'''

































'''genpot'''

#delete empty matches
import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(l))
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)<40:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(l))
'''
132947
98245
'''
#ANALYSIS OF LIGAND BURIAL IN MATCHES
import os
########
prm1='/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params'
sfname='matchanl.sh'
scorefile_name='matchanl.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
#######
count=1
idx=[j for j in range(0,len(matches),100)]  #how many matches per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
             '-s',matchpath,         ##
             '-parser:protocol','~/BSFF/tools/match_analysis_standard.xml',
             '-parser:view','-run:preserve_header',
             '-in:file::extra_res_fa',prm1,
             '-load_PDB_components','False',
             '-scorefile_format','json',
             '-out:file:score_only',scorefile_name]
        cmdd=' '.join(cmd)
        of.write(cmdd+'\n')
    print(count)
    count+=1
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
'''
~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s UM_1_N13S100T55Y42_1_relaxed_relaxed_5tpj_371721_design_5_unrelaxed_model_3_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_cst_4_1.pdb -parser:protocol ~/BSFF/tools/match_analysis_standard.xml -parser:view -run:preserve_header -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params -load_PDB_components False -scorefile_format json -out:file:score_only matchanl.json

'''
import os
sfname='matchanl.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output '+j
    os.system(cmd)
'''
2g is good enough
can lower to like 100 a job or something

NOW FILTER THE WELL BURIED ONES

scp cgalvin@log2.wynton.ucsf.edu:esl/hb3/1np/UM_1_F55T109T86Y93_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_51_1.pdb ~/desktop/ex.pdb
'''



#analysis of scores from json file
sfname='matchanl.json'
import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==2:
        l1=line[:ids[0]+1]
        l2=line[ids[0]+1:]
        scores.append(l1)
        scores.append(l2)
terms=list(scores[0].keys())
print(len(bad))

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

f1=return_filtered(scores,'dsasa','>',0.8)
f2=return_filtered(scores,'ligsasasasa','<',50.0)
print(len(scores))
print(len(f1))
print(len(f2))
'''
96520
52862
40863
'''
import os
filtered_strc=[]
for d in f1:
    filtered_strc.append(d['decoy'])
os.makedirs('buried',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i[:-5]+'.pdb buried/'+i+'.pdb')




























'''
genpot params generation

i should make this paralellizable

mkdir genpot
mkdir cluster_output
'''
import os
matches=[i for i in os.listdir() if i[-3:]=='pdb']
count=1
sfname='ggp.sh'
idx=[j for j in range(0,len(matches),100)]  #how many matches per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    of.write('source ~/anaconda3/etc/profile.d/conda.sh')
    of.write('\n')
    of.write('conda activate pyr37')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        cmd=['time', 'python3', '~/BSFF/tools/generate_genpot_params.py', matchpath, 'esl', '0']
        cmdd=' '.join(cmd)
        of.write(cmdd+'\n')
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    print(count)
    count+=1
'''
qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output _24ggp.sh

python ~/BSFF/tools/generate_genpot_params.py UM_3_L81F58T61Q56T38_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_68_1.pdb esl 0
'''
import os
sfname='ggp'
jobss=[i for i in os.listdir() if sfname in i]
# jobss.remove('_24ggp.sh')
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output '+j
    os.system(cmd)

'''
NOW MOVE JUST THE RELEVANT FILES
(GP PARAMS AND FIXED PDBS)
TO THEIR OWN DIRECTORY


cd genpot
mkdir clean
mv *clean.pdb clean
mv *.params clean
cd clean

arg list too long

import os
l=[i for i in os.listdir() if i[-9:]=='clean.pdb']
print(len(l))
#45673
c=0
for i in l:
    os.system('mv '+i+' clean/'+i)
    c+=1
    print(str(c))
l=[i for i in os.listdir() if i[-7:]=='.params']
c=0
for i in l:
    os.system('mv '+i+' clean/'+i)
    c+=1
    print(str(c))

echo "
import os
l=[i for i in os.listdir() if i[-7:]=='.params']
c=0
for i in l:
    os.system('mv '+i+' clean/'+i)
    c+=1
    print(str(c))
">moveparams.py

echo "
#!/bin/bash
time python3 moveparams.py
">moveparams.sh

qsub -cwd -l mem_free=4G moveparams.sh

/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean
'''


import os
#have to fix the params files to make sure no redundant bonds defined :(
prms=[i for i in os.listdir() if i[-6:]=='params']
c=0
for prm in prms:
    f=open(prm,'r')
    lines=[line for line in f.readlines()]
    f.close()
    batoms=[]
    for ind,line in enumerate(lines):
        if line[:9]=='BOND_TYPE':
            a1=line[10:13].strip(' ')
            a2=line[15:18].strip(' ')
            batoms.append(((a1,a2),ind))
    toremove=[]
    for ix,e1 in enumerate(batoms):
        for e2 in batoms[ix+1:]:
            if e1!=e2:
                if e1[0]==e2[0]:
                    toremove.append(e2[1])
    toremove=sorted(toremove,reverse=True)
    for idx in toremove:
        lines.pop(idx)
    of=open(prm,'w')
    for line in lines:
        of.write(line)
    of.close()
    c+=1
    print(str(c))





#resfile generation
#################################################################################
#################################################################################
#################################################################################
#################################################################################
import os
#################################################################################
allparams=[i for i in os.listdir() if i[-6:]=='params']
prm1=allparams[0]
sfname='srmr.sh'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
45673
'''
count=1
idx=[j for j in range(0,len(matches),10)]  #how many matches per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    of.write('source ~/anaconda3/etc/profile.d/conda.sh')
    of.write('\n')
    of.write('conda activate pyr37')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        cmdd='time python ~/BSFF/tools/match_resfile.py '+matchpath+' '+prm1+' 4.0 genpot'
        of.write(cmdd+'\n')
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    print(count)
    count+=1
os.makedirs('cluster_out',exist_ok=True)
#####
#####
######
import os
sfname='srmr.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+j
    os.system(cmd)
'''
script is identifying upwards of 20 residues as designable even though
using 4 angstrom cutoff, strikes me as very odd and too many res,
maybe i will try design without resfile and
just using enzdes packer task? maybe both methods?
maybe its happening because in the match theres a lotta clashes and a kinda
cluttered binding site?


import os
l=[i for i in os.listdir() if i[-7:]=='resfile']
for i in l:
    os.system('rm '+i)
'''







# # changing cst files to match the pdb of matches names
##and move it to directory with pdbs,params,resfiles
import os
################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
cstdir='/wynton/home/kortemme/cgalvin/esl/1np3hb/1np'
################################################
for i in matches:
    ################################################
    cstid='hybrid_'+i.split('_')[-4]+'.cst'####################################################
    ################################################
    ogcst=open(os.path.join(cstdir,cstid),'r')
    ################################################
    lines=[line for line in ogcst.readlines()]
    ogcst.close()
    newcstname=i.split('.')[0]+'.cst'
    of=open(newcstname,'w')
    for line in lines:
        of.write(line)
    of.close()

'''
import os
l=[i for i in os.listdir() if i[-3:]=='cst']
for i in l:
    os.system('rm '+i)

'''
























#
#
#
#
#
#
#
#
# # FASTDESIGN ON MATCHES
# #1 params for all
# # resfile and cst must have same name as pdb
# import os
# ########
# allparams=[i for i in os.listdir() if i[-6:]=='params']
# prm1=allparams[0]
# sfname='fd1.sh'
# ndesignseachmatch='50'
# outputdirectory='fd'
# scorefile_name='dog_fd.json'
# #########
# matches=[i for i in os.listdir() if i[-3:]=='pdb']
#
# os.makedirs('cluster_output',exist_ok=True)
# os.makedirs(outputdirectory,exist_ok=True)
# c=1
# for xx in range(1):     #this way i can spread across more nodes if neccessary
#     output_prefix='fd'+str(c)
#     sf=open('_'+str(c)+sfname,'w')
#     sf.write('#!/bin/bash')
#     sf.write('\n')
#     sf.write('tasks=(0\n')
#     for match in matches[:-1]:
#         sf.write('       '+str(match.strip('.pdb'))+'\n')
#     sf.write('       '+matches[-1].strip('.pdb')+')')
#     sf.write('\n')
#     cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
#          '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
#          '-parser:protocol','~/BSFF/tools/FD_jump1_genpot.xml',
#          '-parser:view','-run:preserve_header',
#          '-in:file::extra_res_fa',prm1,
#          '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
#          '-ex1','-ex2','-extrachi_cutoff','0',
#          '-out:nstruct',ndesignseachmatch,
#          # '-score:set_weights','hbond_bb_sc','2.0',
#          # '-score:set_weights','hbond_sc','2.0',
#          '-score::weights beta_genpot',
#          '-corrections::gen_potential',
#          '-load_PDB_components','False',
#          '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
#          '-out:path:all',outputdirectory,
#          '-out:prefix',output_prefix,
#          '-scorefile_format','json',
#          '-out:file:scorefile',scorefile_name]
#     sf.write((' ').join(cmd))
#     sf.write('\nqstat -j "$JOB_ID"')
#     sf.close()
#     c+=1
# '''
# print(len(matches))
# 3657
#
#
# ~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.pdb -parser:protocol ~/BSFF/tools/FD_jump1_genpot.xml -parser:view -run:preserve_header -in:file::extra_res_fa UM_10_N107S110T17Y58_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_cst_4_1_0001_ligH_bcc.params -resfile UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.resfile -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.cst -out:path:all fd -out:prefix fd1 -scorefile_format json -out:file:scorefile dog_fd.json
# '''
# import os
# sfname='fd1.sh'
# jobss=[i for i in os.listdir() if sfname in i]
# for j in jobss:
#     cmd='qsub -cwd -t 1-3657 -l mem_free=8G -o cluster_out -e cluster_out '+j
#     os.system(cmd)
#
# '''
# these jobs maxvmem is actually under 1g, i could def safely specify only 4g
# for them and be good
# LOOKS LIKE THEY TAKE ABOUT 30-60 HOURS TO COMPLETE, THATS A PRETTY LONG TIME,
# IDK MAYBE I SHOULD CHILL ON THE 50 DES EACH FOR NEXT ROUND
#
# beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
#
#
#
# I THINK A REASONABLE NEXT STEP IS TO FILTER BASED ON INTERFACEBUNS/LIGHPHOBESASA
# AND THEN MAKE OG PARAMS AND DO DDG/BOLTZ ANALYSIS ON THOSE FILTERED DESIGNS ONLY
#
# ADDITIONALLY WANNA SET UP
# MPNN DESIGN SPACE
# 3BOP
# SPECIAL ROTAMERS (study steroid binders in pdb)
#
# beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
#       user/group     ||           size          ||    chunk files
#      name     |  id  ||    used    |    hard    ||  used   |  hard
# --------------|------||------------|------------||---------|---------
#        cgalvin| 61046||  543.72 GiB| 1000.00 GiB||  2188288|unlimited
# '''
#
#
#
# '''
# FASTDESIGN with 3bop weight 5, intup weight 2
# FASTDESIGN with 3bop weight 5, intup weight 2
# FASTDESIGN with 3bop weight 5, intup weight 2
# '''
# # FASTDESIGN with 3bop weight 5, intup weight 2
# import os
# ########
# allparams=[i for i in os.listdir() if i[-6:]=='params']
# prm1=allparams[0]
# sfname='fd_3bop.sh'
# ndesignseachmatch='25'
# outputdirectory='fd_3bop'
# scorefile_name='dog_fd_3bop.json'
# #########
# matches=[i for i in os.listdir() if i[-3:]=='pdb']
#
# os.makedirs('cluster_output',exist_ok=True)
# os.makedirs(outputdirectory,exist_ok=True)
# c=1
# for xx in range(2):     #this way i can spread across more nodes if neccessary
#     output_prefix='fd'+str(c)
#     sf=open('_'+str(c)+sfname,'w')
#     sf.write('#!/bin/bash')
#     sf.write('\n')
#     sf.write('tasks=(0\n')
#     for match in matches[:-1]:
#         sf.write('       '+str(match.strip('.pdb'))+'\n')
#     sf.write('       '+matches[-1].strip('.pdb')+')')
#     sf.write('\n')
#     cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
#          '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
#          '-parser:protocol','~/BSFF/tools/FD_jump1_genpot_3bop.xml',
#          '-parser:view','-run:preserve_header',
#          '-in:file::extra_res_fa',prm1,
#          '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
#          '-ex1','-ex2','-extrachi_cutoff','0',
#          '-out:nstruct',ndesignseachmatch,
#          # '-score:set_weights','hbond_bb_sc','2.0',
#          # '-score:set_weights','hbond_sc','2.0',
#          '-score::weights beta_genpot',
#          '-corrections::gen_potential',
#          '-load_PDB_components','False',
#          '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
#          '-out:path:all',outputdirectory,
#          '-out:prefix',output_prefix,
#          '-scorefile_format','json',
#          '-out:file:scorefile',scorefile_name]
#     sf.write((' ').join(cmd))
#     sf.write('\nqstat -j "$JOB_ID"')
#     sf.close()
#     c+=1
# '''
# print(len(matches))
# 3657
#
# ~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.pdb -parser:protocol ~/BSFF/tools/FD_jump1_genpot_3bop.xml -parser:view -run:preserve_header -in:file::extra_res_fa UM_10_N107S110T17Y58_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_cst_4_1_0001_ligH_bcc.params -resfile UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.resfile -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.cst -out:path:all fd_3bop -out:prefix fd2 -scorefile_format json -out:file:scorefile dog_fd_3bop.json
#
# '''
# import os
# sfname='fd_3bop.sh'
# jobss=[i for i in os.listdir() if sfname in i]
# for j in jobss:
#     cmd='qsub -cwd -t 1-3657 -l mem_free=4G -l h_rt=25:00:00 -o cluster_out -e cluster_out '+j
#     os.system(cmd)
#
#
# '''
# analysis of fd designs
#     setup ddg/boltz calc
#     submit filtered to cf
#     setup mpnn assisted stabilization
# bsff on steroid carbon skeleton + special rotamers setup
#
# '''
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# coupledmoves
# setup so that 1 params for all
# resfile and cst must have same name as pdb
########
import os
inputs_dir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
# ndesignseachmatch='50'
ndesignseachmatch='50'
sfname='cm.sh'
outputdirectory='low_designs'
scorefile_name='cmdesigns.json'
#
matches1=[i for i in os.listdir('/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res/has_aro_hbond/') if i[-3:]=='pdb']
matches=[os.path.join(inputs_dir,i) for i in matches1]
print(len(matches))



#########
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(1):
    output_prefix='cm'+str(c)
    sf=open('_'+str(c)+sfname,'w')
    sf.write('#!/bin/bash')
    sf.write('\n')
    sf.write('tasks=(0\n')
    for match in matches[:-1]:
        sf.write('       '+str(match.strip('.pdb'))+'\n')
    sf.write('       '+matches[-1].strip('.pdb')+')')
    sf.write('\n')
    cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
         '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
         '-parser:protocol','~/BSFF/tools/cm_genpot_jump1.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-mute protocols.backrub.BackrubMover',
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         # '-score::weights beta_genpot',
         '-corrections::gen_potential',
         # '-score::weights','ligand',
         # '-restore_pre_talaris_2013_behavior',
         # '-score:set_weights','hbond_bb_sc','2.0',
         # '-score:set_weights','hbond_sc','2.0',
         '-score:set_weights','hbnet','0.5',
         # '-score:set_weights','voids_penalty','0.5',
         # '-score:set_weights','buried_unsatisfied_penalty','0.5',
         # '-score:set_weights','approximate_buried_unsat_penalty','5.0',
         '-load_PDB_components','False',
         '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name,
         '-coupled_moves::mc_kt 2.4',
         '-coupled_moves::boltzmann_kt 2.4', '-coupled_moves::ntrials 1000',
         '-coupled_moves::initial_repack false', '-coupled_moves::ligand_mode true',
         # '-coupled_moves::ligand_weight 2',
         '-coupled_moves::fix_backbone false',
         '-coupled_moves::bias_sampling true', '-coupled_moves::bump_check true',
         '-coupled_moves::exclude_nonclashing_positions true']
         # '-coupled_moves::backbone_mover kic',
         # '-coupled_moves::kic_perturber walking']
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1
'''


~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.pdb -parser:protocol ~/BSFF/tools/cm_genpot_jump1.xml -parser:view -run:preserve_header -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -resfile /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.resfile -mute protocols.backrub.BackrubMover -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 50 -corrections::gen_potential -score:set_weights hbnet 0.5 -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.cst -out:path:all low_designs -out:prefix cm1 -scorefile_format json -out:file:scorefile cmdesigns.json -coupled_moves::mc_kt 2.4 -coupled_moves::boltzmann_kt 2.4 -coupled_moves::ntrials 1000 -coupled_moves::initial_repack false -coupled_moves::ligand_mode true -coupled_moves::fix_backbone false -coupled_moves::bias_sampling true -coupled_moves::bump_check true -coupled_moves::exclude_nonclashing_positions true


'''
import os
sfname='cm.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4015 -l mem_free=2G -o cluster_output -e cluster_output '+j
    os.system(cmd)




'''
rm *last.pdb
rm *low.pdb
rm *.fasta
rm *.stats
-bash: /usr/bin/rm: Argument list too long



import os
l=[i for i in os.listdir() if i[-8:]=='last.pdb']
l2=[i for i in os.listdir() if i[-7:]=='low.pdb']
l3=[i for i in os.listdir() if i[-6:]=='.fasta']
c=1
print(len(l))
for i in l:
    os.system('rm '+i)
    c+=1
    print(str(c)+'/'+str(len(l)))
c=1
for i in l2:
    os.system('rm '+i)
    c+=1
    print(str(c)+'/'+str(len(l2)))
c=1
for i in l3:
    os.system('rm '+i)
    c+=1
    print(str(c)+'/'+str(len(l3)))


echo "
import os
l=[i for i in os.listdir() if i[-8:]=='last.pdb']
l2=[i for i in os.listdir() if i[-7:]=='low.pdb']
l3=[i for i in os.listdir() if i[-6:]=='.fasta']
l4=[i for i in os.listdir() if i[-6:]=='.stats']

c=1
for i in l:
    os.system('rm '+i)
c=1
for i in l2:
    os.system('rm '+i)
c=1
for i in l3:
    os.system('rm '+i)
c=1
for i in l4:
    os.system('rm '+i)
">cleancm.py

echo '
#!/bin/bash
time python3 cleancm.py
qstat -j "$JOB_ID"
'>cleancm.sh
qsub -cwd -l mem_free=4G cleancm.sh

'''







#analysis of scores from json file
sfname='cmdesigns.json'

#
# import json
# bad=[]
# scores=[]
# for line in open(sfname,'r'):
#     try:
#         scores.append(json.loads(line))
#     except:
#         bad.append(line)
# for line in bad:
#     ids=[]
#     for ii,i in enumerate(line):
#         if i=='}':
#             ids.append(ii)
#     if len(ids)==2:
#         l1=line[:ids[0]+1]
#         l2=line[ids[0]+1:]
#         scores.append(l1)
#         scores.append(l2)
# terms=list(scores[0].keys())
# print(len(bad))
'''
25082
'''


import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
print(len(bad))
print(len(scores))
ic=0
for d in scores:
    # print(ic)
    # ic+=1
    try:
        if len(d.keys())<5:
            scores.remove(d)
    except:
        scores.remove(d)
print(len(scores))
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==0:
        if len(line)>100:
            if line[-1]=='}':
                if line[0]!='{':
                    ll='{'+line
                    line=ll
                if line[0]=='{' and line[1]=='{':
                    ll=line[1:]
                    line=ll
                try:
                    scores.append(json.loads(line))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==1:
        l1=line[:ids[0]+1].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==2:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==3:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==4:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==5:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==6:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:ids[4]+1].strip('\n')
        l6=line[ids[4]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        if len(l6)>100:
            if l6[-1]=='}':
                if l6[0]!='{':
                    ll='{'+l6
                    l6=ll
                if l6[0]=='{' and l6[1]=='{':
                    ll=l6[1:]
                    l6=ll
                try:
                    scores.append(json.loads(l6))
                except:
                    pass
        bad.remove(line)
terms=list(scores[0].keys())
print(len(bad))
print(len(scores))
'''
bad
scores
scores
bad
scores

37079
152599
152599
18539
169897

'''





#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                try:
                    allscores.append(float(d[term]))
                except:
                    print(d)
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

# plot_dists(terms,scores,'esl_4rm_batch4.pdf')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

f1=return_filtered(scores,'ligbuns','<',0.0)
# f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'lighphobesasa','<',30.0)
f3=return_filtered(f2,'hbtolig','>',3.0)
f4=return_filtered(f3,'liginterface','<',-25.)
# plot_dists(terms,f3,'esl4rmfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))

'''
ligbuns
169897
117273
55444
3552
3551

buns2
169897
28337
10352
1188
1188
'''
import os
filtered_strc=[]
for d in f4:
    nn=d['decoy'][:-5]
    filtered_strc.append(nn)
os.makedirs('filtered',exist_ok=True)
print(len(filtered_strc))
c=0
for i in filtered_strc:
    os.system('cp '+i+'_0001.pdb filtered/'+i+'_0001.pdb')
    c+=1
    print(str(c))
os.system('mv *.pdf filtered/')

























'''
cd enzdes

cd ed2

mkdir ed3
cd ed3

mkdir ed6
cd ed6

i guess i changed batch dict to be 5k each instead of 15? so
ill start with batch 4 here

fuck i need to change the path in the batch dict from esl to esl1

ADDED A FEW NEW OPTIONS HERE TO TEST OUT
    but it looks like if i wanna disable pro/met i need to do so via resfile
    annoying but no way to do through commandline as far as i can tell
'''
#########
########
########
#ENZYME DESIGN
# import os

# #
# # import json
# # f=open('/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/design_batches.json','r')
# # batch_dict=json.load(f)
# # f.close()
# # ########
# # ########
# # ########
# # # matches1=batch_dict['batch_4']
# # matches1=batch_dict['batch_6']
# #
# # #
# # matches=[]
# # for i in matches1:
# #     strc=i.split('/')[-1]
# #     matches.append(os.path.join('/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/',strc))
# # #
# # print(len(matches))

########
########
########
#
########
import os
inputs_dir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
# ndesignseachmatch='50'
ndesignseachmatch='10'
sfname='enzdes2.sh'

#
matches1=[i for i in os.listdir('/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res/has_aro_hbond/') if i[-3:]=='pdb']
matches=[os.path.join(inputs_dir,i) for i in matches1]
print(len(matches))
#
#########
###########
try:
    os.system('rm -r cluster_output')
except:
    pass
os.makedirs('cluster_output',exist_ok=True)
c=1
output_prefix='ed2'+str(c)
sf=open(sfname,'w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('tasks=(0\n')
for match in matches[:-1]:
    sf.write('       '+str(match.strip('.pdb'))+'\n')
sf.write('       '+matches[-1].strip('.pdb')+')')
sf.write('\n')
cmd=['~/main/source/bin/enzyme_design.linuxgccrelease',
     '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
     '-in:file::extra_res_fa',prm1,
     '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
     '-ex1','-ex2','-extrachi_cutoff','0',
     '-enzdes:cst_opt True','-enzdes:bb_min True',
     '-enzdes:chi_min True',
     '-enzdes:cst_design True','-enzdes:design_min_cycles 4', '-enzdes:lig_packer_weight 2.5',
     '-enzdes:cst_min True',
     # '-enzdes:detect_design_interface',
     # '-enzdes:cut1','6.0',
     # '-enzdes:cut2','8.0',
     # '-enzdes:cut3','10.0',
     # '-enzdes:cut4','12.0',
     # '-score:set_weights','hbond_bb_sc','2.0',
     # '-score:set_weights','hbond_sc','2.0',
     '-score:set_weights','hbnet','1.0',
     # '-score:set_weights','voids_penalty','0.5',
     # '-score:set_weights','buried_unsatisfied_penalty','0.5',
     '-score:set_weights','approximate_buried_unsat_penalty','5.0',
     # '-packing:soft_rep_design True',
     '-relax:script InterfaceDesign2019',
     '-no_packstat_calculation True',
     '-out:nstruct',ndesignseachmatch,
     '-packing:linmem_ig 10',
     '-corrections::gen_potential',
     '-load_PDB_components','False',
     '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
     '-out:prefix',output_prefix]
sf.write((' ').join(cmd))
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
c+=1
print(len(matches))
'''

-script <File>
Relax script file
Default: ""


~/main/source/bin/enzyme_design.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_3_L58W65S99W55_1_relaxed_relaxed_5tpj_224592_design_5_unrelaxed_model_1_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_hybrid_42_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt -enzdes:bb_min -enzdes:chi_min -enzdes:cst_design -enzdes:design_min_cycles 4 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min -enzdes:detect_design_interface -enzdes:cut1 6.0 -enzdes:cut2 8.0 -enzdes:cut3 10.0 -enzdes:cut4 12.0 -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_3_L58W65S99W55_1_relaxed_relaxed_5tpj_224592_design_5_unrelaxed_model_1_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_hybrid_42_1_0001_clean.cst -out:prefix ed1

approximate_buried_unsat_penalty

~/main/source/bin/enzyme_design.linuxgccrelease -s  /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -resfile  /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.resfile -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt True -enzdes:bb_min True -enzdes:chi_min True -enzdes:cst_design True -enzdes:design_min_cycles 4 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min True -score:set_weights hbnet 1.0 -score:set_weights approximate_buried_unsat_penalty 5.0 -packing:soft_rep_design True -no_packstat_calculation True -out:nstruct 5 -packing:linmem_ig 10 -corrections::gen_potential -load_PDB_components False -enzdes:cstfile  /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_15_F106W56T13S34_1_relaxed_relaxed_5tpj_169560_design_2_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean.cst -out:prefix ed21


/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/

~/main/source/bin/enzyme_design.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_29_M54Y69S110N38_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_1_rank_1_0001_hybrid_5_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt -enzdes:bb_min -enzdes:chi_min -enzdes:cst_design -enzdes:design_min_cycles 4 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min -enzdes:detect_design_interface -enzdes:cut1 6.0 -enzdes:cut2 8.0 -enzdes:cut3 10.0 -enzdes:cut4 12.0 -out:nstruct 50 -corrections::gen_potential -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_29_M54Y69S110N38_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_1_rank_1_0001_hybrid_5_1_0001_clean.cst -out:prefix ed21
'''
import os
sfname='enzdes2.sh'
cmd='qsub -cwd -t 1-4015 -l mem_free=2G -o cluster_output -e cluster_output '+sfname
os.system(cmd)

'''
beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed2

rm -r ...ed1 and esl

REMEMBER TO EMPLOY STRICTER HB FILTER!
ALSO WANT TO CHANGE RESFILE GENERATION SCRIPT TO JUST EXCLUDE C/P/(M?)
    maybe next batch i will use this instead of automatic int detection?


batch 5 no cpm
    SEEMS LIKE IN GENERAL I GET A FEW CORE DUMPS UP FRONT BUT THEN ITS ALL GOOD
        just gotta be attentive for a minute to delete them all

/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res/has_aro_hbond

'''




#
#
#
#
#
#
# import os
# import json
# outputdir='enzdes'
# if os.path.exists(outputdir):
#     os.chdir(outputdir)
# else:
#     os.makedirs(outputdir)
#     os.chdir(outputdir)
# ########
# inputs_dir='/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean'
# allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
# prm1=allparams[0]
# ndesignseachmatch='10'
# sfname='enzdes.sh'
# ##########
# matches=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-3:]=='pdb']
# print(len(matches))
# count=1
# idx=[j for j in range(0,len(matches),20)]  #how many matches per job
# if len(matches) not in idx:
#     idx.append(len(matches))
# for ei, ix in enumerate(idx[:-1]):
#     sfn='_'+str(count)+sfname
#     of=open(sfn,'w')
#     of.write('#!/bin/bash\n')
#     of.write('\n')
#     for matchpath in matches[ix:idx[ei+1]]:
#         cmd=['~/main/source/bin/enzyme_design.linuxgccrelease',
#              '-s',matchpath,         ##
#              '-in:file::extra_res_fa',prm1,
#              # '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
#              '-ex1','-ex2','-extrachi_cutoff','0',
#              '-enzdes:cst_opt','-enzdes:bb_min',
#              '-enzdes:chi_min',
#              '-enzdes:cst_design','-enzdes:design_min_cycles 4', '-enzdes:lig_packer_weight 2.5',
#              '-enzdes:cst_min',
#              '-enzdes:detect_design_interface',
#              '-enzdes:cut1','4.0',
#              '-enzdes:cut2','6.0',
#              '-enzdes:cut3','8.0',
#              '-enzdes:cut4','10.0',
#              '-out:nstruct',ndesignseachmatch,
#              '-score::weights beta_genpot',
#              '-corrections::gen_potential',
#              '-load_PDB_components','False',
#              '-enzdes:cstfile',matchpath.strip('.pdb')+'.cst']
#         of.write(' '.join(cmd)+'\n')
#     of.write('\nqstat -j "$JOB_ID"')
#     of.close()
#     print(count)
#     count+=1
# try:
#     os.system('rm -r cluster_output')
# except:
#     pass
# os.makedirs('cluster_output',exist_ok=True)
# '''
# ~/main/source/bin/enzyme_design.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_3_M13T82S107S53_1_relaxed_relaxed_5tpj_18512_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_design_7_unrelaxed_model_5_rank_1_0001_hybrid_39_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt -enzdes:bb_min -enzdes:chi_min -enzdes:cst_design -enzdes:design_min_cycles 4 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min -enzdes:detect_design_interface -enzdes:cut1 6.0 -enzdes:cut2 8.0 -enzdes:cut3 10.0 -enzdes:cut4 12.0 -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/UM_3_M13T82S107S53_1_relaxed_relaxed_5tpj_18512_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_design_7_unrelaxed_model_5_rank_1_0001_hybrid_39_1_0001_clean.cst
# '''
# import os
# sfname='enzdes.sh'
# jobs=[i for i in os.listdir() if sfname in i]
# for j in jobs:
#     cmd='qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output '+j
#     os.system(cmd)
# '''
# getting tons of core dumps, is it because of the automatic interface detection or something?
# wtf is going on here?
#     TERMINATING UNTIL I CAN FIGURE OUT WHAT THE FUCKING PROBLEMN IS
#     i will try to make the resfiles old fashioned way and submit with jobs
#         EH BUT MAYBE I SHOULD WAIT FOR ANALYSIS RESULTS FIRST FOR THE 250K DESIGNS
#         I ALREADY MADE, SEE IF THEY LOOK AT ALL WORTHWHILE
# cus the other thing is here that soon i will have expanded scaffold library and then
# maybe i can run back 5 res motif matching with this and get expanded set of designs
# that way
#     or else maybe i can figure out this generalized pi stacking cst idea
#
# '''








'''
LOOKING AT DESIGNS WHERE I DID ARRAY JOB THAT I STOPPED EARLY

/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1

/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed2

'''


import os
#################################################################################
inputs_dir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
sfname='qa.sh'
scorefile_name='enzdesanalysis.json'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
243698
249094
19914

'''
count=1
idx=[j for j in range(0,len(matches),10)]  #how many matches per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
             '-s',matchpath,
             '-parser:protocol','~/BSFF/tools/edanalysis.xml',
             '-parser:view','-run:preserve_header',
             '-in:file::extra_res_fa',prm1,
             '-score::weights beta_genpot',
             '-corrections::gen_potential',
             '-load_PDB_components','False',
             '-scorefile_format','json',
             '-out:file:score_only',scorefile_name]
        of.write((' ').join(cmd))
        of.write('\n')
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    print(count)
    count+=1
os.makedirs('cluster_out',exist_ok=True)
#####
#####
######
import os
sfname='qa.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+j
    os.system(cmd)
'''
/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3


'''





'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE

'''




#analysis of scores from json file
sfname='enzdesanalysis.json'

#
# import json
# bad=[]
# scores=[]
# for line in open(sfname,'r'):
#     try:
#         scores.append(json.loads(line))
#     except:
#         bad.append(line)
# for line in bad:
#     ids=[]
#     for ii,i in enumerate(line):
#         if i=='}':
#             ids.append(ii)
#     if len(ids)==2:
#         l1=line[:ids[0]+1]
#         l2=line[ids[0]+1:]
#         scores.append(l1)
#         scores.append(l2)
# terms=list(scores[0].keys())
# print(len(bad))
'''
25082
'''


import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
print(len(bad))
print(len(scores))
ic=0
for d in scores:
    # print(ic)
    # ic+=1
    try:
        if len(d.keys())<5:
            scores.remove(d)
    except:
        scores.remove(d)
print(len(scores))
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==0:
        if len(line)>100:
            if line[-1]=='}':
                if line[0]!='{':
                    ll='{'+line
                    line=ll
                if line[0]=='{' and line[1]=='{':
                    ll=line[1:]
                    line=ll
                try:
                    scores.append(json.loads(line))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==1:
        l1=line[:ids[0]+1].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==2:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==3:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==4:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==5:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==6:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:ids[4]+1].strip('\n')
        l6=line[ids[4]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        if len(l6)>100:
            if l6[-1]=='}':
                if l6[0]!='{':
                    ll='{'+l6
                    l6=ll
                if l6[0]=='{' and l6[1]=='{':
                    ll=l6[1:]
                    l6=ll
                try:
                    scores.append(json.loads(l6))
                except:
                    pass
        bad.remove(line)
terms=list(scores[0].keys())
print(len(bad))
print(len(scores))
'''
bad
scores
scores
bad
scores
64804
144502
144502
32409
182563

b5
59186
141101
141092
29631
175261

b6
2097
46212
46212
1048
47387
'''





#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                try:
                    allscores.append(float(d[term]))
                except:
                    print(d)
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

# plot_dists(terms,scores,'esl_4rm_batch4.pdf')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'lighphobesasa','<',40.0)
f3=return_filtered(f2,'hbtolig','>',3.0)
f4=return_filtered(f3,'ligoversat','<',0.0)
f5=return_filtered(f4,'bsE_per_res','<',-2.0)
# plot_dists(terms,f3,'esl4rmfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
'''
235413
36795
24301
19609

batch4
182563
30426
20659
16699
16683
12027

b5
175261
24165
12829
7734
7723
5490

b6
47387
3883
1427
533
530
97

des1
19329
1970
527
226
225
51

des2
17044
1843
486
274
271
53
'''
import os
filtered_strc=[]
for d in f4:
    nn=d['decoy'][:-5]
    filtered_strc.append(nn)
os.makedirs('filtered',exist_ok=True)
print(len(filtered_strc))
c=0
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
    c+=1
    print(str(c))
os.system('mv *.pdf filtered/')

# l=[i for i in os.listdir() if i[-3:]=='pdb']
# l[:100]
# filtered_strc[:100]
'''
cd filtered

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/esl_4rm.pdf ~/desktop/esl_4rm.pdf
'''


'''
ADD CORRECT DIRECTORY TO BOTTOM OF CODE HERE (EXCEPT LINE)
cd filtered
pwd
/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3/filtered

echo "
import os
import sys
ligname='esl'
cleandirname='analysis'
#####################
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(pdbs))
os.makedirs(cleandirname,exist_ok=True)
for initial_match in pdbs:
    try:
        fname=initial_match.split('.')[0]
        fname2=os.path.join(cleandirname,initial_match.split('.')[0]+'_lig.pdb')
        f=open(initial_match,'r')
        lines=[]
        alines=[]
        for line in f.readlines():
            if line[0:6]=='HETATM':
                lines.append(line)
            elif line[0:4]=='ATOM':
                alines.append(line)
            elif line[:6]=='REMARK':
                alines.append(line)
        f.close()
        newligpdb=open(fname2,'w')
        for line in lines:
            newligpdb.write(line)
        newligpdb.close()
        os.chdir(cleandirname)
        os.system('obabel -i pdb '+fname+'_lig.pdb'+' -o mol -O '+fname+'_lig.mol')
        os.system('obabel -i mol '+fname+'_lig.mol'+' -o mol -O '+fname+'_ligH.mol -p 7.4')
        os.system('~/main/source/scripts/python/public/molfile_to_params.py -n '+ligname+' -p '+fname+' '+fname+'_ligH.mol')
        newlig=open(fname+'_0001.pdb','r')
        newliglines=[line for line in newlig.readlines() if line[:6]=='HETATM']
        newpdb=open(os.path.join(fname+'_oglig.pdb'),'w')
        for line in alines:
            newpdb.write(line)
        newpdb.write('TER\n')
        for line in newliglines:
            newpdb.write(line)
        newpdb.write('TER\n')
        newpdb.close()
        os.chdir('..')
    except:
        os.chdir('/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res/has_aro_hbond/cm1/low_designs/filtered')
">makeparams.py

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 makeparams.py
qstat -j "$JOB_ID"
'>makeparams.sh
qsub -cwd -l mem_free=4G makeparams.sh







'''

'''
cd analysis
mkdir run
mv *oglig.pdb run
mv *.params run
cd run


import os
l=[i for i in os.listdir() if i[-9:]=='oglig.pdb']
print(len(l))
c=0
for i in l:
    os.system('mv '+i+' run/'+i)
    c+=1
    print(c)
l=[i for i in os.listdir() if i[-7:]=='.params']
print(len(l))
c=0
for i in l:
    os.system('mv '+i+' run/'+i)
    c+=1
    print(c)


cd run
'''
# #run the extra scoring
# import os
# params=[i for i in os.listdir() if i[-6:]=='params']
# paramspath=params[0]
#
# pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
# #
# sf=open('extrascore.sh','w')
# sf.write('#!/bin/bash')
# sf.write('\n')
# sf.write('tasks=(0\n')
# for match in pdbs[:-1]:
#     sf.write('       '+match+'\n')
# sf.write('       '+pdbs[-1]+')')
# sf.write('\n')
# sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]} -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+paramspath+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
# sf.write('\nqstat -j "$JOB_ID"')
# sf.close()
# #
# print(len(pdbs))
'''
qsub -cwd -t 1-1805 -l mem_free=2G extrascore.sh
        took about an hour per job


'''









import os
#################################################################################
params=[i for i in os.listdir() if i[-6:]=='params']
prm1=params[0]
sfname='extrascore.sh'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
batch 1-3
243698

batch 4
11990

5
5479
'''
count=1
idx=[j for j in range(0,len(matches),1)]  #how many per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        of.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s='+matchpath+' -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+prm1+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10  -beta_nov16')
        of.write('\n')
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    print(count)
    count+=1
os.makedirs('cluster_out',exist_ok=True)
#####
#####
'''
/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=ed21UM_5_F51T43S69Y67_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_84_1_0001_clean__DE_13_oglig.pdb -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa ed21UM_10_L57Y30Q43T98_1_relaxed_relaxed_5tpj_384973_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_1_1_0001_clean__DE_9.params -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10  -beta_nov16
'''
######
import os
sfname='extrascore.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+j
    os.system(cmd)
























'''
/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed2/filtered/analysis/run
'''

#analysis of scores from json file
sfname='scores.json'
import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
print(len(bad))
print(len(scores))
ic=0
for d in scores:
    # print(ic)
    # ic+=1
    try:
        if len(d.keys())<5:
            scores.remove(d)
    except:
        scores.remove(d)
print(len(scores))
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==0:
        if len(line)>100:
            if line[-1]=='}':
                if line[0]!='{':
                    ll='{'+line
                    line=ll
                if line[0]=='{' and line[1]=='{':
                    ll=line[1:]
                    line=ll
                try:
                    scores.append(json.loads(line))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==1:
        l1=line[:ids[0]+1].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==2:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==3:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==4:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==5:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        bad.remove(line)
    elif len(ids)==6:
        l1=line[:ids[0]+1].strip('\n')
        l2=line[ids[0]+1:ids[1]+1].strip('\n')
        l3=line[ids[1]+1:ids[2]+1].strip('\n')
        l4=line[ids[2]+1:ids[3]+1].strip('\n')
        l5=line[ids[3]+1:ids[4]+1].strip('\n')
        l6=line[ids[4]+1:].strip('\n')
        if len(l1)>100:
            if l1[-1]=='}':
                if l1[0]!='{':
                    ll='{'+l1
                    l1=ll
                if l1[0]=='{' and l1[1]=='{':
                    ll=l1[1:]
                    l1=ll
                try:
                    scores.append(json.loads(l1))
                except:
                    pass
        if len(l2)>100:
            if l2[-1]=='}':
                if l2[0]!='{':
                    ll='{'+l2
                    l2=ll
                if l2[0]=='{' and l2[1]=='{':
                    ll=l2[1:]
                    l2=ll
                try:
                    scores.append(json.loads(l2))
                except:
                    pass
        if len(l3)>100:
            if l3[-1]=='}':
                if l3[0]!='{':
                    ll='{'+l3
                    l3=ll
                if l3[0]=='{' and l3[1]=='{':
                    ll=l3[1:]
                    l3=ll
                try:
                    scores.append(json.loads(l3))
                except:
                    pass
        if len(l4)>100:
            if l4[-1]=='}':
                if l4[0]!='{':
                    ll='{'+l4
                    l4=ll
                if l4[0]=='{' and l4[1]=='{':
                    ll=l4[1:]
                    l4=ll
                try:
                    scores.append(json.loads(l4))
                except:
                    pass
        if len(l5)>100:
            if l5[-1]=='}':
                if l5[0]!='{':
                    ll='{'+l5
                    l5=ll
                if l5[0]=='{' and l5[1]=='{':
                    ll=l5[1:]
                    l5=ll
                try:
                    scores.append(json.loads(l5))
                except:
                    pass
        if len(l6)>100:
            if l6[-1]=='}':
                if l6[0]!='{':
                    ll='{'+l6
                    l6=ll
                if l6[0]=='{' and l6[1]=='{':
                    ll=l6[1:]
                    l6=ll
                try:
                    scores.append(json.loads(l6))
                except:
                    pass
        bad.remove(line)
terms=list(scores[0].keys())
print(len(bad))
print(len(scores))
'''
bad
scores
scores
bad
scores
13
11699
11699
6
11708

batch5
1
5421
5421
0
5422
'''
#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                try:
                    allscores.append(float(d[term]))
                except:
                    print(d)
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

# plot_dists(terms,scores,'4rm_testrun.pdf')
'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/dog_idel_fd_filt_extrascores.pdf ~/desktop/dog_idel_fd_filt_extrascores.pdf
'''
def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

#
f1=return_filtered(scores,'buns2interface','<',10.0)
f2=return_filtered(f1,'boltz','<',-0.3)
f3=return_filtered(f2,'ddg','<',-24.0)
f4=return_filtered(f3,'contact_molsurf','>',160)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.65)
f7=return_filtered(f6,'lighphobesasa','<',40.0)
f8=return_filtered(f7,'packstat','>',0.6)
f9=return_filtered(f8,'buns_bb_heavy','<',5)
f10=return_filtered(f9,'buns_sc_heavy','<',1)
f11=return_filtered(f10,'ligoversat','<',0)
f12=return_filtered(f11,'oversat','<',0)
f13=return_filtered(f12,'exphyd','<',1200)
f14=return_filtered(f13,'cav','<',140)
#
# plot_dists(terms,f7,'4rmb5bindingmetfilt.pdf')
#
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
print(len(f6))
print(len(f7))
print(len(f8))
print(len(f9))
print(len(f10))
print(len(f11))
print(len(f12))
print(len(f13))
print(len(f14))
'''
first go up to f7
19457
14850
3965
3003
2271
2247
1344
1247

second go looking for expression candidates, up to f14
11708
8851
4151
1827
1594
1587
1194
1116
628
541
248
246
236
201
189


batch4
11708
8851
4151
2337
1973
1959
1443
1442
804
688
570
565
532
483
478


b5
5422
2837
296
172
137
118
102
102
57
40
25
25
25
16
16

'''
import os
filtered_strc=[]
for d in f7:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered2',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i[:-5]+'.pdb filtered2/'+i+'.pdb')
    ########################################
os.system('mv *.pdf filtered2/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2 ~/desktop/esl4rmfilt


DIVERSITY AMONGST FILTERED DESIGNS

cd filtered2
'''
import os
des=[i for i in os.listdir() if i[-3:]=='pdb']
'''
'ed1UM_13_M99T39S17S95_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_8_unrelaxed_model_4_rank_1_0001_hybrid_39_1_0001_clean__DE_4_oglig_0001.pdb'
'''
motifs=[]
scaffs=[]
pairs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-9]
    scaffold=('_').join(s.split('_')[3:-9])
    motifs.append(motif)
    scaffs.append(scaffold)
    pairs.append((motif,scaffold))


# print(len(motifs))
# print(len(scaffs))
# print(len(pairs))
print(len(set(des)))
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(pairs)))
'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2 ~/desktop/esl4rmfilt
42
265
456

333
29
127
178

102
23
69
83

'''
same_scaff_diff_match={}
for scaff in set(scaffs):
    l=[]
    for d in des:
        s=d.split('.')[0]
        scaffold=('_').join(s.split('_')[3:-9])
        if scaffold==scaff:
            match=s.split('_')[2]
            if match not in l:
                l.append(match)
    same_scaff_diff_match[scaff]=l

for key in same_scaff_diff_match.keys():
    if len(same_scaff_diff_match[key])>1:
        print(key)
        print(same_scaff_diff_match[key])
'''

'''

























'''
SUBMIT TO COLABFOLD
/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2


echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate st
time python ~/BSFF/tools/run_colabfold.py input_dir output_dir
qstat -j "$JOB_ID"
'>run_cftest.sh
qsub -cwd -l mem_free=32G -o cluster_output -e cluster_output run_cftest.sh



echo '
# time python run_colabfold.py input_dir output_dir

import sys


input_dir = sys.argv[1]
result_dir = sys.argv[2]

# number of models to use
msa_mode = "single_sequence"
# msa_mode = "MMseqs2 (UniRef only)"
num_models = 5
num_recycles = 6
stop_at_score = 100
use_custom_msa = False
use_amber = False
use_templates = False
do_not_overwrite_results = True
zip_results = True



from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from pathlib import Path


setup_logging(Path(result_dir).joinpath("log.txt"))

queries, is_complex = get_queries(input_dir)


run(queries=queries,
    result_dir=result_dir,
    use_templates=use_templates,
    use_amber=use_amber,
    msa_mode=msa_mode,
    # model_type="auto",
    num_models=num_models,
    num_recycles=num_recycles,
    model_order=[3, 4, 5, 1, 2],
    is_complex=is_complex,
    data_dir=default_data_dir,
    keep_existing_results=do_not_overwrite_results,
    rank_mode="auto",
    pair_mode="unpaired+paired",
    stop_at_score=stop_at_score)
# zip_results=zip_results
'>run_colabfold_6recycles.py
'''
'''
AF2
'''
import os
os.makedirs('cf',exist_ok=True)
l=[i for i in os.listdir() if i[-3:]=='pdb']
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    newp=os.path.join('cf',i)
    of=open(newp,'w')
    for line in lines:
        of.write(line)
    of.close()

'''
cd cf
'''
os.chdir('cf')
l=[i for i in os.listdir() if i[-3:]=='pdb']
# In [2]: len(l)
# Out[2]: 125

'''
make fasta files
'''
from pyrosetta import *
init('-ignore_unrecognized_res')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())
    ofile=open(pdb[:-4]+'.fasta','w')
    ofile.write('>'+pdb+ '\n')
    length=len(sequence)
    remainder=length%60; n_seg=length/60
    indices=[]
    if length-remainder!=0:
        for i in range(0,60,length-remainder):
            indices.append(i)
    else:
        for i in range(0,60,length):
            indices.append(i)
    for i,x in enumerate(indices[:-1]):
        start=x
        end=indices[i+1]
        s=sequence[start:end]+'\n'
        ofile.write(s)
    last_lim=indices[-1]
    last_s=sequence[last_lim:length]
    ofile.write(str(last_s) +'\n')
    ofile.close()

for a in l:
    try:
        fasta(a)
    except:
        print('failure')

os.makedirs('fastas',exist_ok=True)
l2=[i for i in os.listdir() if i[-6:]=='.fasta']
for i in l2:
    os.system('mv '+i+' fastas/'+i)

'''
split fasta directory into subdirectories
cd fastas
'''
os.chdir('fastas')
# import os
nfastasperjob=1
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
l3=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if i[-6:]=='.fasta']

nfastas=len(l3)
directoriestosubmit=[]
c=1
for x in range(0,nfastas,nfastasperjob):
    if x+nfastasperjob<nfastas:
        fastaset=l3[x:x+nfastasperjob]
        currdirname=os.path.join(allfastasdir,('_').join([directory_prefix,str(c)]))
        directoriestosubmit.append(currdirname)
        os.makedirs(currdirname,exist_ok=True)
        for y in fastaset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1
    else:
        fastaset=l3[x:nfastas]
        currdirname=os.path.join(allfastasdir,('_').join([directory_prefix,str(c)]))
        directoriestosubmit.append(currdirname)
        os.makedirs(currdirname,exist_ok=True)
        for y in fastaset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1





# import os
#
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
shf_prefix='_cf'
#
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    submitlist.append((id,os.path.join(id,outputdirname)))
#
c=1
for id,od in submitlist:
    cmdl=['time', 'python3', '~/BSFF/tools/run_colabfold_6recycles.py', id, od]
    cmd=' '.join(cmdl)
    ofn=os.path.join(allfastasdir,'_'.join([shf_prefix,str(c),'.sh']))
    of=open(ofn,'w')
    of.write('#!/bin/bash')
    of.write('\n')
    of.write('source ~/anaconda3/etc/profile.d/conda.sh')
    of.write('\n')
    of.write('conda activate st')
    of.write('\n')
    of.write(cmd)
    of.write('\n')
    of.write('qstat -j "$JOB_ID"')
    of.close()
    c+=1

'''


cd cf/fastas
mkdir cluster_output


'''
# os.makedirs('cluster_output',exist_ok=True)
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('ntf2ogcf_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o cluster_output -e cluster_output '+j
    os.system(cmd)

'''
pwd
/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3/filtered/analysis/run/filtered2/cf/fastas


beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  730.33 GiB| 1000.00 GiB||  4300185|unlimited
'''


























'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='filtdesigns'
allfastasdir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3/filtered/analysis/run/filtered2/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3/filtered/analysis/run/filtered2'

#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
#
data_allstrc={}
seen_names=[]
ogds=[i for i in os.listdir(nonaf2desdir) if i[-3:]=='pdb']
for p in prediction_pdb_paths:
    # print(str(p))
    f=open(p,'r')
    fuzzball_lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[:6]=='HETATM']
    f.close()
    #first gotta get indices of unique residues in fuzzball_lines
    fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
    starts=[];lasts=[]
    for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
            try:
                resnum=int(fuzzball_lines[index][22:29].strip())
            except:
                print(index)
            # print((fuzzball_lines[index][22:29].strip()))
            resname=line[17:20]
            try:
                lastresnum=int(fuzzball_lines[index-1][22:29].strip())
                lastresname=fuzzball_lines[index-1][17:20]
                if resnum!=lastresnum or resname!=lastresname:
                    start=index
                    starts.append(start)
            except:
                start=index
                starts.append(start)
            try:
                nextresname=fuzzball_lines[index+1][17:20]
                next_resnum=int(fuzzball_lines[index+1][22:29].strip())
                if resnum!=next_resnum or resname!=nextresname:
                    last=index+1
                    lasts.append(last)
            except:
                last=len(fuzzball_lines)
                lasts.append(last)
    for index,start in enumerate(starts): #put the indices together for each res
        fuzzball_residue_indices.append((start,lasts[index]))
    fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
    fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
    #
    #
    #
    currstrc_lddts=[]
    for rin in fuzzball_residue_indices:
        rlddt=fuzzball_lines[rin[0]][60:66].strip(' ')
        currstrc_lddts.append(float(rlddt))
    currstrc_plddt=np.mean(currstrc_lddts)
    ############
    currstrc_name='_'.join(p.split('/')[-1].split('_')[:-5])###########################################
    csnogds=currstrc_name+'.pdb'
    if csnogds in ogds:
        p11=pose_from_pdb(p)
        p22=pose_from_pdb(os.path.join(nonaf2desdir,csnogds))
        carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
    ##########
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    #########
    currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
    ##################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
##################################################################
##################################################################
##################################################################
# import json
# with open('af2_data.json','r') as f:
#     data_allstrc=json.load(f)
# nonaf2desdir='/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2'

import numpy as np
#id designs with plddt over threshold
#use best vals:
plddt_threshold=88.0
carmsd_threshold=5.0
aplddts=[]
accepted={}
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        if cplddt>=plddt_threshold:
            if carmsd<=carmsd_threshold:
                if key not in list(accepted.keys()):
                    accepted[key]=[cplddt,carmsd]
                    plddts.append((cplddt,key,k2,carmsd))
                else:
                    if cplddt>accepted[key][0] and carmsd<accepted[key][1]:
                        accepted[key]=[cplddt,carmsd]
                        plddts.append((cplddt,key,k2,carmsd))
    try:
        aplddts.append(plddts[-1])
    except:
        pass
print(len(aplddts))

#plot the best vals
allbest={}
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        plddts.append(cplddt)
        carmsds.append(carmsd)
    bestp=max(plddts)
    bpi=plddts.index(bestp)
    bestc=carmsds[bpi]
    allbest[key]=[bestp,bestc]

import matplotlib.pyplot as plt
p=[]
c=[]
for key in allbest.keys():
    p.append(allbest[key][0])
    c.append(allbest[key][1])
fig = plt.figure()
ax = fig.add_subplot()
#
ax.scatter(c,p)
#
ax.set_title('CARMSD vs. Plddt (best of 5 models)')
ax.set_ylabel('Plddt')
ax.set_xlabel('CA RMSD')
plt.savefig('plddt_vs_carmsd.pdf')
plt.clf()
#

#use averages
# plddt_threshold=85.0
# carmsd_threshold=1.0
# aplddts=[]
# for key in data_allstrc.keys():
#     plddts=[]
#     carmsds=[]
#     for k2 in data_allstrc[key]:
#         cplddt=data_allstrc[key][k2][0]
#         carmsd=data_allstrc[key][k2][2]
#         plddts.append(cplddt)
#         carmsds.append(carmsd)
#     aplddt=np.mean(plddts)
#     acarmsd=np.mean(carmsds)
#     if aplddt>=plddt_threshold:
#         if acarmsd<=carmsd_threshold:
#             aplddts.append((aplddt,key,k2,acarmsd))


#move filtered to own directory
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
ns=[]
for i in aplddts:
    n='_'.join([i[1],'unrelaxed_model',i[2]])
    ns.append(n)
l=[i for i in os.listdir(nonaf2desdir) if i[-3:]=='pdb']
os.makedirs('af2filtered',exist_ok=True)
for i in aplddts:
    n=str(i[1])+'.pdb'
    if n in l:
        np=os.path.join(nonaf2desdir,n)
        newp=os.path.join('af2filtered',n)
        os.system('cp '+np+' '+newp)
os.system('mv *.pdf af2filtered')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed3/filtered/analysis/run/filtered2/cf/fastas/af2filtered ~/desktop/b5af2filt




IMA DO SOMETHING REAL QUICK OF COPYING ALL MATCHES WITHOUT QNR TO ANOTHER DIR

import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
good=[]
for i in l:
    motif=i.split('_')[2]
    chars=list(motif)
    if 'Q' in chars or 'N' in chars or 'R' in chars or 'H' in chars:
        pass
    else:
        if i not in good:
            good.append(i)

len(l)
len(good)


os.makedirs('esl4rm_1_bestmatches')
for i in good:
    os.system('cp '+i+' esl4rm_1_bestmatches/'+i)


pwd
/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res

ay you no what, gonna further id those which have at least one
y/w in the hbond motif

import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
good=[]
for i in l:
    motif=i.split('_')[2]
    chars=list(motif)
    if 'W' not in chars:
        if 'Y' not in chars:
            if 'F' not in chars:
                pass
    else:
        if i not in good:
            good.append(i)

len(good)
#4015

os.makedirs('has_aro_hbond')
for i in good:
    os.system('cp '+i+' has_aro_hbond/'+i)



/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/no_problematic_res/has_aro_hbond
'''



































'''
RUNNING MPNN ON MATCHES TO GENERATE SEQUENCE PROFILE


EXAMPLE OF FIXED POSITIONS DICT
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}

in this case i already applied dsasa filter and found buried matches
so i just need to identify the designable residues
'''
# prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
# ofjsonname='mpnn_params.json'
# ################
# import os
# from pyrosetta import*
# import json
# init('-load_PDB_components False')
# ##########
# # sf = ScoreFunction()
# # from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
# # sf.set_weight(fa_atr, 1)
# # sf.set_weight(fa_rep, .55)
# # sf.set_weight(fa_sol, 1)
# # sf.set_weight(hbond_sc, 1)
# # sf.set_weight(fa_elec, 1)
# # sf.set_weight(hbond_bb_sc,1)
# sf=get_fa_scorefxn()
#
# filters_xml = f'''
#                 <FILTERS>
#                     <DSasa name="dsasa" confidence="0.0"/>
#                 </FILTERS>'''
# dsasa_filter = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(filters_xml).get_filter("dsasa")
# #
# pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
# tfdata={}
#
# # ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'
# # ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn_params.json'
# #########
# # with open(ofjsonname,'r') as f:
# #     ttfdata=json.load(f)
# #
# # for pdb in ttfdata.keys():
# for pdb in pdbs:
#     lig=[prm1]
#     p=Pose()
#     generate_nonstandard_residue_set(p,lig)
#     pose_from_file(p, pdb)
#     p.update_residue_neighbors()
#     dsasa_val=dsasa_filter.report_sm(p)
#     if dsasa_val>0.8:
#         tofreeze=[]
#         #######################################################
#         #######################################################
#         #######################################################
#         #######################################################
#         #######################################################
#         #######################################################
#         sss=pdb.split('_')[2]
#         ss=''
#         for i in sss:
#             if i.isdigit()==True:
#                 ss+=str(i)
#             else:
#                 ss+='_'
#         ds=ss.split('_')
#         dss=[i for i in ds if i.isdigit()==True]
#         for i in dss:
#             tofreeze.append(int(i))
#         ####################
#         ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
#         neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
#         neighborhood_selector_bool = neighborhood_selector.apply(p)
#         neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
#         first_shell_res=list(neighborhood_residues_resnums)
#         ###############
#         for resnum in range(1,p.total_residue()):
#             if resnum not in first_shell_res:
#                 if resnum not in tofreeze:
#                     tofreeze.append(resnum)
#         tfdata[pdb]=tofreeze
#
# print(len(list(tfdata.keys())))
# json.dump(tfdata,open(ofjsonname,'w'))
################
################
################
################
################
################
################
import os
from pyrosetta import*
import json
init('-load_PDB_components False -gen_potential')
##########
allparams=[i for i in os.listdir() if i[-6:]=='params']
prm1=allparams[0]
ofjsonname='mpnn_params.json'
# sf = ScoreFunction()
# from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
# sf.set_weight(fa_atr, 1)
# sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
tfdata={}
for pdb in pdbs:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    tofreeze=[]
    #######################################################
    #######################################################
    #######################################################
    #######################################################
    #######################################################
    #######################################################
    sss=pdb.split('_')[2]
    ss=''
    for i in sss:
        if i.isdigit()==True:
            ss+=str(i)
        else:
            ss+='_'
    ds=ss.split('_')
    dss=[i for i in ds if i.isdigit()==True]
    for i in dss:
        tofreeze.append(int(i))
    ####################
    ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
    neighborhood_selector_bool = neighborhood_selector.apply(p)
    neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    first_shell_res=list(neighborhood_residues_resnums)
    ###############
    for resnum in range(1,p.total_residue()):
        if resnum not in first_shell_res:
            if resnum not in tofreeze:
                tofreeze.append(resnum)
    tfdata[pdb]=tofreeze

print(len(list(tfdata.keys())))
json.dump(tfdata,open(ofjsonname,'w'))
print('json output')
'''
MPNN WITH CONSTRAINTS
'''

#copy pdbs to new directory
#excluding the ligand
os.makedirs('mpnn',exist_ok=True)
l=[i for i in tfdata.keys()]
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    newp=os.path.join('mpnn',i)
    of=open(newp,'w')
    for line in lines:
        of.write(line)
    of.close()

os.chdir('mpnn')
#making the folders with the pdbs
#in this case i have to do 1 per job cus each file will
#have different fixed positions
# import json
#########
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'
# ofjsonname=''
# #########
# with open(ofjsonname,'r') as f:
#     tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/mpnn'
directory_prefix='mpnndesign'
###############
l=[i for i in os.listdir(all_strc_dir) if i[-3:]=='pdb']
nstrc=len(l)
c=1
for x in range(0,nstrc,n_strc_per_job):
    if x+n_strc_per_job<nstrc:
        strcset=l[x:x+n_strc_per_job]
        currdirname=os.path.join(all_strc_dir,('_').join([directory_prefix,str(c)]))
        os.makedirs(currdirname,exist_ok=True)
        for y in strcset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1
    else:
        strcset=l[x:nstrc]
        currdirname=os.path.join(all_strc_dir,('_').join([directory_prefix,str(c)]))
        os.makedirs(currdirname,exist_ok=True)
        for y in strcset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1


#making the fixed position dictionaries and putting them in the
#proper directories
#making the shellfile scripts
'''
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}
'''
import json
import os
#########
#########
ofjsonname='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/mpnn_params.json'
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/mpnn'
directory_prefix='mpnndesign'
shf_prefix='mpnn_'
###############
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    pdbid=[i for i in os.listdir(id) if i[-3:]=='pdb']
    pdbidd=pdbid[0].strip('.pdb')
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    p=[i for i in os.listdir(id) if i[-3:]=='pdb']
    pid=p[0]
    # tf=[str(i[0]) for i in tfdata[pid]]
    uol=[i for i in tfdata[pid]]
    ol=sorted(uol)
    odict={}
    sd={}
    sd["A"]=ol
    odict[pdbidd]=sd
    with open(os.path.join(id,outputdirname,'fixed_pdbs.jsonl'),'w') as of:
        json.dump(odict,of)
    submitlist.append((id,os.path.join(id,outputdirname)))
#
print(len(submitlist))
c=1
for id,od in submitlist:
    ofn=os.path.join(all_strc_dir,'_'.join([shf_prefix,str(c),'.sh']))
    of=open(ofn,'w')
    of.write('#!/bin/bash')
    of.write('\n')
    of.write('source /wynton/home/kortemme/cgalvin/mpnn_test_env/bin/activate')
    of.write('\n')
    of.write('folder_with_pdbs="'+id+'"')
    of.write('\n')
    of.write('output_dir="'+od+'"')
    of.write('\n')
    of.write('if [ ! -d $output_dir ]')
    of.write('\n')
    of.write('then')
    of.write('\n')
    of.write('    mkdir -p $output_dir')
    of.write('\n')
    of.write('fi')
    of.write('\n')
    of.write('path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"')
    of.write('\n')
    of.write('path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"')
    of.write('\n')
    of.write('python /wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains')
    of.write('\n')
    cmdl=['python',
    '/wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/protein_mpnn_run.py',
    '--jsonl_path $path_for_parsed_chains',
    '--out_folder $output_dir',
    '--fixed_positions_jsonl $path_for_fixed_positions',
    '--num_seq_per_target 500',
    '--sampling_temp "0.15"',
    '--batch_size 1']
    cmd=' '.join(cmdl)
    of.write(cmd)
    of.write('\n')
    of.write('qstat -j "$JOB_ID"')
    of.close()
    c+=1

'''
mkdir mpnncout
qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout mpnn__1_.sh

'''
import os
shf_prefix='mpnn_'
jobss=[i for i in os.listdir() if shf_prefix in i and i[-2:]=='sh']
#
jobss.remove('mpnn__1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)

'''
~30 minutes per job only












consolidating mpnn results
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/mpnn'
directory_prefix='mpnndesign'
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
for id in input_directories:
    resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
    if os.path.exists(resultfspath):
        resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
        for y in resl:
            allresl.append(os.path.join(resultfspath,y))
    else:
        print(resultfspath)
'''
print(len(allresl))
1939
analyze the mpnn results
'''
#
add={}
og_data=[]
#
for i in range(len(allresl)):
    design_data=[]
    f=open(allresl[i],'r')
    lines=[line for line in f.readlines()]
    f.close()
    #
    if len(lines)>0:
        ogname=lines[0].split(',')[0].strip(' ').strip('>')
        ogscore=float(lines[0].split(',')[1].strip(' ').strip('score='))
        og_data.append((ogname,ogscore))
        scores=[]
        for i in range(2,len(lines),2):
            score=float(lines[i].split(',')[2].strip(' ').strip('score='))
            seqrec=float(lines[i].split(',')[3].strip(' ').strip('seq_recovery='))
            seq=lines[i+1].strip('\n')
            design_data.append((score,seqrec,seq))
            scores.append(score)
        design_data=sorted(design_data, key=lambda first: first[0])
        add[ogname]=design_data

'''
add[name]=score,seqrec,seq
'''
import os
from pyrosetta import*
import json
init('-gen_potential -load_PDB_components False')

#########################
ofjsonname='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/mpnn_params.json'
pdbsdir='/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean'
allparams=[os.path.join(pdbsdir,i) for i in os.listdir(pdbsdir) if i[-6:]=='params']
prm1=allparams[0]
#########################

with open(ofjsonname,'r') as f:
    tfdata=json.load(f)


os.makedirs('resfiles',exist_ok=True)
seqs_data={}
for strc in list(add.keys()):
    pdb=os.path.join(pdbsdir,strc+'.pdb')
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    ###################
    frozen_res=tfdata[strc+'.pdb']
    des_res=[]
    for i in range(1,p.total_residue()):
        if i not in frozen_res:
            des_res.append(i)
    ###############
    strc_seqs={}
    for des_ind in des_res:
        positional_seq_list=[]
        for mpnn_data_entry in add[strc]:
            score=mpnn_data_entry[0]
            if score<=1.0:
                seq=mpnn_data_entry[2]
                positional_seq_list.append(seq[des_ind-1])
        strc_seqs[des_ind]=positional_seq_list
    seqs_data[strc]=strc_seqs
    #
    p.update_residue_neighbors()
    ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
    neighborhood_selector_bool = neighborhood_selector.apply(p)
    neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    first_shell_res=list(neighborhood_residues_resnums)
    #
    second_shell=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(neighborhood_selector,10.,False)
    ss_bool = second_shell.apply(p)
    ss_resnums=pyrosetta.rosetta.core.select.get_residues_from_subset(ss_bool)
    ss_res=list(ss_resnums)
    #
    for i in des_res:-;
        if i in ss_res:
            ss_res.remove(i)
    #
    ofilename=os.path.join('resfiles',strc+'.resfile')
    ofile=open(ofilename,'w')
    #write header
    ofile.write('USE_INPUT_SC')
    ofile.write('\nstart'+'\n')
    #PIKAA
    for i in range(1,p.total_residue()+1):
        if i in des_res:
            allowableres=list(set(seqs_data[strc][i]))
            s=str(p.pdb_info().pose2pdb(i))+'PIKAA '+('').join(allowableres)
            ofile.write(s+'\n')
        elif i in ss_res:
            s=str(p.pdb_info().pose2pdb(i))+'NATAA'
            ofile.write(s+'\n')
        elif i==p.total_residue():
            s=str(p.pdb_info().pose2pdb(i))+'NATAA'
            ofile.write(s+'\n')
        else:
            s=str(p.pdb_info().pose2pdb(i))+'NATRO'
            ofile.write(s+'\n')
    ofile.close()

'''
ALRIGHT NOW I THINK IVE GOT THIS WHOLE WORKFLOW IN ORDER AND WORKING
but first thought is that there arent very many residues allowed at design positions,
perhaps i shoyld try to increase the temperature for mpnn to get some more diversity?
i will hold this thought in my head, but for now i will indeed go for it with
what i have and see what happens

AAHHHHH IM AN IDIOT THOUGH,
I ACTUALLY DONT WANNA DO THIS NO,
I WANNA DO THIS AFTER DESIGNING FOR GOOD LIGAND INTERACTIONS FIRST
AND THEN FREEZE THOSE INTERACTING RESIDUES ALSO TO TRY TO STABILIZE THE DESIGN
    --->look at scores/af2 conf of fd designs
            redesign for stability with mpnn assist
                    look at scores/af2 conf of resultant des, compare with og fd designs
IN THE MEANTIME SET THE 3BOP FD TO RUN


'''




































'''
RUNNING MPNN ON af2filt designs TO GENERATE SEQUENCE PROFILE

freezing residues that hbond w/ lig, or other residues that those res hb with
(bb-bb not included),
as well as res that have <= -2.0 REU pairwise interaxn w lig
'''
################
################
################
################
################
import os
from pyrosetta import*
import json
init('-load_PDB_components False')
##########
paramsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run'
allparams=[i for i in os.listdir(paramsdir) if i[-6:]=='params']
prm1=os.path.join(paramsdir,allparams[0])
ofjsonname='mpnn_params.json'
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(hbond_sc, 1)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,1)
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
tfdata={}
for pdb in pdbs:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    tofreeze=[]
    ####################
    ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
    neighborhood_selector_bool = neighborhood_selector.apply(p)
    neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    first_shell_res=list(neighborhood_residues_resnums)
    ligand_pose=p.residue(p.total_residue()).clone()
    ###############
    for fsr in first_shell_res:
        network_pose=Pose()
        network_pose.append_residue_by_jump(ligand_pose, 1)
        res_pose1=p.residue(fsr).clone()
        network_pose.append_residue_by_jump(res_pose1, 1)
        pairwise_score=sf(network_pose)
        # print(pairwise_score)
        if pairwise_score<=-2.0:
            tofreeze.append(fsr)
    ############
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    network_pose.update_residue_neighbors()
    rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set,exclude_bb=True)
    s=p.sequence()
    hbond_data={}
    if hbond_set.nhbonds()>0:
        for hbond_index in range(1,hbond_set.nhbonds()+1):
            donres_ind=int(hbond_set.hbond(hbond_index).don_res())
            accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
            donres=s[donres_ind-1]
            accres=s[accres_ind-1]
            acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
            donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
            don_atom_index=int(p.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
            acc_atom_data=str(p.residue(accres_ind).atom_type(acc_atom_index))
            acc_atom=(acc_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
            don_atom_data=str(p.residue(donres_ind).atom_type(don_atom_index))
            don_atom=(don_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
            acc_atom_name=str(p.residue(accres_ind).atom_name(acc_atom_index)).strip(' ')
            don_atom_name=str(p.residue(donres_ind).atom_name(don_atom_index)).strip(' ')
            hbond_data[hbond_index]=[donres_ind,accres_ind,donres,accres,don_atom,acc_atom,don_atom_name,acc_atom_name]
    #identify the network of hbonds about the ligand
    network_members=[]
    hb_w_lig=[]
    ptn_lig_hb=[]
    for keyb in hbond_data.keys():
        l=hbond_data[keyb]
        if l[0]==p.total_residue():
            hb_w_lig.append(l[1])
            network_members.append((l[0],l[1],l[-2],l[-1]))
            ptn_lig_hb.append((l[-2],l[-3],'ligdon'))
        if l[1]==p.total_residue():
            hb_w_lig.append(l[0])
            network_members.append((l[0],l[1],l[-2],l[-1]))
            ptn_lig_hb.append((l[-1],l[-4],'ligacc'))
    hbws=[]
    for i in hb_w_lig:
        for keyc in hbond_data.keys():
            l=hbond_data[keyc]
            if l[0]==i:
                if l[1]!=p.total_residue():
                    hbws.append(l[1])
                    network_members.append((l[0],l[1],l[-2],l[-1]))
            if l[1]==i:
                if l[0]!=p.total_residue():
                    hbws.append(l[0])
                    network_members.append((l[0],l[1],l[-2],l[-1]))
    for rn in hb_w_lig:
        if rn not in tofreeze:
            tofreeze.append(rn)
    for rn in hbws:
        if rn not in tofreeze:
            tofreeze.append(rn)
    tfdata[pdb]=tofreeze


print(len(list(tfdata.keys())))
json.dump(tfdata,open(ofjsonname,'w'))
print('json output')






'''
MPNN WITH CONSTRAINTS
'''

#copy pdbs to new directory
#excluding the ligand
os.makedirs('mpnn',exist_ok=True)
l=[i for i in tfdata.keys()]
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    newp=os.path.join('mpnn',i)
    of=open(newp,'w')
    for line in lines:
        of.write(line)
    of.close()

os.chdir('mpnn')
#making the folders with the pdbs
#in this case i have to do 1 per job cus each file will
#have different fixed positions
# import json
#########
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'
# ofjsonname=''
# #########
# with open(ofjsonname,'r') as f:
#     tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
directory_prefix='mpnndesign'
###############
l=[i for i in os.listdir(all_strc_dir) if i[-3:]=='pdb']
nstrc=len(l)
c=1
for x in range(0,nstrc,n_strc_per_job):
    if x+n_strc_per_job<nstrc:
        strcset=l[x:x+n_strc_per_job]
        currdirname=os.path.join(all_strc_dir,('_').join([directory_prefix,str(c)]))
        os.makedirs(currdirname,exist_ok=True)
        for y in strcset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1
    else:
        strcset=l[x:nstrc]
        currdirname=os.path.join(all_strc_dir,('_').join([directory_prefix,str(c)]))
        os.makedirs(currdirname,exist_ok=True)
        for y in strcset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1


#making the fixed position dictionaries and putting them in the
#proper directories
#making the shellfile scripts
'''
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}
'''
import json
import os
#########
#########
ofjsonname='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn_params.json'
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
directory_prefix='mpnndesign'
shf_prefix='mpnn_'
###############
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    pdbid=[i for i in os.listdir(id) if i[-3:]=='pdb']
    pdbidd=pdbid[0].strip('.pdb')
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    p=[i for i in os.listdir(id) if i[-3:]=='pdb']
    pid=p[0]
    # tf=[str(i[0]) for i in tfdata[pid]]
    uol=[i for i in tfdata[pid]]
    ol=sorted(uol)
    odict={}
    sd={}
    sd["A"]=ol
    odict[pdbidd]=sd
    with open(os.path.join(id,outputdirname,'fixed_pdbs.jsonl'),'w') as of:
        json.dump(odict,of)
    submitlist.append((id,os.path.join(id,outputdirname)))
#
print(len(submitlist))
c=1
for id,od in submitlist:
    ofn=os.path.join(all_strc_dir,'_'.join([shf_prefix,str(c),'.sh']))
    of=open(ofn,'w')
    of.write('#!/bin/bash')
    of.write('\n')
    of.write('source /wynton/home/kortemme/cgalvin/mpnn_test_env/bin/activate')
    of.write('\n')
    of.write('folder_with_pdbs="'+id+'"')
    of.write('\n')
    of.write('output_dir="'+od+'"')
    of.write('\n')
    of.write('if [ ! -d $output_dir ]')
    of.write('\n')
    of.write('then')
    of.write('\n')
    of.write('    mkdir -p $output_dir')
    of.write('\n')
    of.write('fi')
    of.write('\n')
    of.write('path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"')
    of.write('\n')
    of.write('path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"')
    of.write('\n')
    of.write('python /wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains')
    of.write('\n')
    cmdl=['python',
    '/wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/protein_mpnn_run.py',
    '--jsonl_path $path_for_parsed_chains',
    '--out_folder $output_dir',
    '--fixed_positions_jsonl $path_for_fixed_positions',
    '--num_seq_per_target 500',
    '--sampling_temp "0.2"',
    '--batch_size 1']
    cmd=' '.join(cmdl)
    of.write(cmd)
    of.write('\n')
    of.write('qstat -j "$JOB_ID"')
    of.close()
    c+=1

'''
mkdir mpnncout
qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout mpnn__1_.sh

'''
import os
shf_prefix='mpnn_'
jobss=[i for i in os.listdir() if shf_prefix in i and i[-2:]=='sh']
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)

'''
~30 minutes per job only


/bin/sh: git: command not found









consolidating mpnn results
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
directory_prefix='mpnndesign'
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
for id in input_directories:
    resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
    if os.path.exists(resultfspath):
        resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
        for y in resl:
            allresl.append(os.path.join(resultfspath,y))
    else:
        print(resultfspath)
'''
print(len(allresl))
1939
analyze the mpnn results
'''
#
add={}
og_data=[]
#
for i in range(len(allresl)):
    design_data=[]
    f=open(allresl[i],'r')
    lines=[line for line in f.readlines()]
    f.close()
    #
    if len(lines)>0:
        ogname=lines[0].split(',')[0].strip(' ').strip('>')
        ogscore=float(lines[0].split(',')[1].strip(' ').strip('score='))
        og_data.append((ogname,ogscore))
        scores=[]
        for i in range(2,len(lines),2):
            score=float(lines[i].split(',')[2].strip(' ').strip('score='))
            seqrec=float(lines[i].split(',')[3].strip(' ').strip('seq_recovery='))
            seq=lines[i+1].strip('\n')
            design_data.append((score,seqrec,seq))
            scores.append(score)
        design_data=sorted(design_data, key=lambda first: first[0])
        add[ogname]=design_data

'''
add[name]=score,seqrec,seq
'''
import os
from pyrosetta import*
import json
init('-load_PDB_components False')

#########################
ofjsonname='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn_params.json'
pdbsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/'
paramsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/'
allparams=[os.path.join(paramsdir,i) for i in os.listdir(paramsdir) if i[-6:]=='params']
prm1=allparams[0]
#########################

with open(ofjsonname,'r') as f:
    tfdata=json.load(f)



amino_acids={'ALA':'A',
             'ARG':'R',
             'ASN':'N',
             'ASP':'D',
             'CYS':'C',
             'GLU':'E',
             'GLN':'Q',
             'GLY':'G',
             'HIS':'H',
             'ILE':'I',
             'LEU':'L',
             'LYS':'K',
             'MET':'M',
             'PHE':'F',
             'PRO':'P',
             'SER':'S',
             'THR':'T',
             'TRP':'W',
             'TYR':'Y',
             'VAL':'V'}


os.makedirs('resfiles',exist_ok=True)
seqs_data={}
for strc in list(add.keys()):
    pdb=os.path.join(pdbsdir,strc+'.pdb')
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    ###################
    frozen_res=tfdata[strc+'.pdb']
    des_res=[]
    for i in range(1,p.total_residue()):
        if i not in frozen_res:
            des_res.append(i)
    ###############
    strc_seqs={}
    for des_ind in des_res:
        positional_seq_list=[]
        for mpnn_data_entry in add[strc]:
            score=mpnn_data_entry[0]
            if score<=1.0:
                seq=mpnn_data_entry[2]
                positional_seq_list.append(seq[des_ind-1])
        strc_seqs[des_ind]=positional_seq_list
    seqs_data[strc]=strc_seqs
    #
    ofilename=os.path.join('resfiles',strc+'.resfile')
    ofile=open(ofilename,'w')
    #write header
    ofile.write('USE_INPUT_SC')
    ofile.write('\nstart'+'\n')
    #PIKAA
    for i in range(1,p.total_residue()+1):
        if i in strc_seqs.keys():
            allowableres=list(set(seqs_data[strc][i]))
            current_res=p.residue(i).name()[:3]
            olccr=amino_acids[current_res]
            if olccr not in allowableres:
                allowableres.append(olccr)
            s=str(p.pdb_info().pose2pdb(i))+'PIKAA '+('').join(allowableres)
            ofile.write(s+'\n')
        else:
            s=str(p.pdb_info().pose2pdb(i))+'NATAA'
            ofile.write(s+'\n')
    ofile.close()
# print(strc_seqs.keys())
# for key in strc_seqs.keys():
#     print(len(set(strc_seqs[key])))
'''
THERES A QUESTION HERE IF I WANNA PROCEED WITH THIS ALL SEEN RES ALLOWED SCHEME
OR LIMIT TO THOSE SEEN ABOVE A CERTAIN FREQUENCY
        idk maybe for now i just roll with it
okay next is to get pdbs and params in this same dir as resfiles w
proper naming scheme
then run multiple cycles of fd 3bop, followed by filters

'''
import os
#########
pdbsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/'
paramsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/'
allparams=[os.path.join(paramsdir,i) for i in os.listdir(paramsdir) if i[-6:]=='params']
allpdbs=[os.path.join(pdbsdir,i) for i in os.listdir(pdbsdir) if i[-3:]=='pdb']
allresfiles=[i for i in os.listdir() if i[-7:]=='resfile']
###########
spdbs=[]
for a in allresfiles:
    s=a.split('.')[0]
    for i in allpdbs:
        if s in i:
            os.system('cp '+i+' '+os.getcwd())
            spdbs.append(i)
sparams=[]
for a in allresfiles:
    s=a.split('.')[0]
    for i in allparams:
        if s.strip('_oglig_0001.resfile') in i:
            ###########################################################################################
            os.system('cp '+i+' '+os.getcwd()+'/'+i.split('.')[0].split('/')[-1]+'_oglig_0001.params')
            #####################
            sparams.append(i)
print(len(allresfiles))
print(len(spdbs))
print(len(sparams))
'''
okay having an issue here where there are too many params file,
i guess ill remove the extraneoius ones but wtf is going on with this?
                SKETCHY
'''
##
import os
allparams=[i for i in os.listdir() if i[-6:]=='params']
allpdbs=[i for i in os.listdir() if i[-3:]=='pdb']
allresfiles=[i for i in os.listdir() if i[-7:]=='resfile']
print(len(allresfiles))
print(len(allpdbs))
print(len(allparams))
for i in allparams:
    if not os.path.exists(i.strip('.params')+'.pdb'):
        os.system('rm '+i)
allparams=[i for i in os.listdir() if i[-6:]=='params']
allpdbs=[i for i in os.listdir() if i[-3:]=='pdb']
allresfiles=[i for i in os.listdir() if i[-7:]=='resfile']
print(len(allresfiles))
print(len(allpdbs))
print(len(allparams))






'''
FASTDESIGN with 3bop weight 5
no intup
some filters but no ddg/boltz
'''
# FASTDESIGN with 3bop weight 5, intup weight 2
import os
########
allparams=[i for i in os.listdir() if i[-6:]=='params']
prm1=allparams[0]
sfname='fd3bopmpnn.sh'
ndesignseachmatch='5'
outputdirectory='fd_3bop_mpnn'
scorefile_name='esl_fd_3bop_mpnn.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(20):     #this way i can spread across more nodes if neccessary
    output_prefix='fd'+str(c)
    sf=open('_'+str(c)+sfname,'w')
    sf.write('#!/bin/bash')
    sf.write('\n')
    sf.write('tasks=(0\n')
    for match in matches[:-1]:
        sf.write('       '+str(match.strip('.pdb'))+'\n')
    sf.write('       '+matches[-1].strip('.pdb')+')')
    sf.write('\n')
    cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
         '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
         '-parser:protocol','~/BSFF/tools/fd3bop_mpnn.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0','-use_input_sc',
         '-out:nstruct',ndesignseachmatch,
         # '-score:set_weights','hbond_bb_sc','2.0',
         # '-score:set_weights','hbond_sc','2.0',
         '-load_PDB_components','False',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name]
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1
print(len(matches))
'''
print(len(matches))
155

~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s ed9UM_1_L53F106Y57S101S35_1_relaxed_relaxed_5tpj_67405_design_4_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_hybrid_105_1_clean__DE_8_0001_oglig_0001.pdb -parser:protocol ~/BSFF/tools/fd3bop_mpnn.xml -parser:view -run:preserve_header -in:file::extra_res_fa ed13UM_1_L51M99Y18T108S48_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_193_1_clean__DE_2_0001_oglig_0001.params -resfile ed9UM_1_L53F106Y57S101S35_1_relaxed_relaxed_5tpj_67405_design_4_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_hybrid_105_1_clean__DE_8_0001_oglig_0001.resfile -ex1 -ex2 -extrachi_cutoff 0 -use_input_sc -out:nstruct 1 -load_PDB_components False -out:path:all fd_3bop_mpnn -out:prefix fd1 -scorefile_format json -out:file:scorefile esl_fd_3bop_mpnn.json

'''
import os
sfname='fd3bopmpnn.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-155 -l mem_free=4G -o cluster_out -e cluster_out '+j
    os.system(cmd)



'''
only took like 4 hours to do 5desx20 jobs on 155 structures

'''




'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE

/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn

'''




#analysis of scores from json file
sfname='esl_fd_3bop_mpnn.json'


import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==2:
        l1=line[:ids[0]+1]
        l2=line[ids[0]+1:]
        scores.append(l1)
        scores.append(l2)
terms=list(scores[0].keys())
print(len(bad))
#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                try:
                    allscores.append(float(d[term]))
                except:
                    print(d)
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

plot_dists(terms,scores,'esl_fd3bmpnn.pdf')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'lighphobesasa','<',20.0)
f3=return_filtered(f2,'hbtolig','>',3.0)
f4=return_filtered(f3,'contact_molsurf','>',160.0)
f5=return_filtered(f4,'buns_bb_heavy','<',4.0)
f6=return_filtered(f5,'packstat','>',0.65)
f7=return_filtered(f6,'oversat','<',0.0)
f8=return_filtered(f7,'ligoversat','<',0.0)
f9=return_filtered(f8,'buns_sc_heavy','<',0.0)
f10=return_filtered(f9,'exphyd','<',850.0)
f11=return_filtered(f9,'shape_comp','>',0.65)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
print(len(f6))
print(len(f7))
print(len(f8))
print(len(f9))
print(len(f10))
print(len(f11))
'''
15494
11937
11768
11208
10322
6895
1679
1670
1670
974
965
959
'''
import os
filtered_strc=[]
for d in f11:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered_fd3bop2',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered_fd3bop2/'+i+'.pdb')
os.system('mv *.pdf filtered_fd3bop2/')

'''
all
21
12
18
9

filtered
10
8
10
8
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop ~/desktop/eslfd3bmpnnfilt
'''





























#run the extra scoring
import os
paramsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles'
params=[os.path.join(paramsdir,i) for i in os.listdir(paramsdir) if i[-6:]=='params']
paramspath=params[0]

pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('extrascore.sh','w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('tasks=(0\n')
for match in pdbs[:-1]:
    sf.write('       '+match+'\n')
sf.write('       '+pdbs[-1]+')')
sf.write('\n')
sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]} -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+paramspath+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
#
print(len(pdbs))
'''
qsub -cwd -t 1-959 -l mem_free=2G extrascore.sh
        took about an hour per job

/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop
'''




#analysis of scores from json file
sfname='scores.json'
import json
bad=[]
scores=[]
for line in open(sfname,'r'):
    try:
        scores.append(json.loads(line))
    except:
        bad.append(line)
for line in bad:
    ids=[]
    for ii,i in enumerate(line):
        if i=='}':
            ids.append(ii)
    if len(ids)==2:
        l1=line[:ids[0]+1]
        l2=line[ids[0]+1:]
        scores.append(l1)
        scores.append(l2)
terms=list(scores[0].keys())
print(len(bad))
#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                try:
                    allscores.append(float(d[term]))
                except:
                    print(d)
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

plot_dists(terms,scores,'esled2_fd3m_extrascores.pdf')
'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/dog_idel_fd_filt_extrascores.pdf ~/desktop/dog_idel_fd_filt_extrascores.pdf
'''
def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'boltz','<',-0.33)
f3=return_filtered(f2,'ddg','<',-25.0)
f4=return_filtered(f3,'contact_molsurf','>',160)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.66)
f7=return_filtered(f6,'lighphobesasa','<',25.0)
f8=return_filtered(f7,'buns_bb_heavy','<',4.0)
f9=return_filtered(f8,'packstat','>',0.65)
f10=return_filtered(f9,'oversat','<',0.0)
f11=return_filtered(f10,'ligoversat','<',0.0)
f12=return_filtered(f11,'buns_sc_heavy','<',0.0)
f13=return_filtered(f12,'exphyd','<',830.0)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
print(len(f6))
print(len(f7))
print(len(f8))
print(len(f9))
print(len(f10))
print(len(f11))
print(len(f12))
print(len(f13))
'''
953
927
921
812
811
805
789
789
787
495
495
495
492
489
'''
import os
filtered_strc=[]
for d in f13:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered_extrascores2',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i[:-5]+'.pdb filtered_extrascores2/'+i+'.pdb')
    ########################################
os.system('mv *.pdf filtered_extrascores2/')
'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop2/filtered_extrascores2 ~/desktop/esl_fd3bmpnn_extrascores_filt
print(len(set(pairs)))
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(ogscaffs)))
'''

























'''
SUBMIT TO COLABFOLD

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate st
time python ~/BSFF/tools/run_colabfold.py input_dir output_dir
qstat -j "$JOB_ID"
'>run_cftest.sh
qsub -cwd -l mem_free=32G -o cluster_output -e cluster_output run_cftest.sh

'''
'''
AF2
'''
import os
os.makedirs('cf',exist_ok=True)
l=[i for i in os.listdir() if i[-3:]=='pdb']
for i in l:
    f=open(i,'r')
    lines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    newp=os.path.join('cf',i)
    of=open(newp,'w')
    for line in lines:
        of.write(line)
    of.close()

'''
cd cf
'''
os.chdir('cf')
l=[i for i in os.listdir() if i[-3:]=='pdb']
# In [2]: len(l)
# Out[2]: 125

'''
make fasta files
'''
from pyrosetta import *
init('-ignore_unrecognized_res')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())
    ofile=open(pdb[:-4]+'.fasta','w')
    ofile.write('>'+pdb+ '\n')
    length=len(sequence)
    remainder=length%60; n_seg=length/60
    indices=[]
    if length-remainder!=0:
        for i in range(0,60,length-remainder):
            indices.append(i)
    else:
        for i in range(0,60,length):
            indices.append(i)
    for i,x in enumerate(indices[:-1]):
        start=x
        end=indices[i+1]
        s=sequence[start:end]+'\n'
        ofile.write(s)
    last_lim=indices[-1]
    last_s=sequence[last_lim:length]
    ofile.write(str(last_s) +'\n')
    ofile.close()

for a in l:
    try:
        fasta(a)
    except:
        print('failure')

os.makedirs('fastas',exist_ok=True)
l2=[i for i in os.listdir() if i[-6:]=='.fasta']
for i in l2:
    os.system('mv '+i+' fastas/'+i)

'''
split fasta directory into subdirectories
cd fastas
'''
os.chdir('fastas')
# import os
nfastasperjob=1
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
l3=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if i[-6:]=='.fasta']

nfastas=len(l3)
directoriestosubmit=[]
c=1
for x in range(0,nfastas,nfastasperjob):
    if x+nfastasperjob<nfastas:
        fastaset=l3[x:x+nfastasperjob]
        currdirname=os.path.join(allfastasdir,('_').join([directory_prefix,str(c)]))
        directoriestosubmit.append(currdirname)
        os.makedirs(currdirname,exist_ok=True)
        for y in fastaset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1
    else:
        fastaset=l3[x:nfastas]
        currdirname=os.path.join(allfastasdir,('_').join([directory_prefix,str(c)]))
        directoriestosubmit.append(currdirname)
        os.makedirs(currdirname,exist_ok=True)
        for y in fastaset:
            os.system('mv '+y+' '+currdirname+'/'+y.split('/')[-1])
        c+=1





# import os
#
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
shf_prefix='_cf'
#
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    submitlist.append((id,os.path.join(id,outputdirname)))
#
c=1
for id,od in submitlist:
    cmdl=['time', 'python3', '~/BSFF/tools/run_colabfold.py', id, od]
    cmd=' '.join(cmdl)
    ofn=os.path.join(allfastasdir,'_'.join([shf_prefix,str(c),'.sh']))
    of=open(ofn,'w')
    of.write('#!/bin/bash')
    of.write('\n')
    of.write('source ~/anaconda3/etc/profile.d/conda.sh')
    of.write('\n')
    of.write('conda activate st')
    of.write('\n')
    of.write(cmd)
    of.write('\n')
    of.write('qstat -j "$JOB_ID"')
    of.close()
    c+=1

'''
mkdir cluster_output
qsub -cwd -l mem_free=32G -o cluster_output -e cluster_output ntf2ogcf_1_.sh

            okay appears to be working, i will submit all the rest
            manually removing the test one from the jobss list below

cd cf/fastas
mkdir cluster_output

'''
os.makedirs('cluster_output',exist_ok=True)
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('ntf2ogcf_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o cluster_output -e cluster_output '+j
    os.system(cmd)

































'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='filtdesigns'
allfastasdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop2/filtered_extrascores2/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop2/filtered_extrascores2'

#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
#
data_allstrc={}
seen_names=[]
ogds=[i for i in os.listdir(nonaf2desdir) if i[-3:]=='pdb']
for p in prediction_pdb_paths:
    # print(str(p))
    f=open(p,'r')
    fuzzball_lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[:6]=='HETATM']
    f.close()
    #first gotta get indices of unique residues in fuzzball_lines
    fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
    starts=[];lasts=[]
    for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
            try:
                resnum=int(fuzzball_lines[index][22:29].strip())
            except:
                print(index)
            # print((fuzzball_lines[index][22:29].strip()))
            resname=line[17:20]
            try:
                lastresnum=int(fuzzball_lines[index-1][22:29].strip())
                lastresname=fuzzball_lines[index-1][17:20]
                if resnum!=lastresnum or resname!=lastresname:
                    start=index
                    starts.append(start)
            except:
                start=index
                starts.append(start)
            try:
                nextresname=fuzzball_lines[index+1][17:20]
                next_resnum=int(fuzzball_lines[index+1][22:29].strip())
                if resnum!=next_resnum or resname!=nextresname:
                    last=index+1
                    lasts.append(last)
            except:
                last=len(fuzzball_lines)
                lasts.append(last)
    for index,start in enumerate(starts): #put the indices together for each res
        fuzzball_residue_indices.append((start,lasts[index]))
    fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
    fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
    #
    #
    #
    currstrc_lddts=[]
    for rin in fuzzball_residue_indices:
        rlddt=fuzzball_lines[rin[0]][60:66].strip(' ')
        currstrc_lddts.append(float(rlddt))
    currstrc_plddt=np.mean(currstrc_lddts)
    ############
    currstrc_name='_'.join(p.split('/')[-1].split('_')[:-5])###########################################
    csnogds=currstrc_name+'.pdb'
    if csnogds in ogds:
        p11=pose_from_pdb(p)
        p22=pose_from_pdb(os.path.join(nonaf2desdir,csnogds))
        carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
    ##########
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    #########
    currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
    ##################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
##################################################################
##################################################################
##################################################################

import json
with open('af2_data.json','r') as f:
    data_allstrc=json.load(f)

import numpy as np
#id designs with plddt over threshold
#use best vals:
plddt_threshold=85.0
carmsd_threshold=1.5
aplddts=[]
accepted={}
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        if cplddt>=plddt_threshold:
            if carmsd<=carmsd_threshold:
                if key not in list(accepted.keys()):
                    accepted[key]=[cplddt,carmsd]
                    plddts.append((cplddt,key,k2,carmsd))
                else:
                    if cplddt>accepted[key][0] and carmsd<accepted[key][1]:
                        accepted[key]=[cplddt,carmsd]
                        plddts.append((cplddt,key,k2,carmsd))
    try:
        aplddts.append(plddts[-1])
    except:
        pass

print(len(aplddts))



#plot the best vals
allbest={}
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        plddts.append(cplddt)
        carmsds.append(carmsd)
    bestp=max(plddts)
    bpi=plddts.index(bestp)
    bestc=carmsds[bpi]
    allbest[key]=[bestp,bestc]

import matplotlib.pyplot as plt
p=[]
c=[]
for key in allbest.keys():
    p.append(allbest[key][0])
    c.append(allbest[key][1])
fig = plt.figure()
ax = fig.add_subplot()
#
ax.scatter(c,p)
#
ax.set_title('CARMSD vs. Plddt (best of 5 models)')
ax.set_ylabel('Plddt')
ax.set_xlabel('CA RMSD')
plt.savefig('plddt_vs_carmsd.pdf')
plt.clf()
#
'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas/plddt_vs_carmsd.pdf ~/desktop/plddt_vs_carmsd.pdf
'''
#use averages
# plddt_threshold=85.0
# carmsd_threshold=1.0
# aplddts=[]
# for key in data_allstrc.keys():
#     plddts=[]
#     carmsds=[]
#     for k2 in data_allstrc[key]:
#         cplddt=data_allstrc[key][k2][0]
#         carmsd=data_allstrc[key][k2][2]
#         plddts.append(cplddt)
#         carmsds.append(carmsd)
#     aplddt=np.mean(plddts)
#     acarmsd=np.mean(carmsds)
#     if aplddt>=plddt_threshold:
#         if acarmsd<=carmsd_threshold:
#             aplddts.append((aplddt,key,k2,acarmsd))


#move filtered to own directory
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
ns=[]
for i in aplddts:
    n='_'.join([i[1],'unrelaxed_model',i[2]])
    ns.append(n)
l=[i for i in os.listdir(nonaf2desdir) if i[-3:]=='pdb']
os.makedirs('af2filtered2',exist_ok=True)
for i in aplddts:
    n=str(i[1])+'.pdb'
    if n in l:
        np=os.path.join(nonaf2desdir,n)
        newp=os.path.join('af2filtered2',n)
        os.system('cp '+np+' '+newp)

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop/filtered_extrascores/cf/fastas/af2filtered ~/desktop/esled2fd3maf2filt

ONLY 17 VS 15 ORIGINALLY FILTERED DESIGNS HERE



ANALYSIS OF DIVERSITY
'''

import os
des=[i for i in os.listdir() if i[-3:]=='pdb']

motifs=[]
scaffs=[]
pairs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[2]
    scaffold=('_').join(s.split('_')[3:-11])
    motifs.append(motif)
    scaffs.append(scaffold)
    pairs.append((motif,scaffold))


# print(len(motifs))
# print(len(scaffs))
# print(len(pairs))
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(pairs)))
'''
502
261
523
'''
same_scaff_diff_match={}
for scaff in set(scaffs):
    l=[]
    for d in des:
        s=d.split('.')[0]
        scaffold=('_').join(s.split('_')[3:-11])
        if scaffold==scaff:
            match=s.split('_')[2]
            if match not in l:
                l.append(match)
    same_scaff_diff_match[scaff]=l

for key in same_scaff_diff_match.keys():
    if len(same_scaff_diff_match[key])>1:
        print(key)
        print(same_scaff_diff_match[key])
'''
'''
