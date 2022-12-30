got 110 ntf2 strcs submitted right now, lets say ill take top 10 for each bb again to cf,
then ill take top 5 for each, as long as they have plddt>80, so that at most i will
have 550 scaffolds
    (i like the idea of having multiple stable sequences for a given bb)
lets say in this case most i want is 500 motifs
and ima want more motifs for larger and hybrids
3p - 250 (wont match )
3p1np - 500
3p2np - 500
4p - 250
4p1np - 500
4p2np - 500



time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_1/Fragment_1_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_1
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_2/Fragment_2_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_2
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_3/Fragment_3_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_3
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_4/Fragment_4_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_4
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_5/Fragment_5_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_5
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_6/Fragment_6_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e Fragment_6

        DFHBI
    1-carbonyl
    2-secondary imidazole N
    3-alcohol
    4-benzene
    5-cf1
    6-cf2

mkdir 38efragfbs
mv *_residue_contacts 38efragfbs

scp cp38e.py cgalvin@log2.wynton.ucsf.edu:BSFF/
time python3 ~/BSFF/cp38e.py 38e /wynton/home/kortemme/cgalvin/round6/38e/38efragfbs /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb Fragment_1 Fragment_2 Fragment_3

process to copy fragment 1 results to fictional fragment 7 in a separate json that
i will submit for 4respolar

ENUMERATING POALR MOTIFS
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_contacts_clustered_liberal.json n_clust_each_frag ligand_resname nmotifres resultdir_prefix ndesired_motifs dist_cutoff_index
    dist_ind[5]=2.4 (allows d up to 2.5 ang)
    no more 1 arg
    no more 1 trp
    no more than 2 instances any other res
    residues are 'banned' from being in any more motifs once >1/4 the number of desired motifs
    hmmm, should i add a thing to not add too many charges?
            yea fuck it ima try that
            IT WOULD STILL INCREASE MY EFFICENCY HERE to merge clusters differing by 1 param n shit
        added no more than 1 arg/lys, changed cutoff for ind res to half of all motifs ,
        will up to 50 clust each frag


time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_contacts_clustered_liberal.json 40 38e 3 38e_3res_polar 500 5
3p - 250 (wont match )
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_contacts_clustered_liberal.json 40 38e 3 38e_3res_polar 250 5
3p1np - 500
3p2np - 500
4p - 250
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_contacts_clustered_liberal.json 40 38e 3 38e_3res_polar 250 5
4p1np - 500
4p2np - 500


echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/enumerate_polar_motifs_liberal.py 38e_polar_contacts_clustered_liberal.json 40 38e 3 38e_3res_polar 500 5
qstat -j "$JOB_ID"
'>epm3p_$compound_dir.sh
chmod ugo+x epm3p_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output epm3p_$compound_dir.sh



/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb
/wynton/home/kortemme/cgalvin/round6/38e/Transformed_Aligned_PDBs/Fragment_1/Fragment_1_contact_fuzzball.pdb
time python3 ~/BSFF/cluster_nonpolar.py 38e fragfbdir ligand_path frag1 frag2 frag3


time python3 ~/BSFF/cluster_nonpolar.py 38e /wynton/home/kortemme/cgalvin/round6/38e/38efragfbs /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb Fragment_4 Fragment_5 Fragment_6

compound_dir='38e'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/cluster_nonpolar.py 38e /wynton/home/kortemme/cgalvin/round6/38e/38efragfbs /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb Fragment_4 Fragment_5 Fragment_6
qstat -j "$JOB_ID"
'>cnp_$compound_dir.sh
chmod ugo+x cnp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cnp_$compound_dir.sh








#writes the csts to working directory
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random

fragment_fuzzballs_dir='/wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs/a8s_Fragment_4_residue_contacts'
ligparams=['/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params']
polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/a8s/a8s6_new_way_3respol_motifs_polar_csts'

clust_paths=[]
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


'''
hmmmm okay now how do i wanna define np csts for top clusts
    can score and use lowest scoring
    can simply take the first residue
i think the np contacts are sufficiently similar within each cluster that taking the first res
should be fine, so ill try that
'''
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
for clust,freq in nps[:20]:
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
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(sccs[0][0]).clone()
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


c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        f=open(os.path.join(polar_cst_paths,cst),'r')
        lines=[line for line in f.readlines()]
        f.close()
        ofname='hybrid_'+str(c)+'.cst'
        of=open(ofname,'w')
        for line in lines:
            of.write(line)
        rblock=random.choice(np_cst_blocks)
        of.write('CST::BEGIN\n')
        of.write('\n'.join(rblock))
        of.write('\n')
        of.write('CST::END\n')
        of.write('\n')
        of.close()
        print(c)
        c+=1












import os
#now submit for matching with lucs strc
lig_name='a8s'
allmotifsdir='/wynton/home/kortemme/cgalvin/round6/a8s/a8s6_new_way_3respol_motifs_polar_csts/hybridcsts'
shellfile_suffix='_a8s3p1npr6.sh'
paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
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
    for cst_path in motif_paths[key][1:]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        of.write('\n')
        of.write('tasks=(0\n')
        for scaff in scaffolds[:-1]:
            of.write('       '+os.path.join(scaffold_path,scaff.strip('.pdb'))+'\n')
        of.write('       '+os.path.join(scaffold_path,scaffolds[-1].strip('.pdb'))+')')
        of.write('\n')
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
                       '-match::lig_name', lig_name]
        cmd=' '.join(matcher_cmd)
        of.write(cmd+'\n')
        of.write('\nqstat -j "$JOB_ID"')
        print(count)
        count+=1
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)
print('n motifs')
print(len(paths))
print('n scaffolds')
print(str(len(scaffolds)))
'''
clean_1mpd.pdb
1
n motifs
404
n scaffolds
1505

qsub -cwd -t 1-1505 -l mem_free=32G -o cluster_out -e cluster_out _1_a8s3p1npr6.sh
