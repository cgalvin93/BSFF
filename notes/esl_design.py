'''
esl
1 benzene
2 1st 2 6 member rings
3 whole carbon skel except methyl group
4 2 pentadiol hydroxyls, their hs, and the carbons they attached to
5 phenol


differences from last run :
GONNA ALLOW CHARGED RESIDUES IN MOTIFS

REMEMBER NOT TO ALLOW METHIONINE/PROLINE/CYSTEINE DURING DESIGN

GONNA BE USING NEW SCAFFOLD LIBRARY THAT I MADE ON ALL 1700 LUCS NTF2S FROM XINGJIE

USE BENZENE AS FRAG (last time used whole carbon skel)
    will also use whole carbon skel still, and both 6 membered rings
    GONNA USE DIHYDROXYL AS FRAG AND FIND EXPLICIT HBOND ROTAMERS THAT HBOND WITH
    BOTH OH, MAKE A SET OF MOTIFS FROM THESE AS WELL

WHEN REDESIGN FOR STABILITY W MPNN SEQ PROFILE, USE
    all seen in mpnn (as before)
    all seen more than like, 15% of the time
    only top 3-4 most seen
AND COMPARE RESULTS

CAN CONSIDER DERIVING 4PARAM HBOND CSTS FROM STATS FOR EACH THANG
'''




#create a useful directory structure for target
time python ~/desktop/BSFF/new_target.py esl

'''
gotta redo with smaller pentadiol frag
'''


'''
pubchem sdf
ONLY 1 CONFORMER FOR DOG
AVOGADRO
    hydrogens for ph 7.4
    optimize energy with mm94 ff




ex for conformers
AVOGADRO
    hydrogens for ph 7.4
    optimize energy with mm94 ff
    Extensions --> Molecular Mechanics --> Setup Force Field...
        mm94, default
    Extensions --> Molecular Mechanics --> Conformer Search...
        doesnt work?

                        conda activate pyr37
                        time python gencon.py dog.sdf 50 1000 0.1
'''


'''
create params
    in this case there is no carboxylate group neccessitating two differently protonated
    versions, just the one sdf file

compounds=("esl")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}.sdf
done

move params and pdb to Inputs/Rosetta_Inputs:
compounds=("esl")
for I in ${compounds[@]}; do
    mv ${I}.params Inputs/Rosetta_Inputs/${I}.params
    mv ${I}_0001.pdb Inputs/Rosetta_Inputs/${I}_0001.pdb
done
'''







'''
fragment search
with frags.txt files in target parent directories
DELETED ALKENE CHAIN FRAGMENT
MADE DIMETHYL SUBSTIT OFF RING MORE MINIMAL

esl
1 benzene
2 1st 2 6 member rings
3 whole carbon skel except methyl group
4 2 pentadiol hydroxyls, their hs, and the carbons they attached to
5 phenol

CARBONS HAVE TO BE IN ORDER LEAST TO GREATEST

compounds=("esl")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
    cd ..
done

'''



'''
NOW setup to split into multiple jobs for alignment on cluster
conda activate pyr37
compounds=("esl")
for I in ${compounds[@]}; do
    cd ${I}
    mkdir cluster_output
    mv *.pdb Inputs/Fragment_Inputs/
    mv *.mol Inputs/Fragment_Inputs/
    time python ~/desktop/BSFF/process_fragment_search.py
    cd ..
done




getting smiles...
924



okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r esl cgalvin@log2.wynton.ucsf.edu:










CHANGED EXTRACT_ALIGN TO USE RMSD 0.5 IF LESS THAN 10 ATOMS, 1.0 IF GREATER

in target par dir on cluster:
SUBMIT ALIGNMENT JOB
qsub -cwd -t 1-924 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh


'''


'''

SET UP SHELLFILE FOR SUBMITTING FILTER JOB
in target parent directory:
time python ~/BSFF/setup_filter_contacts.py
OR TO RUN ON CLUSTER
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/setup_filter_contacts.py
'>run_setupfilter.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_setupfilter.sh





OUTPUT OF THIS SCRIPT IS THE NUMBER YOU NEED FOR THE FOLLOWING ARRAY JOB

211



chmod ugo+x run_filter_array.sh
qsub -cwd -t 1-211 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh






CONSOLIDATE FILTERED CONTACTS
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/consolidate_filtered_contacts.py
'>run_consoldate.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_consoldate.sh














PROCESS FILTERED CONTACTS - do for each fragment
compounds=("Fragment_1" "Fragment_2" "Fragment_3" "Fragment_4" "Fragment_5")
for I in ${compounds[@]}; do
    time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/esl/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb esl $I
done

NOTHING FOR 4 AND 5

mkdir residue_contact_fuzzballs
mv *_residue_contacts residue_contact_fuzzballs






CLUSTER ALL
frags=["Fragment_1", "Fragment_2", "Fragment_3", "Fragment_4", "Fragment_5"]
import os
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','ILE','LEU','MET','PHE','VAL','GLY']

for fd in frags:
    fds=('_').join(['esl',fd,'residue_contacts'])
    os.chdir(fds)
    l=[]
    for i in polar_residues:
        for x in os.listdir():
            if x[-3:]=='pdb':
                if i in x:
                    l.append(os.path.join(os.getcwd(),x))

    f=open('sc.sh','w')
    f.write('#!/bin/bash')
    f.write('\n')
    f.write('source ~/anaconda3/etc/profile.d/conda.sh')
    f.write('\n')
    f.write('conda activate pyr37')
    f.write('\n')
    f.write('tasks=(0\n')
    for match in l[:-1]:
        f.write('       '+str(match)+'\n')
    f.write('       '+l[-1]+')')
    f.write('\n')
    cmd='time python3 ~/BSFF/spatial_cluster.py ${tasks[$SGE_TASK_ID]} /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb'
    f.write(cmd)
    f.write('\nqstat -j "$JOB_ID"')
    f.close()
    print(fds)
    print(len(l))
    os.chdir('..')

cd esl_Fragment_1_residue_contacts
mkdir clust_out
qsub -cwd -t 1-17 -l mem_free=16G -e clust_out -o clust_out sc.sh
cd ..
cd esl_Fragment_2_residue_contacts
mkdir clust_out
qsub -cwd -t 1-17 -l mem_free=16G -e clust_out -o clust_out sc.sh
cd ..
cd esl_Fragment_3_residue_contacts
mkdir clust_out
qsub -cwd -t 1-17 -l mem_free=16G -e clust_out -o clust_out sc.sh
cd ..
cd esl_Fragment_4_residue_contacts
mkdir clust_out
qsub -cwd -t 1-17 -l mem_free=16G -e clust_out -o clust_out sc.sh
cd ..
cd esl_Fragment_5_residue_contacts
mkdir clust_out
qsub -cwd -t 1-17 -l mem_free=16G -e clust_out -o clust_out sc.sh
cd ..






OKAY - MOTIF GENERATION
    firstly polar i will handle in 2 different ways
    the normal way first
'''
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
now i wanna find clusters that hbond with both pentadiol hydroxyls and
specify cst explicitly, along with adding another cst for the phenol oh


'''


import os
from pyrosetta import *
init('-pack_missing_sidechains False')
#
sf=get_fa_scorefxn()

#
ligparams=['/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params']
#
clusts_data={}
#
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
#
hbond_clusts=[]
for resdir in os.listdir():
    if os.path.isdir(resdir)==True:
        try:
            rn=resdir.split('_')[3]
            if rn in polar_residues:
                os.chdir(resdir)
                clusts=[i for i in os.listdir() if i[-3:]=='pdb']
                clust_data=[]
                for clust in clusts:
                    nhb=0
                    clust_pop=clust.split('_')[-1].split('.')[0]
                    nclust_pop=int(clust_pop)
                    p = Pose()
                    generate_nonstandard_residue_set(p,ligparams)
                    pose_from_file(p, clust)
                    res_scores=[]
                    #
                    ligand_resnum=p.total_residue()
                    for res_ind in range(1,ligand_resnum):
                        contact_pose=Pose()
                        ligand_pose=p.residue(ligand_resnum).clone()
                        res_pose=p.residue(res_ind).clone()
                        contact_pose.append_residue_by_jump(res_pose, 1)
                        contact_pose.append_residue_by_jump(ligand_pose, 1)
                        res_resname=contact_pose.residue(1).name()[0:3]
                        lig_resname=contact_pose.residue(2).name()[0:3]
                        if res_resname==lig_resname:
                            print('uh oh')
                            continue
                        else:
                            pass
                        hbond_set = rosetta.core.scoring.hbonds.HBondSet()
                        contact_pose.update_residue_neighbors()
                        rosetta.core.scoring.hbonds.fill_hbond_set(contact_pose, False, hbond_set)
                        if hbond_set.nhbonds()>1:
                            nhb+=1
                        if nhb>=(0.1*nclust_pop):
                            hbond_clusts.append(os.path.join(resdir,clust))
                            break
                os.chdir('..')
        except:
            pass


os.makedirs('double_hb')
rescounts={}
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
for i in polar_residues:
    c=0
    for j in hbond_clusts:
        rn=j.split('_')[3]
        if rn==i:
            cp=j.split('_')[-1].split('.')[0]
            cpn=int(cp)
            c+=cpn
    rescounts[i]=c
for j in hbond_clusts:
    os.system('cp '+j+' double_hb/'+j.split('/')[-1])

'''
{'SER': 196,
 'THR': 88,
 'GLN': 35,
 'ASN': 70,
 'TYR': 41,
 'TRP': 5,
 'GLY': 0,
 'HIS': 2,
 'ARG': 233,
 'LYS': 26,
 'GLU': 7,
 'ASP': 3}


 OKAY NOW IM GONNA GO IN THIS DOUBLE HB DIRECTORY AND TREAT IT PRETTY MUCH THE SAME WAY
 I WOULD NP CSTS
 '''
#first create random csts for just the phenolic oh,
#then to each i will randomly add double hb csts in the same way i would for np csts
import random
import os
#input ligand name and number of unique polar csts you want
n_desired_motifs=7
lig_resname='esl'
allowed_problematic_motifs=5
#new directory where motifs will be output
polar_cst_paths='dhb'
# polar_cst_paths='hb5'


#for each polar fragment, define whether acceptor is sp2 hybridized or not
#and what are the 3 constrained atoms
relevant_atoms={'Fragment_1':['n','O3','C18','C17']}

# #to put 2 hb on carbonyl
# relevant_atoms={'Fragment_1':['yes','O5','C23','C22'],
#                 'Fragment_2':['n','O1','C1','C3'],
#                 'Fragment_3':['n','O2','C6','C8'],
#                 'Fragment_4':['n','O3','C19','C18'],
#                 'Fragment_5':['yes','O5','C23','C22']}

#define wheteher or not vcharged residues are allowed to be used with frag
charged_allowed={'Fragment_1':'yes'}

#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
q_polar_residues=['LYS','ARG','SER','THR','TYR','TRP','GLY'] #no his
polar_residues=['SER','THR','TYR','TRP','GLY','HIS'] #no his
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
now add double hb csts

'''
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random
######
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep,fa_elec#, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
#####
# sf=get_fa_scorefxn()
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)
# min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
# mm4060 = pyrosetta.rosetta.core.kinematics.MoveMap()
# mm4060.set_bb(False)
# mm4060.set_chi(1,True)
# # mm4060.set_chi(1,True)
# # mm4060.set_jump(1,True)
# min_mover.movemap(mm4060)
# min_mover.score_function(sf)
# # min_mover.max_iter(1)
# min_mover.tolerance(1e-6)
##################################################
fragment_fuzzballs_dirs='/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_4_residue_contacts/double_hb'
ligparams=['/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params']
motifsoutdir='dhbmotifpdbs'
polar_cst_paths='/wynton/home/kortemme/cgalvin/esl/dhb'
############
single_contact_threshold=-1.0
##################################################
##################################################
##################################################
# we make sure we are in the right working directory
if not os.getcwd()==polar_cst_paths:
    os.chdir(polar_cst_paths)
# #make the nop motif pdb output directory
os.makedirs(motifsoutdir,exist_ok=True)


#get the paths of the residue fuzzball pdbs to open
clust_paths=[]
for clust in os.listdir(fragment_fuzzballs_dirs):
   if clust[-3:]=='pdb':
       rn=clust.split('_')[1]
       if rn!='SER':
           if rn!='THR':
               clust_paths.append(os.path.join(fragment_fuzzballs_dirs,clust))


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
for clust,freq in nps[:25]:         ####################################################################
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
In [10]: len(accepted_clusts)
Out[10]: 18
'''

#
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(15):          #this is how many unique hybrid cst i want to make per polar cst
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
    if len(lines)<20:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/dhb/dhbmotifpdbs ~/desktop/dhbmotifpdbs
OKAY INTERESTINGLY WHEN I LOOK AT THESE NONE OF THE SERINE ONES ACTUALLY HBOND WITH
BOTH HYDROXYLS SO IM GONNA REDO BUT WITHOUT THE SERINE CONTACTS







OKAY NOW I ACTUALLY WANT TO ADD NP CSTS TO ALL OF MY EXISTING POLAR CSTS
INCLUDING BOTH THE NORMAL WAY AND THE DOUBLE HB WAY

FOR THE NORMAL WAY, I WILL ADD 1NP
FOR THE DOUBLE HB WAY, I WILL ADD 2 NP CSTS
    (so will be 1polcst, 3 rotamer csts)


'''





'''
ADDING NP CSTS


1 benzene
2 1st 2 6 member rings
3 whole carbon skel except methyl group
4 2 pentadiol hydroxyls, their hs, and the carbons they attached to
5 phenol
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
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_1_residue_contacts','/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_2_residue_contacts','/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_3_residue_contacts','/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_5_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params']
motifsoutdir='npmotifpdbs'
polar_cst_paths='/wynton/home/kortemme/cgalvin/esl/dhb'
############
single_contact_threshold=-2.0
double_contact_score_threshold=-6.0
triple_contact_score_threshold=-9.0
##################################################
##################################################
##################################################
#we make sure we are in the right working directory
if not os.getcwd()==polar_cst_paths:
    os.chdir(polar_cst_paths)
#make the nop motif pdb output directory
os.makedirs(motifsoutdir,exist_ok=True)



dontallow=['ASN','GLN','ALA','GLY','PRO','CYS','ARG','LYS','SER','THR','GLU','ASP','MET']
#get the paths of the residue fuzzball pdbs to open
clust_paths=[]
clust_paths_pops=[]
for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
    tc=0
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    rn=clust.split('.')[0].split('_')[1]
                    if rn not in dontallow:
                        cp=clust.split('.')[0].split('_')[3]
                        ncp=int(cp)
                        tc+=ncp
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    rn=clust.split('.')[0].split('_')[1]
                    if rn not in dontallow:
                        cp=clust.split('.')[0].split('_')[3]
                        ncp=int(cp)
                        clustpath=os.path.join(fragment_fuzzballs_dir,i,clust)
                        clust_paths_pops.append((clustpath,float(ncp)/tc))


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
for clust,freq in nps[:50]:         ####################################################################
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
        for wp in range(int(freq*100)):
            np_cst_blocks.append(cstblock)
        clustnames_blocks.append((pdboutpath,cstblock))



# accepted_clusts
'''
In [13]: len(accepted_clusts)
Out[13]: 31
'''
#
# if os.path.exists('1np'):
#     os.system('rm -r 1np')
#     os.mkdir('1np')
# else:
#     os.mkdir('1np')
# os.chdir('1np')
# #to add 1 np cst to csts
# c=1
# for cst in os.listdir(polar_cst_paths):
#     if cst[-3:]=='cst':
#         bs=[]
#         for i in range(1):          #this is how many unique hybrid cst i want to make per polar cst
#             rblock=random.choice(np_cst_blocks)
#             if rblock not in bs:
#                 bs.append(rblock)
#                 f=open(os.path.join(polar_cst_paths,cst),'r')
#                 lines=[line for line in f.readlines()]
#                 f.close()
#                 ofname='hybrid_'+str(c)+'.cst'
#                 of=open(ofname,'w')
#                 of.write('CST::BEGIN\n')
#                 of.write('\n'.join(rblock))
#                 of.write('\n')
#                 of.write('CST::END\n')
#                 for line in lines:
#                     of.write(line)
#                 of.close()
#                 print(c)
#                 c+=1
#             else:
#                 rblock=random.choice(np_cst_blocks)
#                 if rblock not in bs:
#                     bs.append(rblock)
#                     f=open(os.path.join(polar_cst_paths,cst),'r')
#                     lines=[line for line in f.readlines()]
#                     f.close()
#                     ofname='hybrid_'+str(c)+'.cst'
#                     of=open(ofname,'w')
#                     of.write('CST::BEGIN\n')
#                     of.write('\n'.join(rblock))
#                     of.write('\n')
#                     of.write('CST::END\n')
#                     for line in lines:
#                         of.write(line)
#                     of.close()
#                     print(c)
#                     c+=1
#                 else:
#                     pass
#
#
# l=[i for i in os.listdir() if i[-3:]=='cst']
# print(len(l))
# for i in l:
#     f=open(i,'r')
#     lines=[line for line in f.readlines()]
#     f.close()
#     if len(lines)<40:
#         os.system('rm '+i)
# l=[i for i in os.listdir() if i[-3:]=='cst']
# print(len(l))
#
# os.chdir('..')





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
    if x[2]<=double_contact_score_threshold:
        cstblocks=[]
        for y in clustnames_blocks:
            if x[0]==y[0]:
                cstblocks.append(y[1])
            elif x[1]==y[0]:
                cstblocks.append(y[1])
        cstblockss.append((cstblocks[0],cstblocks[1]))

#cstblockss
'''
In [6]: len(cstblockss)
Out[6]: 54
'''

os.mkdir('2np')
os.chdir('2np')
#should do this in directory where you want them output
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(5): #how many unique per input polar cst
            rblocks=random.choice(cstblockss)
            if rblocks not in bs:
                bs.append(rblocks)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[0]))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblocks[1]))
                of.write('\n')
                of.write('CST::END\n')
                for line in lines:
                    of.write(line)
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
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[0]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblocks[1]))
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
    if len(lines)<30:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))



os.chdir('..')


'''
OKAY I THINK I FINALLY HAVE FOR THE DHB CSTS

'''





























'''
make pos files for new scaff lib
'''
import os
from pyrosetta import *
init()
l=[i for i in os.listdir() if i[-3:]=='pdb']
for i in l:
    p=pose_from_pdb(i)
    nres=p.total_residue()
    pfn=i.split('.')[0]+'.pos'
    of=open(pfn,'w')
    for x in range(2,p.total_residue()):
        of.write(str(x)+'  ')
    of.close()


'''
                            MATCHING
motifs:
/wynton/home/kortemme/cgalvin/esl/dhb/2np

scaffolds:
/wynton/home/kortemme/cgalvin/filtered_extrascores
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
shellfile_suffix='_esldhb.sh'
#
scaffold_path='/wynton/home/kortemme/cgalvin/filtered_extrascores'
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
    idx=[j for j in range(0,len(motif_paths[key]),120)]  #how many csts per job
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
528
n scaffolds
1816

'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('_1_a8s_3p1np.sh')
#
for j in jobss:
    cmd='qsub -cwd -t 1-1816 -l mem_free=16G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''
MATCHING ONLY TAKES ABOUT 1-3 HRS WITH 25 CSTS PER MATCH HERE
they also are taking less than 1g mem so i can definitely specify way less lmao

6-x hrs with occ now

15 hrs 50 ea job 1np3hb, maxvmem still no higher than 3 gb



TAKING LIKE A GOOD 24 HRS AT LEAST W/ 120 CSTS PER JOB USING MY NEW 1800 SCAFFOLD LIB
    probably wouldve been totally fine to just split intoi a couple more jobs to quicken things
    next time

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin





MAYBE DONT ALLOW POLAR RESIDUES IN NP CSTS
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
37556

37488

very few fucked up ones, thats good
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

f1=return_filtered(scores,'dsasa','>',0.7)
f2=return_filtered(scores,'ligsasasasa','<',50.0)
print(len(scores))
print(len(f1))
print(len(f2))
'''
36283
11848
6911

just take the dsasa filt
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



cd buried
mkdir genpot
mkdir cluster_output
'''
import os
matches=[i for i in os.listdir() if i[-3:]=='pdb']
count=1
sfname='ggp.sh'
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

'''
i shouldve done this earlier or baked it into the motif generation itself
but i wanna remove matches with more than 1 positively charged residue
or more than 1 gln/asn

'''
import os
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
for i in matches:
    nq=0
    nnq=0
    motif=i.split('_')[2]
    for j in motif:
        if j=='K' or j=='R':
            nq+=1
        elif j=='N' or j=='Q':
            nnq+=1
    if nq>=2:
        os.system('rm '+i)
    elif nnq>=2:
        os.system('rm '+i)
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
9391

4600
'''



#resfile generation
'''
NOTAA CPRKDE !!!!!
'''
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
4600
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
cd cluster_out


seems mostly ~10-15 designable, ~30-45 repackable
but there are still cases with tons of designable res (seeing one thats 27)


'''







# # changing cst files to match the pdb of matches names
##and move it to directory with pdbs,params,resfiles
import os
################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
cstdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np'
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
used wrong cstdir:
import os
l=[i for i in os.listdir() if i[-3:]=='cst']
for i in l:
    os.system('rm '+i)

'''









'''


ENZYME DESIGN APPLICATION



'''
# #ENZYME DESIGN WITH A RESFILE
import os
########
inputs_dir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean'
########
outputdir='enzdes'
os.makedirs(outputdir)
os.chdir(outputdir)
########
# allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_ligH_bcc.params'
sfname='esl_ed.sh'
ndesignseachmatch='75'
#########
matches=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-3:]=='pdb']
print(len(matches))
try:
    os.system('rm -r cluster_output')
except:
    pass
os.makedirs('cluster_output',exist_ok=True)
c=1
for xx in range(2):
    output_prefix='ed'+str(c)
    sf=open('_'+str(c)+sfname,'w')
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
         '-enzdes:cst_opt','-enzdes:bb_min',
         '-enzdes:chi_min',
         '-enzdes:cst_design','-enzdes:design_min_cycles 4', '-enzdes:lig_packer_weight 2.5',
         '-enzdes:cst_min',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights beta_genpot',
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
~/main/source/bin/enzyme_design.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_ligH_bcc.params -resfile /wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_clean.resfile -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt -enzdes:bb_min -enzdes:chi_min -enzdes:cst_design -enzdes:design_min_cycles 4 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_clean.cst -out:prefix ed1
'''
import os
sfname='esl_ed.sh'

jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4600 -l mem_free=2G -o cluster_output -e cluster_output '+j
    os.system(cmd)






'''
stopping early cus 690k designs maybe a lil od
plus im tryna new thing with disallowing charged res n shit so ill see how it
looks
LOOKING AT DESIGNS

damn theres a whoooole lotta jobs still running after 24 hr

/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

'''


import os
#################################################################################
# inputs_dir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes'
# allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
# prm1=allparams[0]
prm1='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/UM_9_W97L17Q66R13_1_relaxed_relaxed_366802_design_1_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_4_0001_0001_hybrid_297_1_0001_ligH_bcc.params'

sfname='qa.sh'
scorefile_name='enzdesanalysis.json'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
667798
'''
count=1
idx=[j for j in range(0,len(matches),200)]  #how many matches per job
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

man this shit taking waaayyyyy too long, need to crank it down to like 50
per job max or something
'''





'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes

'''




#analysis of scores from json file
sfname='enzdesanalysis.json'


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
'''
155704
338331
338322

77906
429220
429220
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

plot_dists(terms,scores,'esl_dhb.pdf')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    dontwork=[]
    for d in scores:
        try:
            vall=d[term]
            termval=float(vall)
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
            dontwork.append(d)
            # print('\n')
            # print(d)
            # print('\n')
    print(len(dontwork))
    return filtered_scores


f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'lighphobesasa','<',60.0)
f3=return_filtered(f2,'hbtolig','>',3.0)
f4=return_filtered(f3,'ligoversat','<',0.0)
# f5=return_filtered(f4,'bsE_per_res','<',-2.2)

# plot_dists(terms,f3,'esl4rmfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
# print(len(f5))
'''
429258
76210
28632
9920
9916
3519


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

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/esl_dhb.pdf ~/desktop/esl_dhb.pdf
'''
#make regular params
import os
import sys
ligname='esl'
cleandirname='analysis'
#####################
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(pdbs))
os.makedirs(cleandirname,exist_ok=True)
c=0
for initial_match in pdbs:
    print(c)
    c+=1
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

'''
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
        os.chdir('/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered') ############################
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
9878
'''
count=1
idx=[j for j in range(0,len(matches),5)]  #how many per job
if len(matches) not in idx:
    idx.append(len(matches))
for ei, ix in enumerate(idx[:-1]):
    sfn='_'+str(count)+sfname
    of=open(sfn,'w')
    of.write('#!/bin/bash\n')
    of.write('\n')
    for matchpath in matches[ix:idx[ei+1]]:
        of.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s='+matchpath+' -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+prm1+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
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
sfname='extrascore.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+j
    os.system(cmd)



'''
/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run
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
'''
9
9712
9712
4
9714
9714
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

plot_dists(terms,scores,'dhbfiltered1.pdf')
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
f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'boltz','<',-0.3)
f3=return_filtered(f2,'ddg','<',-24.0)
f4=return_filtered(f3,'contact_molsurf','>',150)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.65)
f7=return_filtered(f6,'lighphobesasa','<',50.0)
f8=return_filtered(f7,'ligoversat','<',0)
# f8=return_filtered(f7,'packstat','>',0.55)
# f9=return_filtered(f8,'buns_bb_heavy','<',5)
# f10=return_filtered(f9,'buns_sc_heavy','<',0)
# f12=return_filtered(f11,'oversat','<',0)
# f13=return_filtered(f12,'exphyd','<',1200)
# f14=return_filtered(f13,'cav','<',120)
#
# plot_dists(terms,f14,'4rm_testrun_filt.pdf')
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
# print(len(f9))
# print(len(f10))
# print(len(f11))
# print(len(f12))
# print(len(f13))
# print(len(f14))
'''
9714
7311
3843
1040
977
945
668
667
663
'''
import os
filtered_strc=[]
for d in f8:
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

/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2
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
663
56
182
213
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
1_relaxed_relaxed_180682_design_1_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_4_rank_2_0001_0001_hybrid
['F85I110Q23W39', 'W96L38Q53R92', 'F26L36R67S27']
1_relaxed_relaxed_165937_design_3_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_3_rank_2_0001_0001_hybrid
['F58F101N38W18', 'F101F58N34R62']
1_relaxed_relaxed_23588_design_2_unrelaxed_model_3_rank_1_0001_design_5_unrelaxed_model_1_rank_1_0001_0001_hybrid
['H93L45R100S43', 'F60L17N106R67', 'F60L17N106R63', 'F56F45R91T50']
1_relaxed_relaxed_362198_design_3_unrelaxed_model_3_rank_1_0001_design_5_unrelaxed_model_1_rank_2_0001_0001_hybrid
['H55F110K59S107', 'H55F110K59S31']
1_relaxed_relaxed_150776_design_3_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_3_rank_3_0001_0001_hybrid
['F17L98R44S14', 'F17L98R44S85']
1_relaxed_relaxed_289278_design_2_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_1_rank_5_0001_0001_hybrid
['F111W93Q44T84', 'F111W93Q44T86']
1_relaxed_relaxed_383537_design_1_unrelaxed_model_3_rank_1_0001_design_1_unrelaxed_model_2_rank_2_0001_0001_hybrid
['F87I112K44S96', 'F87I112K44S94']
1_relaxed_relaxed_181441_design_3_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_0001_hybrid
['F17L101K46S88', 'H59F115K70S113']
1_relaxed_relaxed_269053_design_1_unrelaxed_model_3_rank_1_0001_design_5_unrelaxed_model_3_rank_1_0001_0001_hybrid
['F71L103K40Y94', 'F71L103N64Y94']
1_relaxed_relaxed_20437_design_4_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_4_0001_0001_hybrid
['L17L43N93S116', 'H93L43R100S54']
1_relaxed_relaxed_212209_design_1_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_5_rank_2_0001_0001_hybrid
['F17L94N45K83', 'F17L94N45K81']
1_relaxed_relaxed_100376_design_3_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_2_rank_3_0001_0001_hybrid
['F41W56N85T59', 'H55F112K59S109']
1_relaxed_relaxed_264587_design_2_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_0001_hybrid
['F40F51R70S45', 'F51F40R82T58', 'F51F40R82T45']
1_relaxed_relaxed_317143_design_1_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_0001_hybrid
['F98I69K60S83', 'H60F110K67S34']
1_relaxed_relaxed_44883_design_1_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_0001_hybrid
['L17L46N96S44', 'L17L46N96S57']
1_relaxed_relaxed_168039_design_3_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_3_rank_3_0001_0001_hybrid
['F60L17Q43W61', 'F17L100N47K85']
1_relaxed_relaxed_23588_design_2_unrelaxed_model_3_rank_1_0001_design_1_unrelaxed_model_3_rank_3_0001_0001_hybrid
['L17L45N93S43', 'F60L17N106K63', 'F60L17N106R63']
1_relaxed_relaxed_193357_design_4_unrelaxed_model_5_rank_1_0001_design_2_unrelaxed_model_3_rank_4_0001_0001_hybrid
['L17F54R81T57', 'F54L17R81T57']
1_relaxed_relaxed_103933_design_3_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_0001_hybrid
['F95F55N32Y84', 'F95F55Q49Y84']
1_relaxed_relaxed_76383_design_4_unrelaxed_model_5_rank_1_0001_design_3_unrelaxed_model_1_rank_4_0001_0001_hybrid
['F17L95N55Y66', 'F17L95R51S82', 'F17L95R51S14']
1_relaxed_relaxed_23588_design_2_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_0001_hybrid
['L17F60N91T63', 'F60F45N104R63', 'F45F56R89S50', 'L17F60R89S63']
1_relaxed_relaxed_10200_design_2_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_3_rank_2_0001_0001_hybrid
['F109L17N45T97', 'F95F59N52Y63']
1_relaxed_relaxed_212209_design_1_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_3_rank_1_0001_0001_hybrid
['F17L94N45K14', 'F17L94N45K65', 'F17L94K48S92']
1_relaxed_relaxed_264587_design_2_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_5_rank_2_0001_0001_hybrid
['H54F109K58S99', 'H54F109K58S106']
1_relaxed_relaxed_23588_design_2_unrelaxed_model_3_rank_1_0001_design_1_unrelaxed_model_4_rank_1_0001_0001_hybrid
['F29F116Q46S36', 'F29F116Q46S32', 'F29F116Q63S32', 'L17L45N93S116']
1_relaxed_relaxed_189708_design_5_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_0001_hybrid
['F98F59N52R87', 'F98F59N52R96']
1_relaxed_relaxed_176872_design_1_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_4_rank_4_0001_0001_hybrid
['F58W17K28S87', 'F58W69K28S87']
1_relaxed_relaxed_280009_design_5_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_4_rank_3_0001_0001_hybrid
['F29F111Q70S108', 'F29F111Q59S108', 'L13F25N101T59', 'F29F111Q59S37', 'F29F111Q70S37', 'F29F111Q59S33']
1_relaxed_relaxed_300041_design_3_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_1_rank_3_0001_0001_hybrid
['F17L113K75S54', 'L17L40N88S54']
1_relaxed_relaxed_289278_design_2_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_1_rank_1_0001_0001_hybrid
['F17L97K46Y95', 'F52F41K14W46']
1_relaxed_relaxed_83902_design_1_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_3_rank_1_0001_0001_hybrid
['H109L114N33R44', 'H63F39K98S33']
1_relaxed_relaxed_163737_design_2_unrelaxed_model_5_rank_1_0001_design_4_unrelaxed_model_5_rank_5_0001_0001_hybrid
['F51F111N82R57', 'F51F111N82R43']
1_relaxed_relaxed_300911_design_4_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_4_rank_1_0001_0001_hybrid
['F17L102K49S100', 'F60L44K89S59']
1_relaxed_relaxed_34718_design_4_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_5_rank_2_0001_0001_hybrid
['F17L102N47K91', 'F17L102K47Y18']
1_relaxed_relaxed_290590_design_3_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_2_rank_3_0001_0001_hybrid
['F17L100K49S87', 'F17L100K49S71']
1_relaxed_relaxed_383811_design_3_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_5_rank_5_0001_0001_hybrid
['F58W67Q51T85', 'F58W67Q51T87', 'F17L96N45K81']
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
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='filtdesigns'
allfastasdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2'

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
plddt_threshold=85.0
carmsd_threshold=3.0
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
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2/cf/fastas/af2filtered ~/desktop/esl4rmtest_af2filtered


/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas/af2filtered


I WANNA EITHER SEE WHICH OF THESE HAVE ALSO GOOD ROSETTA STABILITY METRICS
ORRRR APPLY A MUCH STRICTER ALPHAFOLD FILTER IF I AM GOING TO DISREGARD THE ROSETTA STABILITY
METRICS
'''








#analysis of scores from json file
sfname='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/scores.json'


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

import os
affiltpdbs=[i for i in os.listdir() if i[-3:]=='pdb']

relscores=[]
for i in affiltpdbs:
    n=i.strip('.pdb')
    for j in scores:
        if j['decoy']==n:
            relscores.append(j)
len(relscores)
'''
9
9712
9712
4
9714
9714

In [4]: len(affiltpdbs)
Out[4]: 131
In [6]: len(relscores)
Out[6]: 131
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

plot_dists(terms,scores,'dhbaffilteredscores.pdf')
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
f1=return_filtered(relscores,'buns2interface','<',0.0)
f2=return_filtered(f1,'boltz','<',-0.3)
f3=return_filtered(f2,'ddg','<',-24.0)
f4=return_filtered(f3,'contact_molsurf','>',150)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.65)
f7=return_filtered(f6,'lighphobesasa','<',50.0)
f8=return_filtered(f7,'ligoversat','<',0)
f9=return_filtered(f8,'packstat','>',0.55)
f10=return_filtered(f9,'buns_bb_heavy','<',5)
f11=return_filtered(f10,'buns_sc_heavy','<',0)
f12=return_filtered(f11,'oversat','<',0)
f13=return_filtered(f12,'exphyd','<',1200)
f14=return_filtered(f13,'cav','<',120)
#
# plot_dists(terms,f14,'4rm_testrun_filt.pdf')
#
print(len(relscores))
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
9714
7311
3843
1040
977
945
668
667
663
'''
import os
filtered_strc=[]
for d in f14:
    # n=d['decoy']+'_0001.pdb'
    filtered_strc.append(d['decoy'])
os.makedirs('filtered3',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i+'.pdb filtered3/'+i+'.pdb')
    ########################################
os.system('mv *.pdf filtered3/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas/af2filtered/filtered3 ~/desktop/esldhbfilt3


DIVERSITY AMONGST FILTERED DESIGNS

cd filtered2

/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2
'''
# import os
# des=[i for i in os.listdir() if i[-3:]=='pdb']
# '''
# 'ed1UM_13_M99T39S17S95_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_8_unrelaxed_model_4_rank_1_0001_hybrid_39_1_0001_clean__DE_4_oglig_0001.pdb'
# '''
# motifs=[]
# scaffs=[]
# pairs=[]
# for d in des:
#     s=d.split('.')[0]
#     motif=s.split('_')[-9]
#     scaffold=('_').join(s.split('_')[3:-9])
#     motifs.append(motif)
#     scaffs.append(scaffold)
#     pairs.append((motif,scaffold))
#
#
# # print(len(motifs))
# # print(len(scaffs))
# # print(len(pairs))
# print(len(set(des)))
# print(len(set(motifs)))
# print(len(set(scaffs)))
# print(len(set(pairs)))
# '''
# scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/filtered2 ~/desktop/esl4rmfilt
# 663
# 56
# 182
# 213
# '''
# same_scaff_diff_match={}
# for scaff in set(scaffs):
#     l=[]
#     for d in des:
#         s=d.split('.')[0]
#         scaffold=('_').join(s.split('_')[3:-9])
#         if scaffold==scaff:
#             match=s.split('_')[2]
#             if match not in l:
#                 l.append(match)
#     same_scaff_diff_match[scaff]=l
#
# for key in same_scaff_diff_match.keys():
#     if len(same_scaff_diff_match[key])>1:
#         print(key)
#         print(same_scaff_diff_match[key])


































'''
RUNNING MPNN ON af2filt designs TO GENERATE SEQUENCE PROFILE

freezing residues that hbond w/ lig, or other residues that those res hb with
(bb-bb not included),
as well as res that have <= -2.0 REU pairwise interaxn w lig

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/dhbfiltered1.pdf ~/desktop/dhbfiltered1.pdf
/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2
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
paramsdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run'
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

# tfdata




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
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
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
# ofjsonname='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn_params.json'
# with open(ofjsonname,'r') as f:
#     tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
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
    '--num_seq_per_target 1000',
    '--sampling_temp "0.2"',
    '--batch_size 1']
    cmd=' '.join(cmdl)
    of.write(cmd)
    of.write('\n')
    of.write('qstat -j "$JOB_ID"')
    of.close()
    c+=1

'''
cd mpnn
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











consolidating mpnn results
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn'
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
ofjsonname='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn_params.json'
pdbsdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2'
paramsdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run'
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


#ALL SEEN RES IN RESFILE DESIGNABLE POSITIONS
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


#ONLY TOP 5 RES
os.makedirs('resfiles2',exist_ok=True)
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
        rescounts=[]
        for res1lc in list(set(positional_seq_list)):
            rc=0
            for res1lc2 in positional_seq_list:
                if res1lc2==res1lc:
                    rc+=1
            rescounts.append((res1lc,rc))
        rescounts=sorted(rescounts, reverse=True, key=lambda nmem: nmem[1])
        if len(rescounts)>=5:
            trcl=[a for a,b in rescounts[:5]]
            strc_seqs[des_ind]=trcl
        elif len(rescounts)>=4:
            trcl=[a for a,b in rescounts[:4]]
            strc_seqs[des_ind]=trcl
        elif len(rescounts)>=3:
            trcl=[a for a,b in rescounts[:3]]
            strc_seqs[des_ind]=trcl
        elif len(rescounts)>=2:
            trcl=[a for a,b in rescounts[:2]]
            strc_seqs[des_ind]=trcl
        else:
            trcl=[a for a,b in rescounts[:1]]
            strc_seqs[des_ind]=trcl
    seqs_data[strc]=strc_seqs
    #
    ofilename=os.path.join('resfiles2',strc+'.resfile')
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


'''
cd resfiles

cd resfiles2
'''
import os
#########
pdbsdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run/filtered2'
paramsdir='/wynton/home/kortemme/cgalvin/esl/dhb/2np/buried/genpot/clean/enzdes/filtered/analysis/run'
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
sfname='fd3bopmpnn2.sh'
ndesignseachmatch='20'
outputdirectory='fd_3bop_mpnn'
scorefile_name='esl_fd_3bop_mpnn2.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(5):     #this way i can spread across more nodes if neccessary
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
# sfname='fd3bopmpnnog.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-640 -l mem_free=4G -o cluster_out -e cluster_out '+j
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
    motif=s.split('_')[-11]
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
before mpnn
2
2
2

after mpnn
1
1
1
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
