



'''

ANALYSIS OF esl BINDER IN THE PDB

had to clean the pdbs manually to remove extraneous homooligomeric chains
also got rid of the antibody binder cus it binds at heterodimeric interface
'''
# clean pdbs
import os
ligname='esl'
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
for i in pdbs:
    f=open(i,'r')
    alines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    f=open(i,'r')
    liglines=[line for line in f.readlines() if line[:6]=='HETATM' and line[17:20]==ligname or line[17:20]==ligname.upper()]
    f.close()
    print(liglines)
    newi=open(i.split('.')[0]+'_clean.pdb','w')
    for line in alines:
        newi.write(line)
    for line in liglines:
        newi.write(line)
    newi.close()
'''
mkdir clean
mv *clean.pdb clean
cd clean
'''
#make regular params
import os
import sys
ligname='esl'
cleandirname='analysis'
#####################
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
print(len(pdbs))
os.makedirs(cleandirname,exist_ok=True)
for initial_match in pdbs:
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
cd analysis
mkdir run
mv *oglig.pdb run
mv *.params run
cd run
'''
#simplify naming scheme
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
for i in pdbs:
    s=i.strip('_oglig.pdb')
    ss=s+'.pdb'
    os.system('mv '+i+' '+ss)
pdbs=[i for i in os.listdir() if i[-6:]=='params']
for i in pdbs:
    s=i.strip('pdb')
    ss=s
    os.system('mv '+i+' '+ss)



#run the extra scoring
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('extrascore.sh','w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('tasks=(0\n')
for match in pdbs[:-1]:
    sf.write('       '+match.strip('.pdb')+'\n')
sf.write('       '+pdbs[-1].strip('.pdb')+')')
sf.write('\n')
sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]}.pdb -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/nat_analysis.xml -in:file::extra_res_fa ${tasks[$SGE_TASK_ID]}.params -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
#
print(len(pdbs))
'''

qsub -cwd -t 1 -l mem_free=8G extrascore.sh

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

plot_dists(terms,scores,'esl_nat_extrascores.pdf')
'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/eslnat/clean/analysis/run/esl_nat_extrascores.pdf ~/desktop/esl_nat_extrascores.pdf
'''









































'''
i think i wanna try estradiol/estrione

check out pdb estriol binders


HOLY SHIT YOU KNOW WHAT, I THINK I FUCKED UP WITH THE POLAR CSTS BY USING THE LIGAND
HYDROGEN ATOM TO DEFINE THE INTERACTION, IT SHOULD BE O-C-C I THINK!
        yeah im gonna redo because of this problem ,
        but it is worth noting that i did it correctly for DOG and still
        ended up having problems with the hydrogen bonds in final designs...
'''







estetrol
estriol
estradiol
estrone (hydroxylated derivatives?)


#create a useful directory structure for target
time python ~/desktop/BSFF/new_target.py esl

'''
ONLY 1 CONFORMER FOR gonane
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



obabel -i mol esl.sdf -o pdb -O esl_ligH.pdb
time antechamber -i esl_ligH.pdb -fi pdb -o esl_bcc.sdf -fo sdf -c bcc -nc 0
~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n esl -p esl esl_bcc.sdf

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

GON
frag 1 - all carbons

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
4

getting smiles...
1


okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r esl cgalvin@log2.wynton.ucsf.edu:












SUBMIT ALIGNMENT JOB IN TARGET PARENT DIRECTORY
qsub -cwd -t 1 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh

                CHANGED CONTACT EXTRACTION RMSD THRESHOLD FROM 0.5 TO 1.0

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
dog
47
chmod ugo+x run_filter_array.sh
qsub -cwd -t 1 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh


#you dont do this when only 1 filter job or something
# CONSOLIDATE FILTERED CONTACTS
# echo '
# #!/bin/bash
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate pyr37
# time python ~/BSFF/consolidate_filtered_contacts.py
# '>run_consoldate.sh
# qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_consoldate.sh


PROCESS FILTERED CONTACTS - do for each fragment
compounds=("Fragment_1")
for I in ${compounds[@]}; do
    time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/esl/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb esl $I
done
mkdir residue_contact_fuzzballs
mv *_residue_contacts residue_contact_fuzzballs




CLUSTER NONPOLAR -
compound_dir='esl'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/cluster_nonpolar.py esl /wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb Fragment_1
qstat -j "$JOB_ID"
'>cnp_$compound_dir.sh
chmod ugo+x cnp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cnp_$compound_dir.sh



scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_1_residue_contacts ~/desktop/esl_clusters















GENERATE MOTIF CSTS -
so the general strategy here is to first create csts for the hbonds needed to
satisfy the ligand polar atoms,
for these we just use idealized hbond geometries
with protein residues as hydrogen donors
randomly selecting residue identities while only allowing at most 1 'problematic' residue
(ie one which introduces extraneous polar atoms into the binding site)

we will make X nonredundant polar csts like this

then we will identify the most statisticaly favored nonpolar contacts, and
evaluate their rosetta energy
and we will randomly add 1 or more of these contacts to our polar csts,
ensuring that the energies are favorable beneath some cutoff (and for
multiple np interactions, we eval the rosetta energy of both to make sure
they dont clash or anything)




'''




'''
IDEALIZED POLAR CSTS
'''


######
import random
import os
#input ligand name and number of unique polar csts you want
n_desired_motifs=100
lig_resname='esl'
allowed_problematic_motifs=20
#new directory where motifs will be output
polar_cst_paths='hb3_occ'
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
charged_allowed={'Fragment_1':'no',
                'Fragment_2':'no',
                'Fragment_3':'no'}

#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
q_polar_residues=['LYS','ARG','SER','THR','GLN','ASN','TYR','TRP','GLY'] #no his
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS'] #no his
problematic_residues=['LYS','ARG','ASP','GLU','GLN','ASN','HIS','TRP']
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



'''
SENSE THESE ARE HYBRID MOTIFS,
IM INCREASING CST SAMPLE NUMBERS AND TOLERANCES,
AND I WILL PUT THE NP CSTS FIRST, SO
THAT EVEN THOUGH THIS WILL EXPONENTIALLY INCREASE TIME AND MEM USAGE,
HOPEFULLY THE NP CSTS WILL ACT AS STRONG ENOUGH 'FILTERS'
THAT THE JOBS OVERALL ARENT USING A CRAZY AMOUNT OF MEMORY OR ANYTHING

'''

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
##################################################
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_1_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl.params']
motifsoutdir='npmotifpdbs'
polar_cst_paths='/wynton/home/kortemme/cgalvin/esl/hb3_occ'
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


sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep#, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
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
for clust,freq in nps[:10]:         ####################################################################
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
[('/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_4_residue_contacts/dog_Fragment_4_MET_contact_fuzzball_clustered/123_MET_cluster_68.pdb',
  68,
  -2.626019590179298),
 ('/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_1_residue_contacts/dog_Fragment_1_LEU_contact_fuzzball_clustered/793_LEU_cluster_55.pdb',
  55,
  -2.3369707866146467),
 ('/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_1_residue_contacts/dog_Fragment_1_PHE_contact_fuzzball_clustered/318_PHE_cluster_50.pdb',
  50,
  -1.7464081170760726),
 ('/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_1_residue_contacts/dog_Fragment_1_MET_contact_fuzzball_clustered/103_MET_cluster_33.pdb',
  33,
  -1.958057130224769)]
 '''

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

os.mkdir('2np')
os.chdir('2np')
#should do this in directory where you want them output
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(3): #how many unique per input polar cst
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
    if len(lines)<40:
        os.system('rm '+i)
l=[i for i in os.listdir() if i[-3:]=='cst']
print(len(l))



os.chdir('..')






#
#
#
#
# #i will go for 3 res motifs in this case!
# clust_cont_scores=[]
# for cind,clust in enumerate(pdboutpaths):
#     p = Pose()
#     generate_nonstandard_residue_set(p,ligparams)
#     pose_from_file(p, clust)
#     for clust2 in pdboutpaths[cind+1:]:
#         p2 = Pose()
#         generate_nonstandard_residue_set(p2,ligparams)
#         pose_from_file(p2, clust2)
#         for clust3 in pdboutpaths[cind+2:]:
#             p3 = Pose()
#             generate_nonstandard_residue_set(p3,ligparams)
#             pose_from_file(p3, clust3)
#             #
#             contact_pose=Pose()
#             res_pose=p.residue(1).clone()
#             res_pose2=p2.residue(1).clone()
#             res_pose3=p3.residue(1).clone()
#             contact_pose.append_residue_by_jump(res_pose, 1)
#             contact_pose.append_residue_by_jump(res_pose2, 1)
#             contact_pose.append_residue_by_jump(res_pose3, 1)
#             ligand_pose=p.residue(2).clone()
#             contact_pose.append_residue_by_jump(ligand_pose, 1)
#             score_this_one=sf(contact_pose)
#             clust_cont_scores.append((clust,clust2,score_this_one,clust3))
#
#
#
# sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[2])
# cstblockss=[]
# for x in sccs:
#     if x[2]<=triple_contact_score_threshold:
#         cstblocks=[]
#         for y in clustnames_blocks:
#             if x[0]==y[0]:
#                 cstblocks.append(y[1])
#             elif x[1]==y[0]:
#                 cstblocks.append(y[1])
#             elif x[3]==y[0]:
#                 cstblocks.append(y[1])
#         cstblockss.append((cstblocks[0],cstblocks[1],cstblocks[2]))
#
# #cstblockss
#
# os.mkdir('3np')
# os.chdir('3np')
# #should do this in directory where you want them output
# c=1
# for cst in os.listdir(polar_cst_paths):
#     if cst[-3:]=='cst':
#         bs=[]
#         for i in range(2): #how many unique per input polar cst
#             rblocks=random.choice(cstblockss)
#             if rblocks not in bs:
#                 bs.append(rblocks)
#                 f=open(os.path.join(polar_cst_paths,cst),'r')
#                 lines=[line for line in f.readlines()]
#                 f.close()
#                 ofname='hybrid_'+str(c)+'.cst'
#                 of=open(ofname,'w')
#                 of.write('CST::BEGIN\n')
#                 of.write('\n'.join(rblocks[0]))
#                 of.write('\n')
#                 of.write('CST::END\n')
#                 of.write('\n')
#                 of.write('CST::BEGIN\n')
#                 of.write('\n'.join(rblocks[1]))
#                 of.write('\n')
#                 of.write('CST::END\n')
#                 of.write('\n')
#                 of.write('CST::BEGIN\n')
#                 of.write('\n'.join(rblocks[2]))
#                 of.write('\n')
#                 of.write('CST::END\n')
#                 for line in lines:
#                     of.write(line)
#                 of.close()
#                 print(c)
#                 c+=1
#             else:
#                 rblocks=random.choice(cstblockss)
#                 if rblocks not in bs:
#                     bs.append(rblocks)
#                     f=open(os.path.join(polar_cst_paths,cst),'r')
#                     lines=[line for line in f.readlines()]
#                     f.close()
#                     ofname='hybrid_'+str(c)+'.cst'
#                     of=open(ofname,'w')
#                     of.write('CST::BEGIN\n')
#                     of.write('\n'.join(rblocks[0]))
#                     of.write('\n')
#                     of.write('CST::END\n')
#                     of.write('\n')
#                     of.write('CST::BEGIN\n')
#                     of.write('\n'.join(rblocks[1]))
#                     of.write('\n')
#                     of.write('CST::END\n')
#                     of.write('\n')
#                     of.write('CST::BEGIN\n')
#                     of.write('\n'.join(rblocks[2]))
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
# l=[i for i in os.listdir() if i[-3:]=='cst']
# print(len(l))
# for i in l:
#     f=open(i,'r')
#     lines=[line for line in f.readlines()]
#     f.close()
#     if len(lines)<50:
#         os.system('rm '+i)
# l=[i for i in os.listdir() if i[-3:]=='cst']
# print(len(l))
#
#
#
# os.chdir('..')

'''

OKAY GONNA LET MY DOG BINDER DESIGN FINISH BEFORE SUBMITTING MATCHING ON THIS
    remember to use enough memory when i do cus the cst sample numbers are all
    very high

    i will at least attempt to match these 3+3 hybrids i guess, but probably not the best
    chances of success
        ACTUALLT HOW ABOUT THIS, I WILL FIRST TRY MATCHING THE 2+3 HYBRIDS, AND DEPENDING
        ON HOW MANY I GET I CAN DECIDE IF I WANNA TRY THE 3+3 OR MAYBE 3+2 OR something
        (or otherwise just proceed with designing the 2+3 matches)
            AHHH SHIT BUT IN THIS CASE PERHAPS I WANNA MAKE MORE 2+3s
            (theres only like 150, i really probably should)

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  651.02 GiB| 1000.00 GiB||  3354874|unlimited

MAYBE DELETE THE OLD DOG DESIGNS BEFORE MATCHING
rm -r fd
rm -r fd_3bop
rm -r dog_cm
rm -r dog_ed



DOG ED2 JOBS TOOK ABOUT 10 HOURS
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
'''

































'''genpot'''

#enter this matches directory
#now get params for each and put in new directory with corrected pdb
import os
standardmatches=[i for i in os.listdir() if i[:2]=='UM']
print(len(standardmatches))
'''
2np3hb
29
IN THIS CASE THIS IS SO FEW MATCHES THAT ILL JUST TRY DESIGNING ON ALL OF THEM

1np3hb
37136

occ 2np3hb
208

'''
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
29

29


37136
36866

208


208
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
4951
3099
'''
import os
filtered_strc=[]
for d in f1:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i[:-5]+'.pdb filtered/'+i+'.pdb')




























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
idx=[j for j in range(0,len(matches),1)]  #how many matches per job
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


LOST 6 MATCHES TO GENPOT PARAMS GENERATION (29 TO 23)
'''


import os
#have to fix the params files to make sure no redundant bonds defined :(
prms=[i for i in os.listdir() if i[-6:]=='params']
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
count=1
idx=[j for j in range(0,len(matches),1)]  #how many matches per job
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
        cmdd='time python ~/BSFF/tools/match_resfile.py '+matchpath+' '+prm1+' 5.0 genpot'
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








# # changing cst files to match the pdb of matches names
##and move it to directory with pdbs,params,resfiles
import os
################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np'
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np'
cstdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np'
################################################
for i in matches:
    ################################################
    cstid='hybrid_'+i.split('_')[-3]+'.cst'
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
# # coupledmoves
# # setup so that 1 params for all
# # resfile and cst must have same name as pdb
# import os
# ########
# allparams=[i for i in os.listdir() if i[-6:]=='params']
# prm1=allparams[0]
# sfname='dog_cm.sh'
# ndesignseachmatch='25'
# outputdirectory='dog_cm'
# scorefile_name='dog_cm.json'
# #########
# matches=[i for i in os.listdir() if i[-3:]=='pdb']
# os.makedirs('cluster_output',exist_ok=True)
# os.makedirs(outputdirectory,exist_ok=True)
# c=1
# for xx in range(2):
#     output_prefix='cm'+str(c)
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
#          '-parser:protocol','~/BSFF/tools/cm_genpot_jump1.xml',
#          '-parser:view','-run:preserve_header',
#          '-in:file::extra_res_fa',prm1,
#          '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
#          '-mute protocols.backrub.BackrubMover',
#          '-ex1','-ex2','-extrachi_cutoff','0',
#          '-out:nstruct',ndesignseachmatch,
#          '-score::weights beta_genpot',
#          '-corrections::gen_potential',
#          # '-score::weights','ligand',
#          # '-restore_pre_talaris_2013_behavior',
#          # '-score:set_weights','hbond_bb_sc','2.0',
#          # '-score:set_weights','hbond_sc','2.0',
#          '-load_PDB_components','False',
#          '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
#          '-out:path:all',outputdirectory,
#          '-out:prefix',output_prefix,
#          '-scorefile_format','json',
#          '-out:file:scorefile',scorefile_name,
#          '-coupled_moves::mc_kt 2.4',
#          '-coupled_moves::boltzmann_kt 2.4', '-coupled_moves::ntrials 1000',
#          '-coupled_moves::initial_repack false', '-coupled_moves::ligand_mode true',
#          # '-coupled_moves::ligand_weight 2',
#          '-coupled_moves::fix_backbone false',
#          '-coupled_moves::bias_sampling true', '-coupled_moves::bump_check true',
#          '-coupled_moves::exclude_nonclashing_positions true',
#          '-coupled_moves::backbone_mover kic',
#          '-coupled_moves::kic_perturber walking']
#     sf.write((' ').join(cmd))
#     sf.write('\nqstat -j "$JOB_ID"')
#     sf.close()
#     c+=1
# '''
# print(len(matches))
# 4331
#
# ~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s  UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.pdb -parser:protocol ~/BSFF/tools/cm_genpot_jump1.xml -parser:view -run:preserve_header -in:file::extra_res_fa UM_10_N107S110T17Y58_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_cst_4_1_0001_ligH_bcc.params -resfile  UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.resfile -mute protocols.backrub.BackrubMover -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 1 -score::weights beta_genpot -corrections::gen_potential -load_PDB_components False -enzdes:cstfile  UM_9_Y90T57N99T50_1_relaxed_relaxed_5tpj_225840_design_7_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_cst_104_1_0001_clean.cst -out:path:all dog_cm -out:prefix cm2 -scorefile_format json -out:file:scorefile dog_cm.json -coupled_moves::mc_kt 2.4 -coupled_moves::boltzmann_kt 2.4 -coupled_moves::ntrials 1000 -coupled_moves::initial_repack false -coupled_moves::ligand_mode true -coupled_moves::fix_backbone false -coupled_moves::bias_sampling true -coupled_moves::bump_check true -coupled_moves::exclude_nonclashing_positions true -coupled_moves::backbone_mover kic -coupled_moves::kic_perturber walking
#
# qsub -cwd -t 1-3657 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh
# '''
# import os
# sfname='dog_cm.sh'
# jobss=[i for i in os.listdir() if sfname in i]
# for j in jobss:
#     cmd='qsub -cwd -t 1-3657 -l mem_free=4G -o cluster_output -e cluster_output '+j
#     os.system(cmd)
#
#
#
#
# '''
# rm *last.pdb
# rm *low.pdb
# rm *.fasta
# rm *.stats
# -bash: /usr/bin/rm: Argument list too long
#
#
#
# import os
# l=[i for i in os.listdir() if i[-8:]=='last.pdb']
# l2=[i for i in os.listdir() if i[-7:]=='low.pdb']
# l3=[i for i in os.listdir() if i[-6:]=='.fasta']
# c=1
# print(len(l))
# for i in l:
#     os.system('rm '+i)
#     c+=1
#     print(str(c)+'/'+str(len(l)))
# c=1
# for i in l2:
#     os.system('rm '+i)
#     c+=1
#     print(str(c)+'/'+str(len(l2)))
# c=1
# for i in l3:
#     os.system('rm '+i)
#     c+=1
#     print(str(c)+'/'+str(len(l3)))
#
#
# echo "
# import os
# l=[i for i in os.listdir() if i[-8:]=='last.pdb']
# l2=[i for i in os.listdir() if i[-7:]=='low.pdb']
# l3=[i for i in os.listdir() if i[-6:]=='.fasta']
# c=1
# for i in l:
#     os.system('rm '+i)
# c=1
# for i in l2:
#     os.system('rm '+i)
# c=1
# for i in l3:
#     os.system('rm '+i)
# ">cleancm.py
#
# echo '
# #!/bin/bash
# time python3 cleancm.py
# qstat -j "$JOB_ID"
# '>cleancm.sh
# qsub -cwd -l mem_free=4G cleancm.sh

'''


ENZYME DESIGN APPLICATION



'''
import os
########
inputs_dir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean'
########
outputdir='enzdes'
os.makedirs(outputdir)
os.chdir(outputdir)
########
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
sfname='esl_ed2.sh'
ndesignseachmatch='10'
#########
matches=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-3:]=='pdb']
try:
    os.system('rm -r cluster_output')
except:
    pass
os.makedirs('cluster_output',exist_ok=True)
c=1
for xx in range(20):
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


'''
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-163 -l mem_free=2G -o cluster_output -e cluster_output '+j
    os.system(cmd)
'''

WHY ONLY 163 INSTEAD OF 208 MATCHES, HOW AM I LOSING SO MAN Y DURING GP PARAMS GEN?



'''
import os
#################################################################################
allparams=[i for i in os.listdir('/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean') if i[-6:]=='params']
prm1=os.path.join('/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean',allparams[0])
sfname='qa.sh'
scorefile_name='esled2.json'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
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
             '-s',matchpath,
             '-parser:protocol','~/BSFF/tools/edanalysis.xml',
             '-parser:view','-run:preserve_header',
             '-in:file::extra_res_fa',prm1,
             '-score::weights beta_genpot',
             '-corrections::gen_potential',
             '-load_PDB_components','False',
             '-scorefile_format','json',
             '-out:file:scorefile',scorefile_name]
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
only taking about 2 hrs with the 350 per job

'''





























'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE
SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE

/wynton/home/kortemme/cgalvin/esl/hb3/2np/genpot/clean/enzdes
'''




#analysis of scores from json file
# sfname='dog_fd.json'
# sfname='dog_fd_3bop.json'
# sfname='dog_cm.json'
# sfname='doged.json'
# sfname='doged2.json'
# sfname='esled.json'
sfname='esled2.json'


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

plot_dists(terms,scores,'esl_ed2.pdf')

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
plot_dists(terms,f3,'doged2filted2.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
'''
2232
280
148
59


LOSING SO MANY SCORES TO JSON FORMATTING ISSUES HERE, REALLY NEED TO ADDRESS THIS
31228
6024
2443
1940
'''
import os
filtered_strc=[]
for d in f3:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
os.system('mv *.pdf filtered/')

'''
cd filtered



scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/hb4/filtered/genpot/clean/dog_ed/filtered ~/desktop/doged
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
for initial_match in pdbs:
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



'''
#run the extra scoring
import os
params=[i for i in os.listdir() if i[-6:]=='params']
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
qsub -cwd -t 1-1805 -l mem_free=2G extrascore.sh
        took about an hour per job

down from 1940 on og params generation

/wynton/home/kortemme/cgalvin/esl/hb3/2np/genpot/clean/enzdes/filtered/analysis/run
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

# plot_dists(terms,scores,'esled2_extrascores.pdf')
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
f3=return_filtered(f2,'ddg','<',-22.0)
f4=return_filtered(f3,'contact_molsurf','>',160)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.66)
f7=return_filtered(f6,'lighphobesasa','<',20.0)
plot_dists(terms,f5,'esled2_extrascoresfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
print(len(f6))
print(len(f7))
'''
59
59
57
39
24
18

1740
1658
565
329
263
257
155
155
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
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/esled2_extrascoresfilt.pdf ~/desktop/esled2_extrascoresfilt.pdf



            ed1 1np3hb
I DONT UNDERSTAND I FILTERED THESE SPECIFICALLY ON HAVING 3 FUCKING HBONDS
AND YET I DONT SEE ANY OF THEM WITH 3 FUCKING HBONDS WHAT THE ACTUAL FUCK
IS THE PROBLEM HERE

IT ALSO LOOKS JUST GENERALLY LIKE THE HBONDS ARE COMPLETELY FUCKED UP
THE LIGAND POTENTIALLY MOVED FAR AWAY FROM THEM
    ON 1 THERES A T THATS SUPPOSED TO HB WITH OH BUT INSTEAD THERES A FUCKING SER
    RIGHT NEXT TO IT DOING THE HB WHILE THE THR HANGS OUT ALONE WTF IS GOING ON HERE?



            ed2 2np3hb
pstat, bb buns could be addtnl filters, but since im planning to
design for further stability anyway prolly mkes sense to look at those after

I THINK I WANNA FREEZE WITH MPNN ANYTHING HBONDING WITH LIGAND (or hbonding with those res)
OR MAKING STRONG INTERACTION WITH IT, NOT NECCESSARILY CST RES EVEN
    i am seeing a case where cst asn is turned out to surface of ptn and a serine
    was designed in to satisfy the lig polar atom left unsat by loss of cst asn
THEN USE 3BOP WHEN I DO DESIGN

LEMME DO SOME ANALYSIS HERE OF WHAT KINDA MATCH/MOTIF/SCAFOLD DIVERSITY I HAVE HERE
'''
import os
des=[i for i in os.listdir() if i[-3:]=='pdb']

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
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(pairs)))
'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/UM_1_L53F106Y57S13S35_1_relaxed_relaxed_5tpj_67405_design_4_unrelaxed_model_5_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_hybrid_104_1.pdb ~/desktop/matchex.pdb
12
18
21
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
1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid
['L51F113S39N30S108', 'L51F113S39N30S35']
1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_4_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_hybrid
['L51M99S39T108N31', 'L51M99Y59T108S27']
1_relaxed_relaxed_5tpj_128842_design_5_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_4_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid
['L109M17Y84T25S53', 'L109M17Y84T54S53']
1_relaxed_relaxed_5tpj_67405_design_4_unrelaxed_model_5_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid
['L53F106T56Q101T35', 'L53F106S56N30S101']
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
allfastasdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2'

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
import numpy as np
#id designs with plddt over threshold
#use best vals:
plddt_threshold=88.0
carmsd_threshold=1.2
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
os.makedirs('af2filtered',exist_ok=True)
for i in aplddts:
    n=str(i[1])+'.pdb'
    if n in l:
        np=os.path.join(nonaf2desdir,n)
        newp=os.path.join('af2filtered',n)
        os.system('cp '+np+' '+newp)

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/cf/fastas/af2filtered2 ~/desktop/esled2af2filt2

OKAY THESE 10 DESIGNS ACTUALLY LOOK PRETTY DARN GOOD, MAYBE WORTH ORDERING

however, lemme see if i can improve overall af2 confidence/rosetta stability metrics
using mpnn informed design and 3bop

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
filters_xml = f'''
                <FILTERS>
                    <DSasa name="dsasa" confidence="0.0"/>
                </FILTERS>'''
dsasa_filter = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(filters_xml).get_filter("dsasa")
#
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
tfdata={}

# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn_params.json'
#########
# with open(ofjsonname,'r') as f:
#     ttfdata=json.load(f)
# #
# # for pdb in ttfdata.keys():
for pdb in pdbs:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    dsasa_val=dsasa_filter.report_sm(p)
    if dsasa_val>0.8:
        tofreeze=[]
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










'''
THESE FINAL DESIGNS LOOK PRETTY GOOD, BUT THE THING IS I WANT GREATER DIVERSITY,
I NEED A SET OF DIFFERENT DESIGNS THAT LOOK THIS GOOD
    now i can do something to maybe take the top best from each match or something above
    some relaxed cutoff, and then put extra dedicated effort to improve the metrics
    that some designs are lacking in
    i can also simply make more designs,
    or try to get more matches, eg by matching only 1np
    i can also try using smaller fragments for bsff again, might ultimately discover some
        higher affinity interaxns this way
    what is the alphafold confidence for ntf2 in the pdb?????maybe its low generally
        because the fold is so strange
            similarly, what are the rosetta metrics, i for instance know that they
            tend to have a few buried unsats from xingjies design of them
    i can also try 3bop during initial design stage again seeing as how i realized that i
    botched it the first time and did not specify the correct scorefxn for the fd application
MAKE PRESENTATION WHILE I THINK ABOUT NEXT STEPS,
    or maybe if anything just start like a 1np3hb matching or something like thgat which will
    take a while and can run in the bg while i prepare that

MAYBE 1NP but i just use like a select 3 w/ best energy/ that i like the most

YOOOOOO I NEED TO LOOK AT CARMSD OF JUST BINDING SITE RESIDUES!

'''



'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS
making a try at doing just bs residues

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-load_PDB_components False')

directory_prefix='filtdesigns'
allfastasdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop/filtered_extrascores/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop/filtered_extrascores'
paramsdir='/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles'
allparams=[os.path.join(paramsdir,i) for i in os.listdir(paramsdir) if i[-6:]=='params']
prm1=allparams[0]
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
for p in prediction_pdb_paths[:5]:
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
        lig=[prm1]
        p22=Pose()
        generate_nonstandard_residue_set(p22,lig)
        pose_from_file(p22, os.path.join(nonaf2desdir,csnogds))
        all_carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p22,p11)
        #########
        p22.update_residue_neighbors()
        p11.update_residue_neighbors()
        #
        ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
        neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector(p22.total_residue())
        # neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 8., False)
        neighborhood_selector_bool = neighborhood_selector.apply(p22)
        neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
        first_shell_res=sorted(list(neighborhood_residues_resnums))
        ##############
        og_bs=Pose()
        for resnum in first_shell_res:
            res_pose=p22.residue(resnum).clone()
            og_bs.append_residue_by_jump(res_pose, 1)
        #####################
        ##############
        af_bs=Pose()
        for resnum in first_shell_res:
            res_pose=p11.residue(resnum).clone()
            af_bs.append_residue_by_jump(res_pose, 1)
        #####################
        carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(og_bs,af_bs)
    ##########
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    #########
    currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
    ##################
    print('\n\n\n')
    print(str(len(first_shell_res)))
    print(str(all_carmsd))
    print(carmsd)
    print('\n\n\n')
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#
#
# import os
# import numpy as np
# import json
# from pyrosetta import *
# init('-load_PDB_components False')
# p11=pose_from_pdb('/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop/filtered_extrascores/cf/fastas/filtdesigns_2/filtdesigns_2_output/fd13ed6UM_1_L54F94Y40W98Y31_1_relaxed_relaxed_5tpj_322479_design_3_unrelaxed_model_3_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_hybrid_5_1_clean__DE_3_0001_oglig_0001_0001_0001_unrelaxed_model_4_rank_5.pdb')
# lig=['/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/ed13UM_1_L51M99Y18T108S48_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_193_1_clean__DE_2_0001_oglig_0001.params']
# p22=Pose()
# generate_nonstandard_residue_set(p22,lig)
# pose_from_file(p22, '/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/genpot/clean/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_3bop_mpnn/filtered_fd3bop/filtered_extrascores/fd13ed6UM_1_L54F94Y40W98Y31_1_relaxed_relaxed_5tpj_322479_design_3_unrelaxed_model_3_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_hybrid_5_1_clean__DE_3_0001_oglig_0001_0001_0001.pdb')
# all_carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p22,p11)
# #########
# p22.update_residue_neighbors()
# p11.update_residue_neighbors()
# #
# ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
# neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 8., False)
# neighborhood_selector_bool = neighborhood_selector.apply(p22)
# neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
# first_shell_res=sorted(list(neighborhood_residues_resnums))
# #
# og_bs=Pose()
# for resnum in first_shell_res:
#     res_pose=p22.residue(resnum).clone()
#     og_bs.append_residue_by_jump(res_pose, 1)



#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
##################################################################
##################################################################
##################################################################

# import numpy as np
#id designs with plddt over threshold
#use best vals:
plddt_threshold=88.0
carmsd_threshold=1.2
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
os.makedirs('af2filtered',exist_ok=True)
for i in aplddts:
    n=str(i[1])+'.pdb'
    if n in l:
        np=os.path.join(nonaf2desdir,n)
        newp=os.path.join('af2filtered',n)
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
ogscaffs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-11]
    scaffold=('_').join(s.split('_')[3:-12])
    ofs=s.split('_')[7]
    motifs.append(motif)
    scaffs.append(scaffold)
    pairs.append((motif,scaffold))
    ogscaffs.append(ofs)

motifs=[]
scaffs=[]
pairs=[]
ogscaffs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-2]
    scaffold=('_').join(s.split('_')[3:-3])
    ofs=s.split('_')[7]
    motifs.append(motif)
    scaffs.append(scaffold)
    pairs.append((motif,scaffold))
    ogscaffs.append(ofs)

# print(len(motifs))
# print(len(scaffs))
# print(len(pairs))
print(len(set(pairs)))
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(ogscaffs)))
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

same_scaff_diff_match={}
for scaff in set(scaffs):
    l=[]
    for d in des:
        s=d.split('.')[0]
        scaffold=('_').join(s.split('_')[3:-3])
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

CHIMERA MOVIE MAKING


movie record ; turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait; movie encode ~/desktop/esl_binder_spin.mp4
'''
























'''
Diffdock
Ptn lig strc prediction https://arxiv.org/pdf/2209.15171.pdf
How to do yeast display? ID a set of strc in current designs to express, ITS TIME
Degrado sm binder meetings https://ucsf.zoom.us/j/95065542919?pwd=TUtaNkhwbXBBWHdYbHhSYU5LN0lhZz09
	every Tuesday 4-5
MAYBE GO BACK TO XJ FILT NTF2, GET ALL OF PROP LENGTH, MPNN REDESIGN, THEN MATCH (THAT SMALL ENUF 4 YEAST DISPLAY)




Same n res, only 86 variable,
Need contiguous terminal set
Of unchanging res






1709 lucs ntf2 scaffs from xj


scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/hb3_occ/2np/matchex ~/desktop/eslmatches



I'M THINMING THE MOVE IS TO EXPAND MY NTF2 LIBRARY




okay while i am working on all the ntf2 strc,
it seems that colabfold predictions for them are gonna take a long ass time,
so in the meantime, maybe it is a good move to think more about what constitutes
idealized hydrogen bonding or nonpolar interactions

would be interesting to redesign some existing binders with solved structures,
espeically in ntf2 folds
    pfam id PF02136
    2Z77 is actually a ntf2 with a estradiol based
        steroid in it, putative steroid isomerase
    EQU equilenin intertesting steroid in 24 pdb entries, including another ntf2 3m8c

MYABE A COMMON KEONT BINDING MOTIF IS TWO RESIDUES WHICH CAN HBOND WITH
EACH OTHER IN APO STATE

OKAY IDENTIFY SET OF NTF2S THAT BIND THINGS
        -af2/rosetta metrics
    and then subset of ntf2s that bind steroids
        -^^^
        -conservation of np/hb interaxns
        -matching my own motifs into these scaffolds,
            or redesigning binding sites and looking at pps/binding metric comparison
            or matching the true motif into the true scaffold and then designing on it?
17 STEROID BINDING NTF2S W RES UNDER 2.5 ANGSTROMS


-I WANT GREATER DIVERSITY IN MY FINAL SET OF DESIGNS,
    try matching 3hb`np and see if i get comparable packing/ddg values n shit but hopefully more diversity
-i want to fix the hydroxyl oritentations (maybe try no 3bop on mpnn informed design? or just
minimizing without 3bop in scorefunction? i think i do before extrascores
but id need to double check:YES THERE IS REF15 MINIMIZATION, maybe need to repack
and not just minimize?)

GOTTA DO A 1-TO-1 COMPARSION WITH TRYING 3HB1NP MATCHING ANYWAY TO
LEARN IF IT PRODUCES COMPARABLE DDGS, SO I MIGHT AS WELL JUST GET THAT strated using
current ntf2 liubrary while i await new library to finish developing

IN TERMS OF HBONDING,
idea of having two targeting lig polar group that can hbond w each other in
apo state is cool ,
also this idea of low barrier hbonds where heteroatoms have similar pka


EVIDENTLY THERE IS A GENERALIZABLE WAY TO DEFINE AROMATIC PI STACKING IN
MATCH CST FILE BY ADDING PSEUDOATOM X AT CENTER OF RING THAT YOU DEFINE AS BONDED
TO ANY OTHER ATOM IN THERE,THEN DO YOUR WHOLE MOLFILE TO PARAMS WITH THAT,
THEN DEFINE CST RELATIVE TO PSEUDOATOM


https://www.rosettacommons.org/content/how-define-pi-pi-stacking-cst-file
 enzdes integration test, which is describing a cysteine esterase active site a
CST::BEGIN
TEMPLATE:: ATOM_MAP: 1 atom_name: X1 C10 C12
TEMPLATE:: ATOM_MAP: 1 residue3: D2N

TEMPLATE:: ATOM_MAP: 2 atom_type: aroC,
TEMPLATE:: ATOM_MAP: 2 residue1: WFY

CONSTRAINT:: distanceAB: 3.50 0.20 50.00 0
CONSTRAINT:: angle_A: 90.00 5.00 50.00 360.00
CONSTRAINT:: angle_B: 90.00 5.00 50.00 360.00
CONSTRAINT:: torsion_A: 90.00 5.00 50.00 180.00
CONSTRAINT:: torsion_B: 180.00 15.00 0.00 120.00
CONSTRAINT:: torsion_AB: 0.00 0.00 0.00 180.00
The easiest way to get the coordinates is to average the six coordinates of the carbon atoms in the aromatic ring you're adding it to (that's a simple sum-and-divide-by-six average, for each of X, Y and Z). You would then manually add this in as an extra atom to the mol/mol2 file that you feed to molfile_to_params.py (Don't forget to manually add a pseudo bond from the virtual atom to the rest of the molecule - molfile_to_params.py doesn't do well with disconnected structures). For the element symbol, use "X" (special cased to be a "virtual" atom). molfile_to_params.py will spit out a PDB file that you can concatenate onto your protein PDB (you should be doing this anyway, even without the virtual atoms, as molfile_to_params.py will typically rename the atoms, and Rosetta will get confused if the atom names don't match up properly).





okay so esl_fd3bmpnn_extrascores_filt i have at least two different match origins
    1 is the one i showed in presentation that look real good
    another looks almost real good but its got some hydroxyls oriented towards each other and a tyr-glu-ligOH weird situation at the end
L54F94Y40W98Y31
L53F106Y57S101S35
L49M92S38T101N28
L51M98S40T107N13







ORDERING INDIVIDUAL DESIGNS FROM TWIST
    make fasta files of ptn or dna sequence
        optimize for k12 e coli
    use same plasmid backbone as xingjie pet28a+
        - will include n terminal his tag
        - need to specify restriction site Ndel/Xhol
better to specify dna in case you need to order primers/make mutations
    cody has good tool for this
need to know which PO account - ask tanja
twist doesnt add stop codon make sure to add
'''
