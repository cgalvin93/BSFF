'''
time python3 ~/BSFF/spatial_cluster.py /wynton/home/kortemme/cgalvin/dog/fragfbs/dog_Fragment_1_residue_contacts/dog_Fragment_1_SER_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb

import os
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP'] #ONLY TRUE POLAR
l=[]
for i in polar_residues:
    for x in os.listdir():
        if x[-3:]=='pdb':
            if i in x:
                l.append(os.path.join(os.getcwd(),x))

f=open('spatial_clustering.sh','w')
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
cmd='time python3 ~/BSFF/spatial_cluster.py ${tasks[$SGE_TASK_ID]} /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb'
f.write(cmd)
f.write('\nqstat -j "$JOB_ID"')
f.close()

print(len(l))

mkdir cluster_output
qsub -cwd -t 1-11 -l mem_free=8G -o cluster_output -e cluster_output spatial_clustering.sh

scp -r cgalvin@log2.wynton.ucsf.edu:round6/a8s/a8sfragfbs/a8s_Fragment_1_residue_contacts/a8s_Fragment_1_ARG_contact_fuzzball_clustered ~/desktop/argex
'''
#
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from collections import defaultdict
import pandas as pd

def pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
    df=pd.DataFrame(l)
    return df

#convert dataframe to pdb file
def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(dataframe.shape[0]):
        line=''.join(dataframe.iloc[i])
        of.write(line)
    of.close()

#convert ligand pdb to dataframe
def lig_pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:6]=='HETATM':
            d={}
            d['recordname']='ATOM'+'  '
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=' X'
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
        else:
            break
    df=pd.DataFrame(l)
    return df

#return the rmsd between two sets of cartesian coordinates
def calc_rmsd(v1,v2):
    diff = np.array(v1) - np.array(v2)
    N = len(v1)
    return np.sqrt((diff * diff).sum() / N)

#
rfb=sys.argv[1]
ligand_path=sys.argv[2]
#
'''
rfb='/wynton/home/kortemme/cgalvin/dog/fragfbs/dog_Fragment_1_residue_contacts/dog_Fragment_1_SER_contact_fuzzball.pdb'
ligand_path='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb'
'''

#
pdfname=rfb.split('.')[0]+'_clustered.pdf'
pdf = PdfPages(pdfname)

df=pdb_to_df(rfb)
print('df loaded')
#get coordinates of residues
residue_coords_dict=defaultdict(list)
for i in range(df.shape[0]):
        print(str(i))
        resnum=df.iloc[i][6]
        atom_name=str(df.iloc[i][2].split()[0])
        ptn_atom_coords=( float(df.iloc[i][8]), float(df.iloc[i][9]), float(df.iloc[i][10]) )
        if atom_name[0]!='H' and atom_name!='OXT' and atom_name!='C' and atom_name!='CA'and atom_name!='O' and atom_name!='N':
            contact_atom_info=(atom_name,float(ptn_atom_coords[0]),float(ptn_atom_coords[1]),float(ptn_atom_coords[2]))
            if contact_atom_info not in residue_coords_dict[resnum]:
                residue_coords_dict[resnum].append(contact_atom_info)
#create clusters based off of sidechain rmsd for sc contacts
print('\nstarting to cluster contacts...')
print(str(rfb))
cluster_dict=defaultdict(list)
already_clustered=[]
for residue_ in residue_coords_dict.keys():
    print(residue_)
    if residue_ not in already_clustered:
        already_clustered.append(residue_)
        current_cluster=[]
        for residue_2 in residue_coords_dict.keys():
            if residue_!=residue_2 and residue_2 not in already_clustered:
                vs1=[];vs2=[]
                for a,b,c,d in residue_coords_dict[residue_]:
                    for e,f,g,h in residue_coords_dict[residue_2]:
                        if a==e:
                            vs1.append((b,c,d))
                            vs2.append((f,g,h))
                        else:
                            pass
                x=len(vs1)
                y=len(vs2)
                if x!=y:
                    print('different number of atoms between residues '+str(residue_)+' and '+str(residue_2))
                    print('cannot calculate rmsd :(')
                else:
                    try:
                        rmsd=calc_rmsd(vs1,vs2)
                        if rmsd<=0.5:
                            already_clustered.append(residue_2)
                            cluster_dict[residue_].append(residue_2)
                        else:
                            pass
                    except:
                        pass
#okay so next thing is to create pdbs of clusters
##############################
#compile and organize cluster populations
cluster_populations=[]
for key in cluster_dict.keys():
    n_members=len(cluster_dict[key])+1
    cluster_populations.append((key,int(n_members)))
cluster_populations=sorted(cluster_populations, reverse=True, key=lambda nmem: nmem[1])
#plot cluster populations
ncm=[y for x,y in cluster_populations]
clabs=[i+1 for i in range(len(cluster_populations))]
plt.bar(clabs, ncm)
plt.xlabel('Cluster ID')
plt.ylabel('Number of Members')
plt.title('Cluster Population Distribution')
pdf.savefig()
plt.clf()
#make a nice home for stats and cluster pdbs
cluster_results_path=rfb.split('.')[0]+'_clustered'
if len(cluster_populations)>0:
    os.makedirs(cluster_results_path,exist_ok=True)
    n_contacts_check=[]
    fullligdf=lig_pdb_to_df(ligand_path)
    #export cluster fuzzballs with full ligand
    rns=[]
    for i in df['resnum']:
        rns.append(int(i))
    for a,b in cluster_populations:
        cluster_df=pd.DataFrame()
        if b>=0:
            n_contacts_check.append(b)
            resnumbers=[int(i) for i in cluster_dict[a]]
            resnumbers.append(int(a))
            for resnum in resnumbers:
                rows=[i for i,e in enumerate(rns) if e == resnum]
                res_df=df.iloc[min(rows):max(rows)+1]
                cluster_df=cluster_df.append(res_df,ignore_index=True)
            cluster_df=cluster_df.append(fullligdf,ignore_index=True)
            resnamee=str(res_df.iloc[0][4])
            ofilename=str(a.strip())+'_'+resnamee+'_cluster_'+str(b)+'.pdb'
            df_to_pdb(cluster_df,os.path.join(cluster_results_path,ofilename))
        else:
            continue
    n_contacts_check_sum=sum(n_contacts_check)
    print('there are '+str(n_contacts_check_sum)+' contacts after clustering')
else:
    print('there are '+str(len(cluster_populations))+' clusters')



pdf.close()







'''

scp -r cgalvin@log2.wynton.ucsf.edu:round6/a8s/a8sfragfbs/a8s_Fragment_1_residue_contacts/a8s_Fragment_1_LYS_contact_fuzzball_clustered ~/desktop/lysclust


DO SO I CAN SPLIT EACH INPUT RESIDUE FUZZBALL ACROSS DIFFERENT CORES
THEN DO SOME EXPLICIT MODELING (MC MOTIF ASSEMBLY)
AND TRY TO MATCH THOSE
    maybe just similar approach as adding np instead of doing mc even
    i can just do random combos of favorable high freq interactions that
    have energy below certain threshold
REU VS POPULATION FOR EACH RESIDUE TYPE
    maybe rosetta get np res right but no polar?
    or just certain res but not others?
FIND POLAR RESIDUE PLACEMENTS THAT CAN SATISFY 2 LIG POL ATOMS SIMULTANEOUSLY?
    'special hb contact' would include things like arg carb bidentate too
FIRST ID BB POSITIONS IN BS FREE TO HBOND? build more sophisticated hb network from there?
USE MPNN GOOD SCORE SEQUENCE PROFILE TO INFORM DESIGN
ALSO TRY ENZDES METHOD INSTEAD OF JUST FD

IDEALIZED CST THING -
    maybe more sophisticated ways to do res stats n shit? idk
    or get ideal/std vals for each param from res-frag stats

SHOULD I BE DOING DIFFERENT LIGAND CONFORMERS?
    hydroxyls can rotate n shit

2 PATHS
    spatial clustering + mc
        special rot?
    idealized hbond csts using res stats to select for each frag
        add np csts?
        special rot?

MANY SPATIAL CLUST JOBS STILL RUNNING AFTER 21 HRS, PROCEED WITH IDEALIZED

'''
import json
with open('dog_polar_contacts_raw.json','r') as f:
    df=json.load(f)

from collections import defaultdict


for k in df.keys():
    pops=defaultdict(list)
    for k2 in df[k].keys():
        pop=0
        for k3 in df[k][k2].keys():
            pop+=len(df[k][k2][k3])
        pops[k2]=pop

for i in df['Fragment_1']['dog_Fragment_1_LYS_contact_fuzzball.pdb']['resdonors'].keys():
    print(df['Fragment_1']['dog_Fragment_1_LYS_contact_fuzzball.pdb']['resdonors'][i])




















#pass as residue type either Hpol for all sc or HNbb for bb gly
#this function i will use just to add np cst at end, polar cst blocks are built manually
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



import random
import os
#
n_desired_motifs=10
lig_resname='dog'
#define the 3 constrained atoms for each fragment
relevant_atoms={'Fragment_1':['yes','O5','C23','C22'],
                'Fragment_2':['n','O1','C1','C3'],
                'Fragment_3':['n','O2','C6','C8'],
                'Fragment_4':['n','O3','C19','C18']}
#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
polar_residues=['LYS','ARG','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY']
problematic_residues=['LYS','ARG','ASP','GLU','GLN','ASN','HIS','TRP']
hydroxyl_fragments=['Fragment_2','Fragment_3','Fragment_4']
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

#for each fragment, create a cst block
motifs=[]
while len(motifs)<n_desired_motifs:
    problem_res=0
    ngly=0
    current_motif=[]
    for frag in relevant_atoms.keys():
        sp2_acc=relevant_atoms[frag][0]
        lig_atoms=relevant_atoms[frag][1:]
        res_resname=random.choice(polar_residues)
        if res_resname=='GLY':
            res_atoms='HNbb,'
            ngly+=1
        else:
            res_atoms='Hpol,'
        if res_resname in problematic_residues:
            problem_res+=1
        if sp2_acc=='yes':
            cstblock=generate_single_constraint_block_base(distanceAB=1.95,angleA=120,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                           residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                           distance_tolerance=0.15, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=20, torsion_AB_tolerance=180, torsion_B_tolerance=180,
                                                           torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                           distance_constraint_sample_number=1)
        else:
            cstblock=generate_single_constraint_block_base(distanceAB=1.95,angleA=120,angleB=180,torsionA=120,torsionAB=0,torsionB=0,
                                                           residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                           distance_tolerance=0.15, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=180, torsion_AB_tolerance=180, torsion_B_tolerance=180,
                                                           torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                           distance_constraint_sample_number=1)
            # if frag in hydroxyl_fragments:
            #     cstblock=generate_single_constraint_block_base(distanceAB=1.95,angleA=120,angleB=180,torsionA=130,torsionAB=0,torsionB=0,
            #                                                    residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
            #                                                    distance_tolerance=0.15, angle_A_tolerance=10, angle_B_tolerance=15,
            #                                                    torsion_A_tolerance=60, torsion_AB_tolerance=180, torsion_B_tolerance=180,
            #                                                    torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
            #                                                    distance_constraint_sample_number=1)
            # else:
            #     cstblock=generate_single_constraint_block_base(distanceAB=1.95,angleA=120,angleB=180,torsionA=0,torsionAB=0,torsionB=0,
            #                                                    residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
            #                                                    distance_tolerance=0.15, angle_A_tolerance=10, angle_B_tolerance=10,
            #                                                    torsion_A_tolerance=180, torsion_AB_tolerance=180, torsion_B_tolerance=180,
            #                                                    torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
            #                                                    distance_constraint_sample_number=1)
        current_motif.append(cstblock)
    if problem_res<=1:
        if ngly<=1:
            if current_motif not in motifs:
                motifs.append(current_motif)

polar_cst_paths='t2_idealized_hbonds'
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
d 1.8-2.1
psi 110-130
phi 160-180
chi 180 +/20, periodicity 2



for the protein donating to the ligand
    angle a
    angle b
    tor a (sp2 acc)
    -------
    no tor b, tor ab

tor a = Res1:Atom3 - Res1:Atom2 - Res1:Atom1 - Res2:Atom1
'torsion_AB' is the dihedral Res1:Atom2 - Res1:Atom1 - Res2:Atom1 - Res2:Atom2
'torsion_B' is the dihedral Res1:Atom1 - Res2:Atom1 - Res2:Atom2 - Res2:Atom3.

lig is res 1, hydroxyl 3 toms = o c c

the orientation of the ligand h relative to donating h is descvribed by tor a
using the lig h to define this torsion, the favorable range is about
20 - -120 = 20 - 240 --> 130 +/- 110
    or actually, when using the actual ptn donating h, its torsion ab i think,
    which is the rotation about the hydrogen bond itself
            OKAY YA KNOW WHAT I HAVE NO FUCKING IDEA BUT IM GONNA TEST IT
ALRIGHT YEA ITS LOOKIN LIKE IT IS THE TOR A but i still gotta figure out ideal angles and tol



if resdonor==True:
    angle_A=ang_ind[psibin]
    angle_B=ang_ind[phibin]
    try:
        torsion_A=tor_ind[int(contact_data[10])-1]
        torsion_A_tolerance=10
        torsion_A_constraint_sample_number=1
    except:
        torsion_A=0
        torsion_A_tolerance=180
        torsion_A_constraint_sample_number=torsion_constraint_sample_number
    torsion_AB=0
    torsion_B=0
    angle_A_tolerance=10
    angle_B_tolerance=10
    torsion_AB_tolerance=180
    torsion_B_tolerance=180

elif resdonor==False:
    angle_A=ang_ind[phibin]
    angle_B=ang_ind[psibin]
    try:
        torsion_B=tor_ind[int(contact_data[10])-1]
        torsion_B_tolerance=10
        torsion_B_constraint_sample_number=1
    except:
        torsion_B=0
        torsion_B_tolerance=180
        torsion_B_constraint_sample_number=torsion_constraint_sample_number
    torsion_AB=0
    torsion_A=0
    angle_A_tolerance=10
    angle_B_tolerance=10
    torsion_AB_tolerance=180
    torsion_A_tolerance=180

d 1.8-2.1
psi 110-130
phi 160-180
chi 180 +/20, periodicity 2

cstblock=generate_single_constraint_block_base(distance_AB,angle_A,angle_B,torsion_A,torsion_AB,torsion_B,
                                               res_resname,lig_resname,lig_atoms,res_atoms,
                                               distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                                               torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                               torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                               distance_constraint_sample_number=0)

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

polar_residues=['LYS','ARG','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY']


  TEMPLATE::   ATOM_MAP: 2 atom_name: 2HZ NZ CE
  TEMPLATE::   ATOM_MAP: 2 residue3: LYS

  TEMPLATE::   ATOM_MAP: 2 atom_name: HG1 OG1 CB
  TEMPLATE::   ATOM_MAP: 2 residue3: THR

  TEMPLATE::   ATOM_MAP: 2 atom_name: HE1 NE1 CD1
  TEMPLATE::   ATOM_MAP: 2 residue3: TRP

  TEMPLATE::   ATOM_MAP: 2 atom_name: H N CA
  TEMPLATE::   ATOM_MAP: 2 residue3: GLY

  TEMPLATE::   ATOM_MAP: 2 atom_name: 1HD2 ND2 CG
  TEMPLATE::   ATOM_MAP: 2 residue3: ASN

   TEMPLATE::   ATOM_MAP: 2 atom_name: 1HE2 NE2 CD
  TEMPLATE::   ATOM_MAP: 2 residue3: GLN

ARG NH2 CZ NE
tyr
'''




'''
                            MATCHING
'''
os.chdir(polar_cst_paths)
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# import os
#now submit for matching with lucs strc
#
lig_name='dog'
# lig_name='38e'
#
# paramspath='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
paramspath='/wynton/home/kortemme/cgalvin/dog.params'
# paramspath='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
#
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p_1np'
# allmotifsdir='/wynton/home/kortemme/cgalvin/dog/t2_idealized_hbonds'
allmotifsdir=os.getcwd()
#
shellfile_suffix='_dog_ideal.sh'#
###########################################################################
############################################################################
############################################################################
############################################################################
scaffold_path='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'
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
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
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
    idx=[j for j in range(0,len(motif_paths[key]),10)]
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
50
n scaffolds
418

~/main/source/bin/match.linuxgccrelease -s /wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/relaxed_5tpj_143446_design_6_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001.pdb -extra_res_fa /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params -match:geometric_constraint_file /wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/cst_10.cst -match::scaffold_active_site_residues  /wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/relaxed_5tpj_143446_design_6_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001.pos -match:output_format PDB -match:match_grouper SameSequenceGrouper -match:consolidate_matches -match:output_matches_per_group 1 -use_input_sc -ex1 -ex2 -extrachi_cutoff 0 -enumerate_ligand_rotamers false -match::lig_name dog



qsub -cwd -t 1-418 -l mem_free=40G -o cluster_out -e cluster_out _1_dog_ideal.sh

scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/matchex ~/desktop/matchex
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('_1_a8s_3p1np.sh')
#
for j in jobss:
    cmd='qsub -cwd -t 1-418 -l mem_free=35G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''
scp -r cgalvin@Log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/t3_idealized_hbonds/matchex ~/desktop/matchex3

i think this torsion does actually need to take on discrete values to count as a good hbond,
orienting the oxygen lone pair towards the hydrogen

-102.7
133.3
wrong h = -22.5, -146.6

if i get params file for multiple conformers will it allow the hydroxyls to rotate away
${ROSETTA}/main/source/src/python/apps/public/molfile_to_params.py -n LIG -p LIG --conformers-in-one-file LIG_conf.sdf

gonna try to do matching and design with the conformer containing params file
'''



























import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
len(l)

'''
Out[1]: 4362
'''

#clean matches and put in own directory
# import os
os.makedirs('design',exist_ok=True)
# pdbs=[i for i in os.listdir() if i[:2]=='UM']
pdbs=l
print(str(len(pdbs)))
bad=[]
c=1
for pdb in pdbs:
    print(str(c))
    c+=1
    with open(pdb,'r') as f:
        l=[line for line in f.readlines()]
        if len(l)<50:
            bad.append(pdb)
print(str(len(bad)))
'''
31

'''
pdbs2=[i for i in pdbs if i not in bad]
firsts=pdbs2
#quarantining this for now cus my files were getting fucked up
#try one more time actually to see if it works after deleting bad files
l=[]
for pdb in firsts:
    s='egrep "^ATOM|HETATM|REMARK 666" '+pdb+' >design/'+pdb
    # print(s)
    # os.system(s)
    l.append(s)
for i,c in enumerate(l):
    os.system(c)
    print(str(i))



os.chdir('design')





#make a resfile for each match structure
import os
#################################################################################
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'

sfname='srmr.sh'
#################################################################################
os.makedirs('cluster_out',exist_ok=True)
matches=[i for i in os.listdir() if i[-3:]=='pdb']
# idx=[j for j in range(0,len(matches),100)]
# if len(matches) not in idx:
#     idx.append(len(matches))
c=1
# for ei, ix in enumerate(idx[:-1]):
#     selection=matches[ix:idx[ei+1]]
#     sf=open('_'+str(c)+sfname,'w')
#     sf.write('#!/bin/bash')
#     sf.write('\n')
#     sf.write('source ~/anaconda3/etc/profile.d/conda.sh')
#     sf.write('\n')
#     sf.write('conda activate pyr37')
#     sf.write('\n')
#     sf.write('tasks=(0\n')
#     for match in selection[:-1]:
#         sf.write('       '+str(match)+'\n')
#     sf.write('       '+selection[-1]+')')
#     sf.write('\n')
#     cmd='time python ~/BSFF/tools/match_resfile.py ${tasks[$SGE_TASK_ID]} '+prm1+' 5.0 no'
#     sf.write(cmd)
#     sf.write('\nqstat -j "$JOB_ID"')
#     sf.close()
#     c+=1
sf=open('_'+str(c)+sfname,'w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('source ~/anaconda3/etc/profile.d/conda.sh')
sf.write('\n')
sf.write('conda activate pyr37')
sf.write('\n')
sf.write('tasks=(0\n')
for match in matches[:-1]:
    sf.write('       '+str(match)+'\n')
sf.write('       '+matches[-1]+')')
sf.write('\n')
cmd='time python ~/BSFF/tools/match_resfile.py ${tasks[$SGE_TASK_ID]} '+prm1+' 5.0 no'
sf.write(cmd)
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
c+=1
print(len(matches))
'''
print(len(matches))
4331
qsub -cwd -t 1-4331 -l mem_free=1G -o cluster_out -e cluster_out _1srmr.sh

cst_10.cst
'''

# # changing cst files to match the pdb of matches names
import os
################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np'
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np'
cstdir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds'


################################################
for i in matches:
    ################################################
    cstid='cst_'+i.split('_')[-2]+'.cst'
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







# FASTDESIGN ON MATCHES
# setup so that 1 params for all
# resfile and cst must have same name as pdb
import os
########
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
sfname='dog_ideal_fd.sh'
ndesignseachmatch='10'
outputdirectory='dog_ideal_fd'
scorefile_name='dog_ideal_fd.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(1):
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
         '-parser:protocol','~/BSFF/tools/FD_jump1_nofilters.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights','ligand',
         '-score:set_weights','hbond_bb_sc','2.0',
         '-score:set_weights','hbond_sc','2.0',
         '-load_PDB_components','False',
         '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name]
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1

'''
print(len(matches))
4331

qsub -cwd -t 1-2214 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh
qsub -cwd -t 1-18 -l mem_free=16G -o cluster_output -e cluster_output r638efd.sh
'''


# import os
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4331 -l mem_free=8G -o cluster_output -e cluster_output '+j
    os.system(cmd)


























#analysis of scores from json file
sfname='dog_ideal_fd.json'
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

plot_dists(terms,scores,'dog_idel_fd.pdf')

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
f4=return_filtered(f3,'contact_molsurf','>',180)
plot_dists(terms,f4,'dogidealfilt.pdf')
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
'''
13302
691
691
648
'''
import os
filtered_strc=[]
for d in f4:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
os.system('mv *.pdf filtered/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered ~/desktop/dogidealfilt

WEIRDNESS WITH HYDROXYL GROUPS FACING EACH OTHER
WANT TO GET BOLTZ/DDG SCORES AND AF CONFIDENCE/RMSD
'''

#run the extra scoring
paramspath='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('score_dog_filt.sh','w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('tasks=(0\n')
for match in pdbs[:-1]:
    sf.write('       '+match+'\n')
sf.write('       '+pdbs[-1]+')')
sf.write('\n')
sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]} -parser:protocol /wynton/home/kortemme/cgalvin/design_and_analysis/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+paramspath+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
#
print(len(pdbs))
'''
qsub -cwd -t 1-648 -l mem_free=8G score_dog_filt.sh
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

plot_dists(terms,scores,'dog_idel_fd_filt_extrascores.pdf')
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
f2=return_filtered(f1,'boltz','<',-0.25)
f3=return_filtered(f2,'ddg','<',-23.0)
print(len(f1))
print(len(f2))
print(len(f3))
'''
135
57
42
'''
import os
filtered_strc=[]
for d in f3:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered2',exist_ok=True)
for i in filtered_strc:
    os.system('cp /wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/'+i[:-5]+'.pdb filtered2/'+i+'.pdb')
os.system('mv *.pdf filtered2/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/filtered2 ~/desktop/dogidealfiltered2
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
mkdir cf
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


import os
#
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
shf_prefix='_dogidealcf'
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
allfastasdir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/filtered2/cf/fastas'
# allfastasdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas'
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
nonaf2desdir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/filtered2'

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
# json.dump(data_allstrc,open('a8s3p1np_filt_af2.json','w'))
json.dump(data_allstrc,open('dogideal_filt_af2.json','w'))
##################################################################
##################################################################
##################################################################

# import numpy as np
#id designs with plddt over threshold
plddt_threshold=85.0
carmsd_threshold=2.0
aplddts=[]
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        plddts.append(cplddt)
        carmsds.append(carmsd)
    aplddt=np.mean(plddts)
    acarmsd=np.mean(carmsds)
    if aplddt>=plddt_threshold:
        if acarmsd<=carmsd_threshold:
            aplddts.append((aplddt,key,k2,acarmsd))

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
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/filtered/filtered2/cf/fastas/af2filtered ~/desktop/dogidealaf2filt

----------------
MAYBE I CAN USE THE INFORMATION IN MPNN LOW SCORING SEQUENCE PROFILE FOR A GIVEN
MATCHED BACKBONE TO INFORM DESIGN SUCH THAT IT IS BIASED TOWARDS
THE TOLERATED SEQUENCE PROFILE OF THAT BACKBONE, PRODUCING DESIGNS
WHICH ARE MORE LIKELY TO HAVE HIGH AF CONFIDENCE AND LOW RMSD TO INPUT
    for a given match, take it and use mpnn to design the residues which will be
    flagged as designable in resfile for rosetta fd/cm/whatever,
    of the resultant sequences, identify those with good mpnn scores below some threshold,
    analyze this favorable mpnn sequence profile to identify favored residues at each position
        layer deeper and identify something akin to coupled mutations even?idk how to get deeper with it
        really but i feel like i could
    carry out fd with a resfile that allows only the favored residues at each designable position
i want then that the matched backbones, when redesigned by mpnn, and the design sequences submitted to
alphafold, i get consistenly low rmsd between input and output, that is mpnn considers the
bb itself good
    could keep fucking with the current bb set iteratively until i get an avg rmsd of
    favorable mpnn outputs that is below a certain threshold, say 1 angstrom
---------------
PROBABLY I CAN FIX ONE MORE OF THE MATCHER TORSIONS TO CONSTRAIN THE HYDROXYL-HYDROXYL
HBONDS IN THE PROPER ORIENTATION
    struggling with this, but trying to include ligand conformers in matching+design
    and see if that will allow hydroxyls to like, turn the right way?
---------------
SUSPECTING THAT FD DOESNT GENERATE ENOUGH SEQUENCE DIVERSITY
    analuyze design seqs to confirm
    run cm (w gp too cus fuck it why not) in addition to fd
    somehow special rot?
        still thing of just np, is it worth, does rosetta need the help, how can
            i demonstrate, etc.....


TRYING COUPLEDMOVES WITH LIGAND CONFORMERS

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/t2_idealized_hbonds/UM_1_H4S11S76G73_1_relaxed_5tpj_41695_design_3_unrelaxed_model_1_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_cst_2_1.pdb ~/desktop/sdf.pdb
/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design


'''

#




# coupledmoves
# setup so that 1 params for all
# resfile and cst must have same name as pdb
import os
########
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
# prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
prm1='/wynton/home/kortemme/cgalvin/dog.params'
sfname='dog_ideal_cm.sh'
ndesignseachmatch='5'
outputdirectory='dog_ideal_cm'
scorefile_name='dog_ideal_cm.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
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
         '-parser:protocol','~/BSFF/tools/cm_standard_jump1.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-mute protocols.backrub.BackrubMover',
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights','ligand',
         '-restore_pre_talaris_2013_behavior',
         '-score:set_weights','hbond_bb_sc','2.0',
         '-score:set_weights','hbond_sc','2.0',
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
         '-coupled_moves::exclude_nonclashing_positions true',
         '-coupled_moves::backbone_mover kic',
         '-coupled_moves::kic_perturber walking']
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1
'''
print(len(matches))
4331

qsub -cwd -t 1-2214 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh
'''
# import os
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4331 -l mem_free=8G -o cluster_output -e cluster_output '+j
    os.system(cmd)






'''
get filterered strc from this run and see if they
still have weird oh-ph hbond orientations or if using the conformer params allows this
to be rescued
    if so, gonna try mpnn biased design and see if i can improve af2 conf of final designs

'''


#analysis of scores from json file
sfname='dog_ideal_cm.json'
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

plot_dists(terms,scores,'dog_ideal_cm.pdf')

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
f3=return_filtered(f2,'sc','>',0.6)
# f3=return_filtered(f2,'hbtolig','>',3.0)
# f4=return_filtered(f3,'contact_molsurf','>',180)
plot_dists(terms,f3,'dogidealfiltcm.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
# print(len(f4))
'''
13302
691
691
648
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
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_cm/filtered ~/desktop/dogidealfiltcm


OKAY INTERESTING RESULTS, IT LOOKS LIKE IT DID
IN FACT FIX THE HYDROXYL ORIENTATION THING, AND IN ADDITION IVE GOT TONNSSSSS OF
POLAR RESIDUES AND MUCH DENSER HBOND NETWORKS, BUT ITS TOO POLAR
AND IVE GOT BAD SC SCORES N SHIT,
    try fd with the conformer params
    try not upweighting the hbond terms maybe, just ptn-lig interaxns?
    try the mpnn thing
        id matches w/ good burial first
    try with special rotamers? this gonna be even more complicated with the extra conformers
ALSO IF IM GONNA RUN WITH THE EXTRA CONFORMER DESIGN I SHOULD CHECK THE ENERGY OF THE
CONFORMERS AND ONLY USE THOSE BELOW SOME CUTOFF OR SOMETHING
    maybe i can add something to the script to calculate the energy of
    the input conformer (which should be lowest energy)
    and then only spiut out additional conformers which have an energy within some range
    of this input value

LEMME TRY FD WITH ADDITIONAL CONFORMERS
'''


# FASTDESIGN ON MATCHES w extra conformers
# setup so that 1 params for all
# resfile and cst must have same name as pdb
import os
########
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
# prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
prm1='/wynton/home/kortemme/cgalvin/dog.params'
sfname='dog_ideal_fd2.sh'
ndesignseachmatch='10'
outputdirectory='dog_ideal_fd2'
scorefile_name='dog_ideal_fd2.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(1):
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
         '-parser:protocol','~/BSFF/tools/FD_jump1_nofilters.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights','ligand',
         '-restore_pre_talaris_2013_behavior',
         # '-score:set_weights','hbond_bb_sc','2.0',
         # '-score:set_weights','hbond_sc','2.0',
         '-load_PDB_components','False',
         '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name]
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1

'''
print(len(matches))
4331

qsub -cwd -t 1-2214 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh
qsub -cwd -t 1-18 -l mem_free=16G -o cluster_output -e cluster_output r638efd.sh
'''
# import os
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4331 -l mem_free=8G -o cluster_output -e cluster_output '+j
    os.system(cmd)








#analysis of scores from json file
sfname='dog_ideal_fd2.json'
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

plot_dists(terms,scores,'dog_ideal_fd2.pdf')

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
# f3=return_filtered(f2,'sc','>',0.6)
f3=return_filtered(f2,'hbtolig','>',4.0)
f4=return_filtered(f3,'contact_molsurf','>',180)
plot_dists(terms,f3,'dogidealfiltfd2.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
'''
43166
12644
574
574
539

i should add the remove match csts -> relax/min to end of fd script as well
'''
import os
filtered_strc=[]
for d in f4:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
os.system('mv *.pdf filtered/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd2/filtered ~/desktop/dogidealfd2filt

this is without torsion a constraint
using extra ligand conformers
hb weights x2
lig weight x2
ligand scorefunction

i wanna try the mpnn biased design next, but let me see whats up with the boltz scores/
af confidence of this set
'''

#run the extra scoring
paramspath='/wynton/home/kortemme/cgalvin/dog.params'
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('score_dog_filt.sh','w')
sf.write('#!/bin/bash')
sf.write('\n')
sf.write('tasks=(0\n')
for match in pdbs[:-1]:
    sf.write('       '+match+'\n')
sf.write('       '+pdbs[-1]+')')
sf.write('\n')
sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]} -parser:protocol /wynton/home/kortemme/cgalvin/design_and_analysis/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+paramspath+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
#
print(len(pdbs))
'''
qsub -cwd -t 1-539 -l mem_free=8G score_dog_filt.sh

ideally i want some way to only upqeight hbond terms for residues which are already
hbonding with the ligand
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

plot_dists(terms,scores,'dog_idel_fd2_filt_extrascores.pdf')
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
f2=return_filtered(f1,'boltz','<',-0.25)
f3=return_filtered(f2,'ddg','<',-23.0)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
'''
135
57
42
'''
import os
filtered_strc=[]
for d in f3:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered2',exist_ok=True)
for i in filtered_strc:
    os.system('cp /wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd/'+i[:-5]+'.pdb filtered2/'+i+'.pdb')
os.system('mv *.pdf filtered2/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd2/filtered/filtered2 ~/desktop/dogidealfd2filt2

FIX DIHYDROXYL
    maybe if i use the ligand hydroxyl h as the first atom, there will be a
    way to describe the hbonds correctly
I MAKE 2 ETHANOLS IN AVOGADRO AND MINIMIZE ENERGY TO SEE IDEAL HBOND GEOMETRY
    lignd is always res1, lemme look at the cst dihedrals and angles

right way
    distance 1.744
    angle a 119
    angle b 175
    torsion a   139.2
    torsion ab  -71.2
    torsion b   -19
wrong way
    distance
    angle a     119
    angle b     175
    torsion a   139.2
    torsion ab  -71.2
    torsion b   -19
                OKAY THIS CONFIRMS THAT THERE IS NO WAY TO SPECIFY THE
                POSITION OF THE LIGAND HYDROGEN THIS WAY, IT NEEDS TO BE INCLUDED
                AS AN ATOM

okay what if i relax the structures with extra conformers annnnd using a
scorefunction with no hbond terms (i guess just faatr/rep)
    the problem  is this could fuck up my other hbond interactions right?because it
    wont like the overlapping of the vdw radii
SHIIIIT WELL IF I INCLUDE THE H AS AN AT0M IN CST HOW CAN I SPECIFY
THE PARAMTER RANGES THAT CORRESPOND TO FAVORABLE HBONDS

from rossettacom mons
"When you generate the params file for your ligand with the molfile_to_params.py
script, it should be smart enough to recognize rotatable hydrogen-bondable hydr
ogens, and add lines to sample those hydrogen positions in the params file. (Yo
u can look for "PROTON_CHI" lines in the params file for the rotatable hydrogen
s.) The conformers in the conformer files only need to be the heavy atom confor
mers."
    -https://www.rosettacommons.org/node/4024

is this possibly a consequence of using ligand weights? should i use ref15?

OKAY SO I SHOULDNT NEED THE CONFORMERS PARAMS
MAYBE TRY NOT UPWEIGHTING HBONDS
INCLUDING PRETALARISRESTORE OPTION WITH FD (was forgetting this before so idk what was happening)
INCLUDE A NO CST FASTRELAX BEFORE APPLYING FILTERS/METRICS
    if i wanna promote more hbonding i will try to use hbnet or 3bop or something
'''

# FASTDESIGN ON MATCHES w extra conformers
# setup so that 1 params for all
# resfile and cst must have same name as pdb
import os
########
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
# prm1='/wynton/home/kortemme/cgalvin/dog.params'
sfname='dog_ideal_fd3.sh'
ndesignseachmatch='10'
outputdirectory='dog_ideal_fd3'
scorefile_name='dog_ideal_fd3.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(1):
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
         '-parser:protocol','~/BSFF/tools/FD_jump1_nofilters.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights','ligand',
         '-restore_pre_talaris_2013_behavior',
         # '-score:set_weights','hbond_bb_sc','2.0',
         # '-score:set_weights','hbond_sc','2.0',
         '-load_PDB_components','False',
         '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name]
    sf.write((' ').join(cmd))
    sf.write('\nqstat -j "$JOB_ID"')
    sf.close()
    c+=1
print(len(matches))
# import os
# l=[i for i in os.listdir() if i[-5:]=='fasta']
# l1=[i for i in os.listdir() if i[-7:]=='low.pdb']
# l2=[i for i in os.listdir() if i[-8:]=='last.pdb']
# for i in l:
#     os.system('rm '+i)
# for i in l1:
#     os.system('rm '+i)
# for i in l2:
#     os.system('rm '+i)
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-4331 -l mem_free=8G -o cluster_output -e cluster_output '+j
    os.system(cmd)

'''

~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s UM_3_S66R62T20Y40_1_relaxed_5tpj_74354_design_7_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_cst_40_1.pdb -parser:protocol ~/BSFF/tools/FD_jump1_nofilters.xml -parser:view -run:preserve_header -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params -resfile UM_3_S66R62T20Y40_1_relaxed_5tpj_74354_design_7_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_cst_40_1.resfile -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 10 -score::weights ligand -restore_pre_talaris_2013_behavior -load_PDB_components False -enzdes:cstfile UM_3_S66R62T20Y40_1_relaxed_5tpj_74354_design_7_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_cst_40_1.cst -out:path:all dog_ideal_fd3 -out:prefix fd1 -scorefile_format json -out:file:scorefile dog_ideal_fd3.json


scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd3/matchex ~/desktop/dogidealfd3
'''
#analysis of scores from json file
sfname='dog_ideal_fd3.json'
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

plot_dists(terms,scores,'dog_ideal_fd3.pdf')

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
# f3=return_filtered(f2,'sc','>',0.6)
f3=return_filtered(f2,'hbtolig','>',4.0)
f4=return_filtered(f3,'contact_molsurf','>',180)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))

'''
39496
14240
698
564
548
'''
plot_dists(terms,f4,'dogidealfiltfd3.pdf')
import os
filtered_strc=[]
for d in f4:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
os.system('mv *.pdf filtered/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd3/filtered ~/desktop/dogidealfd3filt
'''

#run the extra scoring
paramspath='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('score_dog_filt.sh','w')
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
qsub -cwd -t 1-548 -l mem_free=4G score_dog_filt.sh
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

plot_dists(terms,scores,'dog_idel_fd3_filt_extrascores.pdf')


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
f3=return_filtered(f2,'hbtolig','>',4.0)
f4=return_filtered(f3,'contact_molsurf','>',180)
f5=return_filtered(f4,'boltz','<',-0.3)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
'''
469
412
402
333
333
11
'''
plot_dists(terms,f5,'dog_idel_fd3_filt_extrascoresfilt.pdf')

import os
filtered_strc=[]
for d in f5:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered2',exist_ok=True)
for i in filtered_strc:
    os.system('cp /wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd3/'+i[:-5]+'.pdb filtered2/'+i+'.pdb')
os.system('mv *.pdf filtered2/')

'''

scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/dog_ideal_fd3/filtered/filtered2 ~/desktop/dogidealfd3extrafilt2
'''


























'''
RUNNING MPNN ON MATCHES TO GENERATE SEQUENCE PROFILE

THEN IM GONNA USE 3BOP WHEN I DESIGN ALSO




OKAY, FIRST I NEED TO IDENTIFY MATCHES WITH GOOD LIGAND BURIAL,
ROUGHLY, I WILL DO THIS USING THE DSASA FILTER


NOW RUN MPNN WITH FIRST SHELL RESIDUES FROZEN
maybe the simplest is just to do all residues within a certain distance cutoff
of ligand, that for sure has to include the hbond residues,
hmmm and then maybe check residue energy? shit but then i need to explicitly check hbonds


/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered
/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered


'''
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
ofjsonname='mpnn_params.json'
################
import os
from pyrosetta import*
import json
init('-load_PDB_components False')
##########
# sf = ScoreFunction()
# from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
# sf.set_weight(fa_atr, 1)
# sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)
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
#
# for pdb in ttfdata.keys():
for pdb in pdbs:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    dsasa_val=dsasa_filter.report_sm(p)
    if dsasa_val>0.8:
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
'''
In [3]: print(len(list(tfdata.keys())))
   ...:
2625



MPNN WITH CONSTRAINTS
'''

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
ofjsonname='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/mpnn_params.json'
#########
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn'
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn'
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/mpnn'
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


#making the shellfile scripts
# import json
# import os
#########
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'

# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn_params.json'
#########
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn'
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn'
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/mpnn'
directory_prefix='mpnndesign'
###############
shf_prefix='fdogiltdesmpnn'
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    p=[i for i in os.listdir(id) if i[-3:]=='pdb']
    pid=p[0]
    # tf=[str(i[0]) for i in tfdata[pid]]
    tf=[str(i) for i in tfdata[pid]]
    tfs=' '.join(tf)
    submitlist.append((id,os.path.join(id,outputdirname),tfs))
#
c=1
for id,od,tfs in submitlist:
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
    of.write('fixed_positions="'+tfs+'"')
    of.write('\n')
    of.write('python /wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains')
    of.write('\n')
    of.write('python /wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --position_list "$fixed_positions"')
    of.write('\n')
    cmdl=['python',
    '/wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/protein_mpnn_run.py',
    '--jsonl_path $path_for_parsed_chains',
    '--out_folder $output_dir',
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
qsub -cwd -l mem_free=16G -o mpnncout -e mpnncout filtdesmpnn_1_.sh

'''
mkdir mpnncout
# import os
jobss=[i for i in os.listdir() if shf_prefix in i]
#
# jobss.remove('filtdesmpnn_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)



'''
consolidating mpnn results
        some didnt work
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/mpnn'
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
1981
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

now i do some analysis of the sequences

I GUESS I WANNA COUNT THE NUMBER OF EACH AMINO ACID AT EACH POSITION FOR DESIGNABLE POSITIONS
'''
import os
from pyrosetta import*
import json
init('-load_PDB_components False')

#########################
prm1='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
ofjsonname='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/mpnn_params.json'
pdbsdir='/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/design/'
#########################

with open(ofjsonname,'r') as f:
    tfdata=json.load(f)


pdbs=[i for i in os.listdir() if i[-3:]=='pdb']


for strc in list(add.keys())[:1]:
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
    for des_ind in frozen_res:
        positional_seq_list=[]
        for mpnn_data_entry in add[strc]:
            score=mpnn_data_entry[0]
            if score<=1.0:
                seq=mpnn_data_entry[2]
                positional_seq_list.append(seq[des_ind-1])
        strc_seqs[des_ind]=positional_seq_list


'''
OKAY WEVE GOT A PROBLOEM IT APPEARS THAT PROTEINMPNN DESIGNED AT EVERY POSITION.....

the make fied dict helper script doesnt appear to be working
an example of an output dict is:

{"UM_5_S57K17T94T13_1_relaxed_5tpj_322479_design_3_unrelaxed_model_3_rank_1_0001_design_1_unrelaxed_model_1_rank_1_0001_cst_5_1": {"A": []}}

maybe i can simply make these dictionaries myself and then change shellscripts to
get rid of the helper script execution line and simply provide the
path to the fixed positions dict instead....



'''















'''

                    NEXT ROUND

FURTHER IDEALIZE NTF2 SCAFF LIBRARY
NONPOLAR FRAGS FOR MATCHING
IDEALIZED HBONDS FOR INDIVIDUAL RESATOM:LIGATOM PAIRS
MPNN INFORMED DESIGN (pssm)
GENPOT
AUTOMATE RETRIEVAL/ANALYSIS/COMPARISON OF PDB BINDERS FOR TARGET
    and also edit bsff scripts generally, may have rotation student soon
examples of same scaff diff match to show different binding modes in same site


so i guess the chronological move would be to
    idealize ntf2s
    restart bsff run with nonpolar frags included
    match hybrid motifs which include specific idealized hbonds for atom:atom pairs
    pssm/genpot/cm etc. based design
    analysis and comparison with binders in pdb ...
    
'''
