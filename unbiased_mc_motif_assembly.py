#time ipython mc_motif_assembly.py fuzzball_for_assembly fuzzball_for_assembly_json params
#time ipython ~/desktop/prj/bsff/bsff_scripts/mc_motif_assembly.py fuzzball_for_assembly.pdb fuzzball_for_assembly_frequencies.json Inputs/Rosetta_Inputs/dex.params

#time python ~/desktop/prj/bsff/bsff_scripts/mc_motif_assembly.py -d a8s
#time python ~/desktop/prj/bsff/bsff_scripts/mc_motif_assembly.py -d 'a8s' -w 1.0 -t 1000 -r 1000 -n 4000 -p test
#time python ~/desktop/prj/bsff/bsff_scripts/mc_motif_assembly.py -d 'a8s' -w 1.0 -t 10 -r 10 -n 10 -p test


#PARAMETERS
#set frequency bias,
#number of simulated annealing trajectories,
#number of trials at each temperature (ie mc moves),
#number of 'extra' motifs to output as pdb (always outputs final motif from each trajectory,
####this will output N additional motifs with best energy from all accepted mc moves)
#prefix gets attached to output folder with motif pdbs as well as text file with solutions and scores
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-d','--parent_directory', help='compound parent directory', required=True)
parser.add_argument('-t','--n_trajectories', help='number of monte carlo trajectories', required=False)
parser.add_argument('-r','--n_trials', help='number of trials for each monte carlo trajectory', required=False)
parser.add_argument('-n','--n_additional_motifs', help='number of motifs to output in addition to last from each trajectory', required=False)
parser.add_argument('-p','--output_prefix', help='prefix for output folder and motif pdbs', required=False)
#
parser.add_argument('-a','--fa_atr', help='rosetta faatr weight', required=False)
parser.add_argument('-b','--fa_rep', help='rosetta farep weight', required=False)
parser.add_argument('-c','--fa_elec', help='rosetta faelec weight', required=False)
parser.add_argument('-e','--fa_sol', help='rosetta fasol weight', required=False)
parser.add_argument('-f','--hbbbsc', help='rosetta hbbbsc weight', required=False)
parser.add_argument('-g','--hbsc', help='rosetta hbsc weight', required=False)
#
args = vars(parser.parse_args())
if args['parent_directory']:
    parent_directory=args['parent_directory']
else:
    print('must specify compound parent directory')
if args['n_trajectories']:
    n_trajectories=int(args['n_trajectories'])
else:
    n_trajectories=1000
if args['n_trials']:
    n_trials=int(args['n_trials'])
else:
    n_trials=1000
if args['n_additional_motifs']:
    n_additional_motifs=int(args['n_additional_motifs'])
else:
    n_additional_motifs=4000
if args['output_prefix']:
    output_prefix=args['output_prefix']
else:
    output_prefix=args['parent_directory']
#rosetta energy term weights
if args['fa_atr']:
    faatrweight=float(args['fa_atr'])
else:
    faatrweight=1.0
if args['fa_rep']:
    farepweight=float(args['fa_rep'])
else:
    farepweight=0.55
if args['fa_elec']:
    faelecweight=float(args['fa_elec'])
else:
    faelecweight=1.0
if args['fa_sol']:
    fasolweight=float(args['fa_sol'])
else:
    fasolweight=1.0
if args['hbbbsc']:
    hbbbscweight=float(args['hbbbsc'])
else:
    hbbbscweight=1.0
if args['hbsc']:
    hbscweight=float(args['hbsc'])
else:
    hbscweight=1.0


#imports
import sys
import json
import os
import random
import math
from pyrosetta import *
init('-ignore_unrecognized_res -load_PDB_components False')
#load scorefunction
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, faatrweight)
sf.set_weight(fa_rep, farepweight)
sf.set_weight(fa_sol, fasolweight)
sf.set_weight(hbond_sc, hbscweight)
sf.set_weight(fa_elec, faelecweight)
sf.set_weight(hbond_bb_sc,hbbbscweight)
# sf.set_weight(lk_ball_iso,-0.38)

#load fuzzball and frequency json
fuzzball=os.path.join(parent_directory,'clean_fuzzball_for_assembly.pdb')
freqjson=os.path.join(parent_directory,'fuzzball_for_assembly_frequencies.json')
paramspath=os.path.join(parent_directory,'Inputs','Rosetta_Inputs')
params_name=[i for i in os.listdir(paramspath) if i[-6:]=='params']
params=[os.path.join(paramspath,params_name[0])]

'''
#for development:
fuzzball='fuzzball_for_assembly.pdb'
freqjson='fuzzball_for_assembly_frequencies.json'
params=['Inputs/Rosetta_Inputs/dex.params']
'''
#
fuzzball_pose=Pose()
generate_nonstandard_residue_set(fuzzball_pose,params)
pose_from_file(fuzzball_pose, fuzzball)
ligand_resnum=fuzzball_pose.total_residue()
fuzzball_size=fuzzball_pose.size()
#
frequency_dict=json.load(open(freqjson,'r'))

#function to initiate a motif with random residues
lig_pose=fuzzball_pose.residue(ligand_resnum).clone()
def init_random_motif():
    motif_pose=Pose()
    res_in_motif=[]
    while len(res_in_motif)<3:
        random_resnum=random.randint(1,ligand_resnum-1)
        if random_resnum not in res_in_motif:
            res_in_motif.append(random_resnum)
            random_res_pose=fuzzball_pose.residue(random_resnum).clone()
            current_jump = motif_pose.size()
            motif_pose.append_residue_by_jump(random_res_pose, current_jump)
    current_jump = motif_pose.size()
    motif_pose.append_residue_by_jump(lig_pose, current_jump)
    #get the score of this initial motif
    raw_score=sf(motif_pose)
    sumfreqs=0
    for i in res_in_motif:
        frequency=frequency_dict[str(i)][0]
        ref=frequency_dict[str(i)][1]
        stat=stat_weight*(math.log(frequency/ref))
        sumfreqs+=stat
    weighted_score=raw_score+sumfreqs
    return motif_pose,res_in_motif,weighted_score

# motif_pose,res_in_motif,weighted_score=init_random_motif()

#now simulated annealing
temperatures = [50, 10, 5, 1.5, 0.9, 0.6, 0.3]
solutions=[]
all_accepted=[]
e_over_trajectories={}
accepted_seqs=[]
for z in range(n_trajectories):
    motif_pose,res_in_motif,weighted_score=init_random_motif()
    e_over_traj=[]
    e_over_traj.append(weighted_score)
    for kT in temperatures:
        for trial in range(n_trials):
            trial_pose = motif_pose.clone()
            res_in_trial=res_in_motif.copy()
            initial_score=weighted_score
            random_position=random.choice([1,2,3])
            res_in_trial.pop(random_position-1)
            random_resnum=random.randint(1,ligand_resnum-1)
            res_in_trial.insert(random_position-1,random_resnum)
            random_res_from_fuzzball = fuzzball_pose.residue(random_resnum).clone()
            trial_pose.replace_residue(random_position, random_res_from_fuzzball, False)
            raw_score=sf(trial_pose)
            sumfreqs=0
            for i in res_in_trial:
                frequency=frequency_dict[str(i)][0]
                ref=frequency_dict[str(i)][1]
                stat=stat_weight*(math.log(frequency/ref))
                sumfreqs+=stat
            trial_score=raw_score+sumfreqs
            if trial_score<initial_score:
                motif_pose=trial_pose.clone()
                weighted_score=trial_score
                res_in_motif=res_in_trial.copy()
                all_accepted.append((res_in_motif[0],res_in_motif[1],res_in_motif[2],weighted_score))
                e_over_traj.append(weighted_score)
                accepted_seqs.append((motif_pose.sequence()[:3],weighted_score))
            else:
                deltaE=trial_score-initial_score
                P=1-(math.exp(-deltaE/kT))
                randfloat=random.random()
                if randfloat>=P:
                    motif_pose=trial_pose.clone()
                    weighted_score=trial_score
                    res_in_motif=res_in_trial.copy()
                    all_accepted.append((res_in_motif[0],res_in_motif[1],res_in_motif[2],weighted_score))
                    e_over_traj.append(weighted_score)
                    accepted_seqs.append((motif_pose.sequence()[:3],weighted_score))
                else:
                    continue
    # solutions.append((res_in_motif[0],res_in_motif[1],res_in_motif[2]))
    e_over_trajectories[str(z)]=e_over_traj
    # print('trajectory '+str(z+1)+': '+str(weighted_score))

#append 'n additional motifs' to solutions list
saccepted=sorted(all_accepted,key=lambda val: val[3])
solnsonlyres=[]
solnsonlyE=[]
for i in saccepted:
    solnsonlyres.append((i[0],i[1],i[2]))
    solnsonlyE.append(i[3])

for i in solnsonlyres[:n_additional_motifs]:
    if i not in solutions:
        solutions.append(i)

#output solutions
ofile=open(parent_directory+'/'+output_prefix+'_motif_solutions.txt','w')
jsonsolns={}
jsonsolns['Solutions']=solutions
json.dump(jsonsolns,ofile)
ofile.close()

#output solution motif pdb files so they can be used for fastmatch
motif_output_directory=parent_directory+'/'+output_prefix+'_motif_pdbs'
os.makedirs(motif_output_directory,exist_ok=True)
motif_paths=[]
def output_motif(residues):
    motif_res=residues[:3]
    # print(motif_res)
    motif=Pose()
    for res in motif_res:
        residue_pose = fuzzball_pose.residue(res).clone()
        current_jump = motif.size()
        motif.append_residue_by_jump(residue_pose, current_jump)
    current_jump = motif.size()
    motif.append_residue_by_jump(lig_pose, current_jump)
    # print(motif.sequence())
    # sf.show(motif)
    motif_name=(str(motif_res[0])+'_'+str(motif_res[1])+'_'+str(motif_res[2]))
    if os.path.exists(os.path.join(motif_output_directory,motif_name+'.pdb'))==False:
        if len(motif_paths)<n_additional_motifs:
            motif.dump_pdb(os.path.join(motif_output_directory,motif_name+'.pdb'))
            motif_paths.append(os.path.join(motif_output_directory,motif_name+'.pdb'))

#have to clean motif pdbs so its clear motif res arent connected in chain
for motif in solutions:
    output_motif(motif)

clean_motifs_path=os.path.join(motif_output_directory,'clean_motifs')
os.makedirs(clean_motifs_path,exist_ok=True)
for filename in motif_paths:
    f=open(filename,'r')
    lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[0:3]=='TER' or line[0:6]=='HETNAM' or line[0:6]=='HETATM']
    f.close()
    ofile=open(clean_motifs_path+'/'+filename.split('/')[-1],'w')
    cleanlines=[]
    for line in lines:
        if line[0:4]=='ATOM':
            cleanlines.append(line)
        if line[0:6]=='HETATM':
            if line[17:20]=='HOH':
                pass
            else:
                newline='ATOM  '+line[6:]
                cleanlines.append(line)
        if line[:3]=='TER':
            cleanlines.append(line)
    for i,line in enumerate(cleanlines):
        try:
            resnum=int(line[22:29].strip())
            nextresnum=int(cleanlines[i+1][22:29].strip())
            if resnum==nextresnum:
                ofile.write(line)
            else:
                ofile.write(line)
                ofile.write('TER\n')
        except:
            ofile.write(line)
    ofile.close()


for i in os.listdir(motif_output_directory):
    if i[-3:]=='pdb':
        os.system('rm '+os.path.join(motif_output_directory,i))









import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
# from collections import defaultdict

pdfname=parent_directory+'_'+output_prefix+'_mc_stats.pdf'
pdf = PdfPages(pdfname)


# fig = plt.figure()
# ax = fig.add_subplot()
# ax.hist(solnsonlyE,color='blue',alpha=0.5,bins=40)
# ax.axvline(np.mean(solnsonlyE),c='g',linestyle='dashed')
# ax.set_xlabel('Energy (biased REU)')
# ax.set_ylabel('Count')
# ax.set_title('Energy Distribution; All Accepted Moves')
# ax.set_xlim(-25,25)
# ax.text(0.9,0.9,'mu = '+str(np.mean(solnsonlyE))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='black', fontsize=8)
# pdf.savefig()
# plt.clf()
# plt.close()

top_10k_index=10000
top_30k_index=30000
top_50k_index=50000


fig = plt.figure()
ax = fig.add_subplot()
thexst=[i for i in range(1,len(solnsonlyE)+1)]
ax.scatter(thexst,solnsonlyE)
ax.axvline(top_10k_index,c='g',linestyle='dashed')
ax.axvline(top_30k_index,c='g',linestyle='dashed')
ax.axvline(top_50k_index,c='g',linestyle='dashed')
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
ax.set_xlabel('Motif Index')
ax.set_ylabel('Energy (biased REU)')
ax.set_ylim(-25,25)
ax.set_title('Ordered Energy Distribution; All Accepted Moves')
ax.text(0.9,0.9,'mu = '+str(np.mean(solnsonlyE))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='black', fontsize=8)
pdf.savefig()
plt.clf()
plt.close()


import random
randindices=[]
for i in range(100):
    randindices.append(random.randint(0,len(list(e_over_trajectories.keys()))-1))
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xlabel('Step Number')
ax.set_ylabel('Energy (biased REU)')
ax.set_title('Energy over Trajectories')
ax.set_ylim(-25,25)
for rindex in randindices:
    rkey=list(e_over_trajectories.keys())[rindex]
    eotl=e_over_trajectories[rkey]
    thexs=[i for i in range(1,len(eotl)+1)]
    ax.plot(thexs,eotl,linewidth=0.5)
pdf.savefig()
plt.clf()
plt.close()


import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']

#convert 3 letter aa code sequence to 1 letter code sequence
def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3]
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq


accepted_seqs=sorted(accepted_seqs,key=lambda val: val[1])
select_seqs=accepted_seqs[:50000]
#plot amino acid distribution across all positions
aadict={}
sanitylist=[]
for resname in amino_acids:
    aa=threeto1(resname)
    count=0
    for s,tscore in accepted_seqs:
        for i in s:
            if i==aa:
                count+=1
    aadict[resname]=count/(len(accepted_seqs)*3)
    sanitylist.append(count/(len(accepted_seqs)*3))
# print(str(sum(sanitylist)))
scx,scy=list(aadict.keys()),aadict.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
fig = plt.figure()
ax = fig.add_subplot()
ax.bar(scx, scy)
ax.set_xticklabels(scx,rotation='vertical')
ax.set_xlabel('Residue')
ax.set_ylabel('Frequency')
ax.set_title('Amino Acid Frequencies for Motifs; All Accepted Moves')
pdf.savefig()
plt.clf()
plt.close()

#plot amino acid distribution across all positions
aadict={}
sanitylist=[]
for resname in amino_acids:
    aa=threeto1(resname)
    count=0
    for s,score in select_seqs:
        for i in s:
            if i==aa:
                count+=1
    aadict[resname]=count/(len(select_seqs)*3)
    sanitylist.append(count/(len(select_seqs)*3))
# print(str(sum(sanitylist)))
scx,scy=list(aadict.keys()),aadict.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
fig = plt.figure()
ax = fig.add_subplot()
ax.bar(scx, scy)
ax.set_xticklabels(scx,rotation='vertical')
ax.set_xlabel('Residue')
ax.set_ylabel('Frequency')
ax.set_title('Amino Acid Frequencies for Motifs; Top 50k Energy')
pdf.savefig()
plt.clf()
plt.close()


#lists to hold keys for residues of given contact chemistry
aliphatic=[]
aromatic=[]
polar=[]
charged_acidic=[]
charged_basic=[]
glycines=[]
prolines=[]
methionines=[]
for s,tscore in accepted_seqs:
    for contact_resname in s:
        if contact_resname=='A' or contact_resname=='I' or contact_resname=='L' or contact_resname=='V':
           aliphatic.append(contact_resname)
        elif contact_resname=='F' or contact_resname=='Y' or contact_resname=='W':
            aromatic.append(contact_resname)
        elif contact_resname=='S' or contact_resname=='T' or contact_resname=='N' or contact_resname=='Q' or contact_resname=='C':
             polar.append(contact_resname)
        elif contact_resname=='E' or contact_resname=='D':
            charged_acidic.append(contact_resname)
        elif contact_resname=='R' or contact_resname=='K' or contact_resname=='H':
            charged_basic.append(contact_resname)
        elif contact_resname=='G':
            glycines.append(contact_resname)
        elif contact_resname=='P':
            prolines.append(contact_resname)
        elif contact_resname=='M':
            methionines.append(contact_resname)
#bar plot contact chemistry frequencies
total_res=float(len(accepted_seqs)*3)
contact_chem_categories=['aliph','arom','polar','q-',
                         'q+','gly','met']
contact_chem_freqs=[]
naliphatic=len(aliphatic)/total_res;contact_chem_freqs.append(naliphatic)
naromatic=len(aromatic)/total_res;contact_chem_freqs.append(naromatic)
npolar=len(polar)/total_res;contact_chem_freqs.append(npolar)
ncharged_acidic=len(charged_acidic)/total_res;contact_chem_freqs.append(ncharged_acidic)
ncharged_basic=len(charged_basic)/total_res;contact_chem_freqs.append(ncharged_basic)
nglycines=len(glycines)/total_res;contact_chem_freqs.append(nglycines)
nmethionines=len(methionines)/total_res;contact_chem_freqs.append(nmethionines)
sortlist=[]
for i,e in enumerate(contact_chem_freqs):
    sortlist.append((contact_chem_categories[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
contact_chem_categories=[];contact_chem_freqs=[]
for a,b in sortlist:
    contact_chem_categories.append(a);contact_chem_freqs.append(b)
fig = plt.figure()
ax = fig.add_subplot()
ax.bar(contact_chem_categories, contact_chem_freqs)
# plt.xticks(rotation='vertical')
ax.set_xlabel('Contact Chemistry')
ax.set_ylabel('Frequency')
ax.set_title('Contact Chemistry Frequencies Across all Motif Positions; All Accepted Moves')
pdf.savefig()
plt.clf()


all_resnum=[]
for s in all_accepted:
    all_resnum.append(s[0])
    all_resnum.append(s[1])
    all_resnum.append(s[2])
#plot res number distribution
n_diff_res=len(set(all_resnum))
resnum_dict={}
for x in set(all_resnum):
    count=0
    for i in all_resnum:
        if i==x:
            count+=1
    resnum_dict[str(x)]=count
scx,scy=list(resnum_dict.keys()),resnum_dict.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
fig = plt.figure()
ax = fig.add_subplot()
ax.barh(scx, scy)
# plt.xticks(rotation='vertical')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False)
ax.set_ylabel('Residue Number')
ax.set_xlabel('Frequency')
ax.set_title('Residue Number Distribution in Motifs; '+str(n_diff_res)+' unique residues')
pdf.savefig()
plt.clf()
plt.close()

all_resnum=[]
for s in saccepted[:50000]:
    all_resnum.append(s[0])
    all_resnum.append(s[1])
    all_resnum.append(s[2])
#plot res number distribution
n_diff_res=len(set(all_resnum))
resnum_dict={}
for x in set(all_resnum):
    count=0
    for i in all_resnum:
        if i==x:
            count+=1
    resnum_dict[str(x)]=count
scx,scy=list(resnum_dict.keys()),resnum_dict.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
fig = plt.figure()
ax = fig.add_subplot()
ax.barh(scx, scy)
# plt.xticks(rotation='vertical')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False)
ax.set_ylabel('Residue Number')
ax.set_xlabel('Frequency')
ax.set_title('Residue Number Distribution in Motifs, Top 50k; '+str(n_diff_res)+' unique residues')
pdf.savefig()
plt.clf()
plt.close()


pdf.close()
