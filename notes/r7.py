'''
                NTF2 MPNN
starting w the r2 best
scaffold_path='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'

'''
#making the folders with the pdbs
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'
directory_prefix='mpnnjobs'
shf_prefix='mpnn'
##############
import os
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
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(i)==True and i[:len(directory_prefix)]==directory_prefix]
#
submitlist=[]
for id in input_directories:
    outputdirname=id.split('/')[-1]+'_output'
    os.makedirs(os.path.join(id,outputdirname),exist_ok=True)
    submitlist.append((id,os.path.join(id,outputdirname)))
#
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
    of.write('python /wynton/home/kortemme/cgalvin/ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains')
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

os.makedirs('mpnncout',exist_ok=True)
jobss=[i for i in os.listdir() if i[-3:]=='.sh' and shf_prefix in i]
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=4G -o mpnncout -e mpnncout '+j
    os.system(cmd)





















'''
consolidating mpnn results
'''
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'
directory_prefix='mpnnjobs'
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
cb=0
for id in input_directories:
    try:
        resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
        resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
        for y in resl:
            allresl.append(os.path.join(resultfspath,y))
    except:
        cb+=1
        print(str(cb))

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

#create the directory of fasta files for designs to submit to colabfold
selectfastasoutputdir=os.path.join(all_strc_dir,'tocf')
os.makedirs(selectfastasoutputdir,exist_ok=True)
for key in add.keys():
    selectseqs=[i[2] for i in add[key][:10]]
    c=1
    for seq in set(selectseqs):
        ofn='_'.join([key,'design',str(c)])
        ofnn=os.path.join(selectfastasoutputdir,ofn+'.fasta')
        of=open(ofnn,'w')
        of.write('>'+ofn+'\n')
        of.write(seq)
        of.close()
        c+=1
'''
split fasta directory into subdirectories
'''
nfastasperjob=5
directory_prefix='ntf2_rd3_designs'
allfastasdir=selectfastasoutputdir
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

'''
SUBMIT TO COLABFOLD
'''
#
shf_prefix='cf_ntf2_rd3'
#
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
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

os.chdir(allfastasdir)
jobss=[i for i in os.listdir() if i[-3:]=='.sh' and shf_prefix in i]
#
#

os.makedirs('cluster_output',exist_ok=True)
for j in jobss:
    cmd='qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output '+j
    os.system(cmd)

















'''
ANALYZING COLABFOLD RESULTS
/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/tocf

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')
#for the af2 predicxtions
directory_prefix='ntf2_rd3_designs'
#whats the main directory that contains all the af2 job subdirs
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/tocf'
#main directory containing the initial strcs that mpnn was run on
nonaf2desdir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'
#mpnn prefix
mpnn_prefix='mpnnjobs'


og_pdb_paths=[]
#get the og strc pdb paths
input_directories=[os.path.join(nonaf2desdir,i) for i in os.listdir(nonaf2desdir) if os.path.isdir(os.path.join(nonaf2desdir,i))==True and i[:len(mpnn_prefix)]==mpnn_prefix]
for id in input_directories:
    outputdirname=id
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            og_pdb_paths.append(os.path.join(outputdirname,i))

#get the af2 job pdb paths
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))

'''
prediction_pdb_paths[0]
...relaxed_5tpj_403371_design_6_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_design_10_unrelaxed_model_2_rank_1.pdb

og_pdb_path[0]
...relaxed_5tpj_391508_design_2_unrelaxed_model_3_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001.pdb
'''
#
data_allstrc={}
seen_names=[]
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
    currstrc_name='_'.join(p.split('/')[-1].split('_')[:-7])    ###############################################
    csnogds=currstrc_name+'.pdb'
    for nj in og_pdb_paths:
        if csnogds in nj:
            p11=pose_from_pdb(p)
            p22=pose_from_pdb(nj)
            carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
    #########
    currstrc_model=str(p.split('/')[-1].split('_')[-3])
    # currstrc_design_number=str(p.split('/')[-1].split('_')[-6])
    ##########
    name_with_design_id=currstrc_name='_'.join(p.split('/')[-1].split('_')[:-5])
    if name_with_design_id not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    ##################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
# json.dump(data_allstrc,open('a8s3p1np_filt_af2.json','w'))
json.dump(data_allstrc,open('round2relaxedfiltaf2_data.json','w'))
##################################################################
##################################################################
##################################################################

'''
echo ^
getaf2data.py

echo '#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python getaf2data.py
qstat -j "$JOB_ID"'>submit_getaf2data.sh

qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output submit_getaf2data.sh


'''






















'''
FILTERING BASED ON PKLDDT/CARMSD
data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]
'''
directory_prefix='ntf2_rd3_designs'
#whats the main directory that contains all the af2 job subdirs
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/tocf'
#whats the name of the json file witbh all the data
rjson='round2relaxedfiltaf2_data.json'


#
import numpy as np
import os
######
import json
df=open(rjson,'r')
data_allstrc=json.load(df)
df.close()
#######


#id designs with plddt and carmsd below threshold
plddt_threshold=90.0
carmsd_threshold=1.0
#list holding ones that meet threshhold
aplddts=[]
#if average plddt and carmsd across 5 models meets threshold,
#take the highest plddt strc to use for matching (but only if its carmsd also meets threshold)
#
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        plddts.append(cplddt)
        carmsds.append(carmsd)
        # if cplddt>=plddt_threshold:
        #     if carmsd<=carmsd_threshold:
        #         aplddts.append((aplddt,key,k2,acarmsd))
    aplddt=np.mean(plddts)
    mplddt=max(plddts)
    # acarmsd=min(carmsds)
    acarmsd=np.mean(carmsds)
    if aplddt>=plddt_threshold:
        if acarmsd<=carmsd_threshold:
            for k2 in data_allstrc[key]:
                cplddt=data_allstrc[key][k2][0]
                carmsd=data_allstrc[key][k2][2]
                if cplddt==mplddt:
                    if carmsd<=carmsd_threshold:
                        aplddts.append((aplddt,key,k2,acarmsd))
#how many designs did af2 prediction complete on
print(len(list(data_allstrc.keys())))
#how many made it
print(len(aplddts))
#how many unique backbones (input scaffold for mpnn, so excluding last design number)
l=[]
for i in aplddts:
    n=i[1]
    nn=n.split('_')[:-2]        ##############################
    l.append(('_').join(nn))
print(len(list(set(l))))
#how many unique og scaffold ids
l=[]
for i in aplddts:
    n=i[1]
    nn=n.split('_')[2]        ##############################
    l.append(nn)
print(len(list(set(l))))
'''
2939
515
148
65
'''

#get the af2 job pdb paths
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))

#move filtered to own directory to run relax on
toget=[]
for i in aplddts:
    n='_'.join([i[1],'unrelaxed_model',i[2]])
    for p in prediction_pdb_paths:
        if n in p:
            toget.append(p)
os.makedirs('af2filtered',exist_ok=True)
for f in toget:
    newf=os.path.join('af2filtered',f.split('/')[-1])
    os.system('cp '+f+' '+newf)



'''
RELAXING AF2 MODELS
in directory with my best strc...

constrain to start coords? yeah probably, i assume cf2 gonna be pretty close to correct
dont want it to blow up for any weird reasons ya know
'''
import os
os.makedirs('relaxed',exist_ok=True)
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('rel_cf.sh','w') ####################enter name
ofile.write('#!/bin/bash')
ofile.write('\n')
ofile.write('tasks=(0\n')
for strc in pdbnames[:-1]:
    ofile.write('       '+strc+'\n')
ofile.write('       '+pdbnames[-1]+')')
ofile.write('\n~/main/source/bin/relax.default.linuxgccrelease -ignore_unrecognized_res -s ${tasks[$SGE_TASK_ID]} -ex1 -ex2 -extrachi_cutoff 0 -no_optH false -flip_HNQ -relax:constrain_relax_to_start_coords -relax:coord_cst_stdev 0.5 -nstruct 1 -in:file:fullatom -relax:fast -out:prefix relaxed_ -out:path:all relaxed')
ofile.write('\nqstat -j "$JOB_ID"')
ofile.close()
print(str(len(pdbnames)))
'''
mkdir cluster_output
qsub -cwd -t 1-515 -l mem_free=4G -o cluster_output -e cluster_output rel_cf.sh
'''




























#gonna want to get rmsd between all members of same og scaff and plot them all
#also between members of different families and plot them all
#i wanna see all the values
#within fams should be low, between should be high (or at least higher)
#
#IN DIRECTORY WITH RELAXDED STRUCTURES
import os
from collections import defaultdict
import json
relstrc=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(relstrc))
#how many unique input scaffold for mpnn, so excluding last design number
l=[]
for i in relstrc:
    nn=i.split('_')[:-8]        ##############################
    l.append(('_').join(nn))
print(len(list(set(l))))
#how many unique og scaffold ids
l2=[]
scaffold_clusters=defaultdict(list)
for i in relstrc:
    nn=i.split('_')[3]        ##############################
    l2.append(nn)
    scaffold_clusters[nn].append(i)
print(len(list(set(l2))))
'''
515
148
65
'''
#now get rmsds within families
from pyrosetta import *
init('-ignore_unrecognized_res')

#this was to do literally all pairwise rmsds but i think thats overkill
# scaffold_clusters_rmsds={}
# for key in list(scaffold_clusters.keys())[:1]:
#     all_pairwise_carmsds=[]
#     members=[i for i in scaffold_clusters[key]]
#     for ind,ent in enumerate(members[:-1]):
#         p1=ent
#         for ent2 in members[ind+1:]:
#             p2=ent2
#             p11=pose_from_pdb(p1)
#             p22=pose_from_pdb(p2)
#             carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
#             all_pairwise_carmsds.append((p1,p2,carmsd))
#     scaffold_clusters_rmsds[key]=all_pairwise_carmsds
#
# json.dump(scaffold_clusters_rmsds,open('scaffold_clusters_rmsds.json','w'))

scaffold_clusters_rmsds={}
for key in list(scaffold_clusters.keys()):
    all_pairwise_carmsds=[]
    members=[i for i in scaffold_clusters[key]]
    for ind,ent in enumerate(members[:1]):
        p1=ent
        for ent2 in members[ind+1:]:
            p2=ent2
            p11=pose_from_pdb(p1)
            p22=pose_from_pdb(p2)
            carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
            all_pairwise_carmsds.append((p1,p2,carmsd))
    scaffold_clusters_rmsds[key]=all_pairwise_carmsds
json.dump(scaffold_clusters_rmsds,open('scaffold_clusters_rmsds.json','w'))

#get rmsds across og scaffold families
scaffold_rmsds={}
sl=list(scaffold_clusters.keys())
for key in sl:
    all_pairwise_carmsds=[]
    members=[i for i in scaffold_clusters[key]]
    p1=members[0]
    p11=pose_from_pdb(p1)
    for key2 in sl:
        if key2!=key:
            members2=[i for i in scaffold_clusters[key2]]
            p2=members2[0]
            p22=pose_from_pdb(p2)
            carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
            all_pairwise_carmsds.append((p1,p2,carmsd))
    scaffold_rmsds[key]=all_pairwise_carmsds
json.dump(scaffold_rmsds,open('scaffold_rmsds.json','w'))

'''
'363097'


len(scaffold_clusters_rmsds['363097'])
len(scaffold_clusters['363097'])

In [30]: len(scaffold_clusters_rmsds['363097'])
    ...:
Out[30]: 153

In [31]: len(scaffold_clusters['363097'])
    ...:
Out[31]: 18

18 choose 2 = 153, looks good here


OKAY SO LETS SET UP TO RUN ON CLUSTER NODE
        ACTUALLY WAIT DO I REALLY EVEN CARE ENOUGH TO DO LITERALLY ALL PAIRWIE RMSDS?
        MAYBE ITS FINE TO JUST DO PAIRWISE RMSDS OF ONE REPRESENTATIVE STRUCTURE WITH
        ALL OTHER MEMBERS
            this is probably good enough to just give me a decent idea of how different
            the structures are yeah....
AFTER CHANING IT TO BE LIKE THIS
In [35]: len(scaffold_clusters_rmsds['363097'])
    ...:
Out[35]: 17

In [36]: len(scaffold_clusters['363097'])
    ...:
Out[36]: 18

okay so that looks GOOD

echo ^
getrmsd.py

echo '#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python getrmsd.py
qstat -j "$JOB_ID"'>getrmsd.sh

qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output getrmsd.sh



'''
'''
PLOT THE RMSDS OF REPRESENTATIVE STRUCTURES WITHIN/ACROSS FAMILIES
'''
scj='scaffold_clusters_rmsds.json'
sj='scaffold_rmsds.json'
import json


########
df=open(scj,'r')
scaffold_clusters_rmsds=json.load(df)
df.close()
df=open(sj,'r')
scaffold_rmsds=json.load(df)
df.close()
#######


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
##################################################################
outfilename='ntf2round3rmsds.pdf'
##################################################################
pdf = PdfPages(outfilename)
#get list of avg plddt values for each strc to plot
withinrmsds=[]
acrossrmsds=[]
for key in scaffold_clusters_rmsds.keys():
    famvals=[i[2] for i in scaffold_clusters_rmsds[key]]
    acrvals=[i[2] for i in scaffold_rmsds[key]]
    for x in famvals:
        withinrmsds.append(float(x))
    for y in acrvals:
        acrossrmsds.append(float(y))
#plot them
fig = plt.figure()
ax = fig.add_subplot()
#
ax.hist(withinrmsds,color='r',alpha=0.3)#bins=int(len(allscores)/20)
ax.hist(acrossrmsds,color='b',alpha=0.3)#bins=int(len(allscores)/20)
#
ax.axvline(np.mean(withinrmsds),c='r',linestyle='dashed')
ax.axvline(np.mean(acrossrmsds),c='b',linestyle='dashed')
ax.text(0.9,0.85,'mean RMSD within scaffold family = '+str(np.mean(withinrmsds))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
ax.text(0.9,0.9,'mean RMSD across scaffold family = '+str(np.mean(acrossrmsds))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
#
ax.set_title('CA RMSD within/across NTF2 Scaffold "Families"')
ax.set_ylabel('Frequency')
ax.set_xlabel('RMSD')
pdf.savefig()
plt.clf()
#
pdf.close()

'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/af2filtered/relaxed/ntf2round3rmsds.pdf ~/desktop/ntf2round3rmsds.pdf
'''


'''
MAKE POS FILES for matching
'''
import os
from pyrosetta import *
init('-ignore_unrecognized_res')
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



            ROUND 7 WITH DOG


'''

#create a useful directory structure for target
time python ~/desktop/BSFF/new_target.py dog

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

compounds=("dog")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}.sdf
done

move params and pdb to Inputs/Rosetta_Inputs:
compounds=("dog")
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

        DOG
    1-hydroxycyclohexane
    2-adjacent cyclohexane
    3-third cyclohexane
    4-pentane
    5-lactone

compounds=("dog")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
    cd ..
done

'''



'''
NOW setup to split into multiple jobs for alignment on cluster
conda activate pyr37
compounds=("dog")
for I in ${compounds[@]}; do
    cd ${I}
    mkdir cluster_output
    mv *.pdb Inputs/Fragment_Inputs/
    mv *.mol Inputs/Fragment_Inputs/
    time python ~/desktop/BSFF/process_fragment_search.py
    cd ..
done

getting smiles...
135



okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r dog cgalvin@log2.wynton.ucsf.edu:












SUBMIT ALIGNMENT JOB
qsub -cwd -t 1-135 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh


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
qsub -cwd -t 1-47 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh



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
    time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/dog/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb dog $I
done

mkdir residue_contact_fuzzballs
mv *_residue_contacts residue_contact_fuzzballs













CLUSTER NONPOLAR -
compound_dir='dog'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/cluster_nonpolar.py dog /wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb Fragment_1 Fragment_2 Fragment_3 Fragment_4 Fragment_5
qstat -j "$JOB_ID"
'>cnp_$compound_dir.sh
chmod ugo+x cnp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cnp_$compound_dir.sh




















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
lig_resname='dog'
allowed_problematic_motifs=50
#new directory where motifs will be output
# polar_cst_paths='hb4'
polar_cst_paths='hb5'


#for each polar fragment, define whether acceptor is sp2 hybridized or not
#and what are the 3 constrained atoms
# relevant_atoms={'Fragment_1':['yes','O5','C23','C22'],
#                 'Fragment_2':['n','O1','C1','C3'],
#                 'Fragment_3':['n','O2','C6','C8'],
#                 'Fragment_4':['n','O3','C19','C18']}

#to put 2 hb on carbonyl
relevant_atoms={'Fragment_1':['yes','O5','C23','C22'],
                'Fragment_2':['n','O1','C1','C3'],
                'Fragment_3':['n','O2','C6','C8'],
                'Fragment_4':['n','O3','C19','C18'],
                'Fragment_5':['yes','O5','C23','C22']}

#define wheteher or not vcharged residues are allowed to be used with frag
charged_allowed={'Fragment_1':'yes',
                'Fragment_2':'no',
                'Fragment_3':'no',
                'Fragment_4':'no',
                'Fragment_5':'yes'}

#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
q_polar_residues=['LYS','ARG','SER','THR','GLN','ASN','TYR','TRP','GLY'] #no his
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY'] #no his
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






#for each fragment, create a cst block
motifs=[]
nproblematic=0
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
                                                           distance_tolerance=0.2, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=20, torsion_AB_tolerance=180, torsion_B_tolerance=180,
                                                           torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                           distance_constraint_sample_number=1)
        else:
            cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=180,torsionA=180,torsionAB=0,torsionB=0,
                                                           residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                           distance_tolerance=0.2, angle_A_tolerance=10, angle_B_tolerance=15,
                                                           torsion_A_tolerance=180, torsion_AB_tolerance=180, torsion_B_tolerance=180,
                                                           torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                           distance_constraint_sample_number=1)
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
fragment_fuzzballs_dirs=['/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_1_residue_contacts','/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_2_residue_contacts','/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_3_residue_contacts','/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_4_residue_contacts','/wynton/home/kortemme/cgalvin/dog/residue_contact_fuzzballs/dog_Fragment_5_residue_contacts']
ligparams=['/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params']
# polar_cst_paths='/wynton/home/kortemme/cgalvin/dog/idealized_hbond_csts_allpolar1hb'
# polar_cst_paths='/wynton/home/kortemme/cgalvin/dog/idealized_hbond_csts_carbonyl2hb'
motifsoutdir='npmotifpdbs'
############
single_contact_threshold=-1.5
double_contact_score_threshold=-4.0
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
                                                       distance_constraint_sample_number=0)
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
        for i in range(2):          #this is how many unique hybrid cst i want to make per polar cst
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
        for i in range(2): #how many unique per input polar cst
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



































'''
                            MATCHING
'''
###########################################################################
###########################################################################
###########################################################################
###########################################################################
import os
lig_name='dog'
paramspath='/wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params'
allmotifsdir=os.getcwd()
#
# shellfile_suffix='_dog4_1.sh'#
# shellfile_suffix='_dog4_2.sh'#
# shellfile_suffix='_dog5_1.sh'#
# shellfile_suffix='_dog5_2.sh'#
# shellfile_suffix='_dog4.sh'
shellfile_suffix='_dog5.sh'
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
    idx=[j for j in range(0,len(motif_paths[key]),25)]  #how many csts per job
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
96
n scaffolds
515

~/main/source/bin/match.linuxgccrelease -s /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001.pdb -extra_res_fa /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params -match:geometric_constraint_file /wynton/home/kortemme/cgalvin/dog/idealized_hbond_csts_allpolar1hb/1np/hybrid_7.cst -match::scaffold_active_site_residues /wynton/home/kortemme/cgalvin/ntf2_round3_best/relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001.pos -match:output_format PDB -match:match_grouper SameSequenceGrouper -match:consolidate_matches -match:output_matches_per_group 1 -use_input_sc -ex1 -ex2 -extrachi_cutoff 0 -enumerate_ligand_rotamers false -match::lig_name dog -load_PDB_components false

qsub -cwd -t 1-515 -l mem_free=30G -o cluster_out -e cluster_out _1_dog_ideal.sh

scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/dog/test_idealized_hbonds/matchex ~/desktop/matchex
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('_1_a8s_3p1np.sh')
#
for j in jobss:
    cmd='qsub -cwd -t 1-515 -l mem_free=26G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''
WELL THE HYBRID MOTIFS BARELY PRODUCE ANY MATCHES
    0 FOR 5-1
    7 for 4-1
GONNA TRY MATCHING JUST THE 4S AND 5S
    but ill make 100 instead of just 50
'''



























import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
len(l)

'''
100 motifs each
hb4
Out[1]: 1077

hb5
Out[1]: 8

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
