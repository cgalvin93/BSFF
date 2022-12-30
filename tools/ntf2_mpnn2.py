'''
running mpnn on all of xjs filtered ntf2 scaffolds

/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds
'''
import os
l=[i for i in os.listdir()]
print(len(l))
#1709
for i in l:
    os.system('gunzip '+i)
#
os.makedirs('pdbs',exist_ok=True)
l=[i for i in os.listdir() if i[-3:]=='pdb']
for i in l:
    os.system('cp '+i+' pdbs/'+i)

'''
MPNN

/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs
'''
#making the folders with the pdbs
#i will likewise submit ten per job
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs'
directory_prefix='mpnn'
shf_prefix='ntf2design'
#
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
# import os
#
#
# n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs'
# directory_prefix='mpnn'
# shf_prefix='ntf2design'
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
    '--num_seq_per_target 1000',
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
qsub -cwd -l mem_free=40G -o mpnncout -e mpnncout ntf2ogmpnn_1_.sh

'''
import os
os.makedirs('mpnncout',exist_ok=True)
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)





















#okay, now lets see how long these things take
'''
consolidating mpnn results
'''
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs'
directory_prefix='mpnn'
n_top_seqs_for_cf=5
#
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
badc=0
for id in input_directories:
    try:
        resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
        resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
        for y in resl:
            allresl.append(os.path.join(resultfspath,y))
    except:
        badc+=1
        print(badc)


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
selectfastasoutputdir='cf'
os.makedirs(selectfastasoutputdir,exist_ok=True)
for key in add.keys():
    selectseqs=[i[2] for i in add[key][:n_top_seqs_for_cf]]
    c=1
    for seq in set(selectseqs):
        ofn='_'.join([key,'design',str(c)])
        ofnn=os.path.join(selectfastasoutputdir,ofn+'.fasta')
        of=open(ofnn,'w')
        of.write('>'+ofn+'\n')
        of.write(seq)
        of.close()
        c+=1






os.chdir('cf')
'''
split fasta directory into subdirectories
'''
import os
nfastasperjob=5
directory_prefix='ntf2_designs'
#allfastasdir=cf
allfastasdir=os.getcwd()
####


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
# import os
#
import os
directory_prefix='ntf2_designs'
#allfastasdir=cf
allfastasdir=os.getcwd()
shf_prefix='cfntf2'
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

'''
mkdir cluster_output
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cf_ntf2_r1d_1_.sh
'''
import os
os.makedirs('cluster_out',exist_ok=True)
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('cf_ntf2_r1d_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''


beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

'''






















'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS


import os
#delete empty matches
directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf'
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/'
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
c=0
print(len(prediction_pdb_paths))
#30k
for i in prediction_pdb_paths:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    f.close()
    if len(lines)<40:
        os.system('rm '+i)
    c+=1
    print(str(c))

directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf'
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/'
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))

print(len(prediction_pdb_paths))
'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf'
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/'

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
'''
In [12]: ogds[0]
Out[12]: '134395.pdb'
In [13]: prediction_pdb_paths[0]
Out[13]: '/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/ntf2_designs_2/ntf2_designs_2_output/229026_design_1_unrelaxed_model_2_rank_1.pdb'
'''
print('\n\n\n\n')
print(len(prediction_pdb_paths))
print('\n\n\n\n')
'''
30111
'''
co=0
for p in prediction_pdb_paths[:50]:
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
    # currstrc_name='_'.join(p.split('/')[-1].split('_')[0])###########################################
    currstrc_name=p.split('/')[-1].split('_')[0]###########################################
    csnogds=currstrc_name+'.pdb'
    if csnogds in ogds:
        try:
            p11=pose_from_pdb(p)
            p22=pose_from_pdb(os.path.join(nonaf2desdir,csnogds))
            carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
        except:
            continue
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
    co+=1
    print(str(co))
#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
print('json output')
##################################################################
##################################################################
##################################################################
'''
echo "
import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')
from collections import defaultdict

directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf'
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/'

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

print('\n\n\n\n')
print(len(prediction_pdb_paths))
print('\n\n\n\n')

co=0
for p in prediction_pdb_paths:
    # print(str(p))
    try:
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
        # currstrc_name='_'.join(p.split('/')[-1].split('_')[0])###########################################
        currstrc_name=p.split('/')[-1].split('_')[0]###########################################
        csnogds=currstrc_name+'.pdb'
        if csnogds in ogds:
            try:
                p11=pose_from_pdb(p)
                p22=pose_from_pdb(os.path.join(nonaf2desdir,csnogds))
                carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
            except:
                continue
        ##########
        if currstrc_name not in seen_names:
             seen_names.append(currstrc_name)
             data_allstrc[currstrc_name]=defaultdict(list)
        else:
            pass
        #########
        currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
        currstrc_design=str(p.split('/')[-1].split('_')[-6])
        ##################
        data_allstrc[currstrc_name][currstrc_design].append([currstrc_plddt,currstrc_model,carmsd])
        co+=1
        print(str(co))
    except:
        pass
#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
print('json output')
##################################################################
##################################################################
##################################################################
">getafdata.py

echo "
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 getafdata.py
">getafdata.sh

qsub -cwd -l mem_free=4G getafdata.sh

'''
import json
with open('af2_data.json','r') as f:
    data_allstrc=json.load(f)
'''
In [5]: len(data_allstrc.keys())
Out[5]: 1342

data_allstrc[currstrc_name][currstrc_design].append([currstrc_plddt,currstrc_model,carmsd])

'''
import numpy as np
#get the highest confidence structure from each input
plddt_threshold=75.0
carmsd_threshold=3.0
aplddts=[]
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        for model in data_allstrc[key][k2]:
            cplddt=model[0]
            modeln=model[1]
            carmsd=model[2]
            if cplddt>=plddt_threshold:
                if carmsd<=carmsd_threshold:
                    plddts.append((cplddt,key,k2,modeln,carmsd))
    l=sorted(plddts, reverse=True, key=lambda first: first[0])
    if len(l)>0:
        aplddts.append(l[0])
print(len(aplddts))
'''
1268
'''
#plot all vals
p=[]
c=[]
for key in data_allstrc.keys():
    for k2 in data_allstrc[key]:
        for model in data_allstrc[key][k2]:
            cplddt=model[0]
            carmsd=model[2]
            p.append(cplddt)
            c.append(carmsd)

import matplotlib.pyplot as plt
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
'''
#move filtered to own directory
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
import os
#
directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf'
#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
'''
229026_design_1_unrelaxed_model_2_rank_1.pdb
(86.66570175438598, '229026', '1', '2', 2.0971522331237793)
'''
os.makedirs('af2filtered',exist_ok=True)
for i in aplddts:
    n='_'.join([str(i[1]),'design',str(i[2]),'unrelaxed_model',str(i[3])])
    for j in prediction_pdb_paths:
        if n in j:
            newp=os.path.join('af2filtered',j.split('/')[-1])
            os.system('cp '+j+' '+newp)
os.system('mv *.pdf af2filtered')
'''
In [3]: len(os.listdir())
Out[3]: 1269
'''


'''

next :
    RELAX -> MPNN -> CF
'''



'''
RELAXING AF2 MODELS OF BEST ROUND1 MPNN DESIGNS
in directory with my best strc...
'''
import os
os.mkdir('relaxed')
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('rel.sh','w') ####################enter name
ofile.write('#!/bin/bash')
ofile.write('\n')
ofile.write('tasks=(0\n')
for strc in pdbnames[:-1]:
    ofile.write('       '+strc+'\n')
ofile.write('       '+pdbnames[-1]+')')
ofile.write('\n~/main/source/bin/relax.default.linuxgccrelease -ignore_unrecognized_res -s ${tasks[$SGE_TASK_ID]} -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -in:file:fullatom -relax:fast -out:prefix relaxed_ -out:path:all relaxed')
ofile.write('\nqstat -j "$JOB_ID"')
ofile.close()
print(str(len(pdbnames)))
'''
1268
mkdir cluster_output
qsub -cwd -t 1-1268 -l mem_free=2G -o cluster_output -e cluster_output rel.sh

'''

'''
MPNN

/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed
'''
#making the folders with the pdbs
#i will likewise submit ten per job
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed'
directory_prefix='mpnn'
shf_prefix='ntf2des'
#
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
# import os
#
#
# n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs'
# directory_prefix='mpnn'
# shf_prefix='ntf2design'
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
    '--num_seq_per_target 1000',
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
qsub -cwd -l mem_free=40G -o mpnncout -e mpnncout ntf2ogmpnn_1_.sh

'''
import os
os.makedirs('mpnncout',exist_ok=True)
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)





















#okay, now lets see how long these things take
'''
consolidating mpnn results
'''
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed'
directory_prefix='mpnn'
n_top_seqs_for_cf=5
#
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
badc=0
for id in input_directories:
    try:
        resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
        resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
        for y in resl:
            allresl.append(os.path.join(resultfspath,y))
    except:
        badc+=1
        print(badc)
print(badc)
'''

98

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

#create the directory of fasta files for designs to submit to colabfold
selectfastasoutputdir='cf'
os.makedirs(selectfastasoutputdir,exist_ok=True)
for key in add.keys():
    selectseqs=[i[2] for i in add[key][:n_top_seqs_for_cf]]
    c=1
    for seq in set(selectseqs):
        ofn='_'.join([key,'design',str(c)])
        ofnn=os.path.join(selectfastasoutputdir,ofn+'.fasta')
        of=open(ofnn,'w')
        of.write('>'+ofn+'\n')
        of.write(seq)
        of.close()
        c+=1






os.chdir('cf')
'''
split fasta directory into subdirectories
'''
import os
nfastasperjob=5
directory_prefix='ntf2_designs'
#allfastasdir=cf
allfastasdir=os.getcwd()
####


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
# import os
#
import os
directory_prefix='ntf2_designs'
#allfastasdir=cf
allfastasdir=os.getcwd()
shf_prefix='cfntf2'
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
cd cf
'''
import os
os.makedirs('cluster_out',exist_ok=True)
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('cf_ntf2_r1d_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o cluster_out -e cluster_out '+j
    os.system(cmd)

'''
/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin

'''






















'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS


'''
directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf'
#in this case i moved all the relaxed structures to mpnn input folders, but i need
#to copy them all to one directory to be the nonaf2desdir

# import os
# inputdirs=[i for i in os.listdir() if os.path.isdir(i)==True and i[:4]=='mpnn']
#
# newdir='allrelaxedinputs'
# os.makedirs(newdir,exist_ok=True)
#
# for id in inputdirs:
#     pdbs=[os.path.join(id,i) for i in os.listdir(id) if i[-3:]=='pdb']
#     for p in pdbs:
#         os.system('cp '+p+' '+newdir+'/'+p.split('/')[-1])
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/allrelaxedinputs'


'''
echo "
import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')
from collections import defaultdict

directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf'
nonaf2desdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/allrelaxedinputs'

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

print('\n\n\n\n')
print(len(prediction_pdb_paths))
print('\n\n\n\n')

co=0
for p in prediction_pdb_paths:
    # print(str(p))
    try:
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
        # currstrc_name='_'.join(p.split('/')[-1].split('_')[1])###########################################
        currstrc_name=p.split('/')[-1].split('_')[1]###########################################
        # csnogds=currstrc_name+'.pdb'
        for ogdn in ogds:
            if currstrc_name in ogdn:
                try:
                    p11=pose_from_pdb(p)
                    p22=pose_from_pdb(os.path.join(nonaf2desdir,ogdn))
                    carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
                except:
                    continue
        ##########
        if currstrc_name not in seen_names:
             seen_names.append(currstrc_name)
             data_allstrc[currstrc_name]=defaultdict(list)
        else:
            pass
        #########
        currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
        currstrc_design=str(p.split('/')[-1].split('_')[-6])
        ##################
        data_allstrc[currstrc_name][currstrc_design].append([currstrc_plddt,currstrc_model,carmsd])
        co+=1
        print(str(co))
    except:
        pass
#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
print('json output')
##################################################################
##################################################################
##################################################################
">getafdata.py

echo "
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 getafdata.py
">getafdata.sh

qsub -cwd -l mem_free=4G getafdata.sh

'''
import json
with open('af2_data.json','r') as f:
    data_allstrc=json.load(f)
'''
In [2]: len(data_allstrc.keys())
Out[2]: 629

data_allstrc[currstrc_name][currstrc_design].append([currstrc_plddt,currstrc_model,carmsd])

'''
import numpy as np
#get the highest confidence structure from each input
# plddt_threshold=85.0
# carmsd_threshold=1.5
# aplddts=[]
# for key in data_allstrc.keys():
#     plddts=[]
#     carmsds=[]
#     for k2 in data_allstrc[key]:
#         for model in data_allstrc[key][k2]:
#             cplddt=model[0]
#             modeln=model[1]
#             carmsd=model[2]
#             if cplddt>=plddt_threshold:
#                 if carmsd<=carmsd_threshold:
#                     plddts.append((cplddt,key,k2,modeln,carmsd))
#     l=sorted(plddts, reverse=True, key=lambda first: first[0])
#     if len(l)>0:
#         aplddts.append(l[0])
# print(len(aplddts))
'''
558
'''

#get all passin g threshold
plddt_threshold=85.0
carmsd_threshold=1.5
aplddts=[]
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        for model in data_allstrc[key][k2]:
            cplddt=model[0]
            modeln=model[1]
            carmsd=model[2]
            if cplddt>=plddt_threshold:
                if carmsd<=carmsd_threshold:
                    plddts.append((cplddt,key,k2,modeln,carmsd))
    l=sorted(plddts, reverse=True, key=lambda first: first[0])
    if len(l)>0:
        for ps in l:
            aplddts.append(ps)
print(len(aplddts))
'''
6387
'''
#plot all vals
p=[]
c=[]
for key in data_allstrc.keys():
    for k2 in data_allstrc[key]:
        for model in data_allstrc[key][k2]:
            cplddt=model[0]
            carmsd=model[2]
            p.append(cplddt)
            c.append(carmsd)

import matplotlib.pyplot as plt
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
#move filtered to own directory
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
import os
#
directory_prefix='ntf2_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf'
#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
'''
relaxed_191605_design_3_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_4_rank_1.pdb
(86.66570175438598, '229026', '1', '2', 2.0971522331237793)
'''
os.makedirs('af2filtered',exist_ok=True)
for i in aplddts:
    desid=str(i[1])
    designn=str(i[2])
    modeln=str(i[3])
    found=0
    # n='_'.join(['relaxed',str(i[1]),'design',str(i[2]),'unrelaxed_model',str(i[3])])
    for j in prediction_pdb_paths:
        if j.split('/')[-1].split('_')[1]==desid:
            if j.split('/')[-1].split('_')[-6]==designn:
                if j.split('/')[-1].split('_')[-3]==modeln:
                    # print(j)
                    # print(i)
        # if n in j:
                    newp=os.path.join('af2filtered',j.split('/')[-1])
                    # print(newp)
                    os.system('cp '+j+' '+newp)
                    found+=1
                    break
    if found==0:
        print(i)
os.system('mv *.pdf af2filtered')
'''
import os
len(os.listdir())
In [1]: import os
   ...: len(os.listdir())
Out[1]: 2048

/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf/af2filtered


NEXT IS TO GET ROSETTA METRICS
ACTUALLY MAYBE I SHOULD INCLUDE MULTIPLE DESIGNS FROM EACH BB IF IM ALSO GONNA
FILTER THEM THIS WAY, AND I CAN BE SURE TO FILTER TO ONLY BEST ROSETTA METRICS
FOR EACH IF I WANT TO

'''































'''
RELAXING AF2 MODELS OF BEST ROUND1 MPNN DESIGNS
in directory with my best strc...
'''
import os
os.mkdir('relaxed')
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('rel.sh','w') ####################enter name
ofile.write('#!/bin/bash')
ofile.write('\n')
ofile.write('tasks=(0\n')
for strc in pdbnames[:-1]:
    ofile.write('       '+strc+'\n')
ofile.write('       '+pdbnames[-1]+')')
ofile.write('\n~/main/source/bin/relax.default.linuxgccrelease -ignore_unrecognized_res -s ${tasks[$SGE_TASK_ID]} -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -in:file:fullatom -relax:fast -out:prefix relaxed_ -out:path:all relaxed')
ofile.write('\nqstat -j "$JOB_ID"')
ofile.close()
print(str(len(pdbnames)))
'''
6387
mkdir cluster_output
qsub -cwd -t 1-6387 -l mem_free=2G -o cluster_output -e cluster_output rel.sh







/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf/af2filtered/relaxed
'''
import os
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('r2rm.sh','w') ####################enter name
ofile.write('#!/bin/bash')
ofile.write('\n')
ofile.write('tasks=(0\n')
for strc in pdbnames[:-1]:
    ofile.write('       '+strc+'\n')
ofile.write('       '+pdbnames[-1]+')')
ofile.write('\n~/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol ~/BSFF/tools/scaffold_metrics.xml -parser:view -run:preserve_header -ignore_unrecognized_res -s ${tasks[$SGE_TASK_ID]} -scorefile_format json -out:file:score_only r2scores.json')
ofile.write('\nqstat -j "$JOB_ID"')
ofile.close()
print(str(len(pdbnames)))
'''

mkdir cluster_output
qsub -cwd -t 1-6387 -l mem_free=2G -o cluster_output -e cluster_output r2rm.sh



'''



#analysis of scores from json file
sfname='r2scores.json'
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
print(len(bad))


'''
hphobe_sasa / total_sasa
cav
rescount
bsa
'''
#ccalculate normalized metrics
bc=0
for d in scores:
    try:
        hps=float(d['hphobe_sasa'])
        ts=float(d['total_sasa'])
        rc=float(d['rescount'])
        bnpsa=float(d['bsa'])
        cav=float(d['cav'])
        #
        x=bnpsa/rc
        y=hps/ts
        z=cav/rc
        #
        d['bnpsapr']=x
        d['hsasa_over_tsasa']=y
        d['normcav']=z
    except:
        bc+=1
        print(d)
print(bc)

terms=list(scores[0].keys())
print(terms)




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

plot_dists(terms,scores,'ntf2_rd2_rosetta_metrics.pdf')
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

'''
    d['bnpsapr']=x
    d['hsasa_over_tsasa']=y
    d['normcav']=z
'''
f1=return_filtered(scores,'buns_bb_heavy','<',4.0)
f2=return_filtered(f1,'ps','>',0.6)
f3=return_filtered(f2,'oversat','<',0.0)
f4=return_filtered(f3,'buns_sc_heavy','<',0.0)
f5=return_filtered(f4,'bnpsapr','>',50.)
f6=return_filtered(f5,'hsasa_over_tsasa','<',0.58)
f7=return_filtered(f6,'normcav','<',11.)
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
print(len(f5))
print(len(f6))
print(len(f7))
'''
6206
3935
3097
2931
1934
1934
1816
1816
'''
import os
filtered_strc=[]
for d in f7:
    filtered_strc.append(d['decoy'])

scaffids=[]
for s in filtered_strc:
    scaffid=s.split('_')[2]
    if scaffid not in scaffids:
        scaffids.append(scaffid)
    else:
        pass
len(scaffids)
'''
Out[11]: 428
'''

os.makedirs('filtered_extrascores',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i[:-5]+'.pdb filtered_extrascores/'+i+'.pdb')
    ########################################
os.system('mv *.pdf filtered_extrascores/')
'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf/af2filtered/relaxed/ntf2_rd2_rosetta_metrics.pdf ~/desktop/ntf2_rd2_rosetta_metrics.pdf

/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf/af2filtered/relaxed/filtered_extrascores

zip -r 5tpj_redesigned.zip 5tpj_redesigned

scp cgalvin@log2.wynton.ucsf.edu:5tpj_redesigned.zip ~/desktop/5tpj_redesigned.zip

scp cgalvin@log2.wynton.ucsf.edu:oglucs_5tpj.zip ~/desktop/oglucs_5tpj.zip

'''



#figs for bp223
import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
len(l)
#
ogdir='/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/'
l2=[i for i in os.listdir(ogdir) if i[-3:]=='pdb']
# In [24]: len(l2)
# Out[24]: 1709
#
l3=[i.split('.')[0] for i in l2]
l4=[i.split('_')[2] for i in l]
filteredout=[]
for i in l3:
    if i not in l4:
        filteredout.append(i)
#
len(filteredout)
'''
1281
In [24]: len(l2)
Out[24]: 1709
In [25]: 1709-1281
Out[25]: 428
In [26]: 428/1709
Out[26]: 0.25043885313048564



scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/lucs_scaffolds/5tpj_scaffolds/pdbs/cf/af2filtered/relaxed/cf/plddt_vs_carmsd.pdf ~/desktop/bp223_afconfvsrmsd.pdf

pLDDT
'''
