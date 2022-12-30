

'''
                    NTF2 AF2
mkdir ntf2
'''
import os
l=[i for i in os.listdir() if i[:4]=='5tpj' and i[-3:]=='pdb']
for i in l:
    os.system('cp '+i+' ntf2/'+i)

'''
cd ntf2
'''
import os
l=[i for i in os.listdir() if i[:4]=='5tpj' and i[-3:]=='pdb']
# In [2]: len(l)
# Out[2]: 125

'''
make fasta files
'''
from pyrosetta import *
init('-ignore_unrecognized_res')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())[:-1]
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

os.makedirs('ogfastas',exist_ok=True)
l2=[i for i in os.listdir() if i[-6:]=='.fasta']
for i in l2:
    os.system('mv '+i+' ogfastas/'+i)

'''
split fasta directory into subdirectories
'''
import os
nfastasperjob=10
directory_prefix='ntf2og'
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/ogfastas'
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

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate st
time python ~/BSFF/tools/run_colabfold.py input_dir output_dir
qstat -j "$JOB_ID"
'>run_cftest.sh
qsub -cwd -l mem_free=32G -o cluster_output -e cluster_output run_cftest.sh

'''
import os
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/ogfastas'
directory_prefix='ntf2og'
shf_prefix='ntf2ogcf'
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
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
jobss.remove('ntf2ogcf_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=40G -o cluster_output -e cluster_output '+j
    os.system(cmd)








'''
                NTF2 MPNN


'''
#making the folders with the pdbs
#i will likewise submit ten per job
n_strc_per_job=10
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2'
directory_prefix='ntf2ogstrc'
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
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2'
directory_prefix='ntf2ogstrc'
shf_prefix='ntf2ogmpnn'
#
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
    '--sampling_temp "0.25"',
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
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
jobss.remove('ntf2ogmpnn_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=40G -o mpnncout -e mpnncout '+j
    os.system(cmd)


#okay, now lets see how long these things take
'''
mpnn jobs only took about an hour!
af2 still running after like 24 hr,
    seems like its done about 7/10 strc
    BUT I CAN EASILY JUST SUBMIT THESE AS SINGLE SEQ PER JOB BRO...
    or like at least 5
ALSO IN BOTH CASES NOT USING A LOT OF MEMORY
    af2 jobs only ~7gb at most,
    mpnn jobs only 1.5 gb .....



AIGHT BET, SO
    i wanna consolidate my mpnn results,
    analyze them,
    and submit my top sequences to colabfold,
        doing less seq per job (at most 5 and that should be 24hr tops)

consolidating mpnn results
'''
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2'
directory_prefix='ntf2ogstrc'
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
for id in input_directories:
    resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
    resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
    for y in resl:
        allresl.append(os.path.join(resultfspath,y))
'''
analyze the mpnn results and select sequences for af2 modeling

    are the score dists different for different input bbs?
    des score dist vs og score dist vs lowest score dist
    seq recov dist
    seq entropy of designs
take top ten scoring strc for each to submit to colabfold
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
selectfastasoutputdir='/wynton/home/kortemme/cgalvin/ntf2/rd1mpnnresults/tocf'
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
import os
nfastasperjob=2
directory_prefix='ntf2_rd1_designs'
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd1mpnnresults/tocf'
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
import os
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd1mpnnresults/tocf'
directory_prefix='ntf2_rd1_designs'
shf_prefix='cf_ntf2_r1d'
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
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
jobss.remove('cf_ntf2_r1d_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output '+j
    os.system(cmd)

'''
okay so next thing is going to be consolidating af2 results for both
og and these select round1 designs, and comparing plddt distributions,
then maybe selecting top X to relax with rosetta
    and submit relaxed strc for another round of mpnn


beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin


a couple of rpund1 af2s are still running, but the og results are done,
lets take
a look
    outputs are the structures themselves with lddt in the b factor column

/wynton/home/kortemme/cgalvin/ntf2/ogfastas/ntf2og_1/ntf2og_1_output/pdbs...








OG AF2 ANALYSIS
'''
import os
import numpy as np
import json

directory_prefix='ntf2og'
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/ogfastas'
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
for p in prediction_pdb_paths:
    print(str(p))
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
    currstrc_name=p.split('/')[-1].split('_')[2]
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    currstrc_model=str(p.split('/')[-1].split('_')[5])
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('ntf2_og_af2.json','w'))
##################################################################
##################################################################
##################################################################
'''
rjson='ntf2_og_af2.json'
df=open(rjson,'r')
data_allstrc=json.load(df)
df.close()
'''
#plot avg plddt (across 5 models) for each structure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
##################################################################
outfilename='ntf2_og_af2.pdf'
##################################################################
pdf = PdfPages(outfilename)
#get list of avg plddt values for each strc to plot
aplddts=[]
for key in data_allstrc.keys():
    plddts=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        plddts.append(cplddt)
    aplddt=np.mean(plddts)
    aplddts.append(aplddt)
#plot them
fig,ax=plt.subplots()
ax.hist(aplddts)#bins=int(len(allscores)/20)
ax.set_title('Average plddt for each Sequence across 5 Models')
ax.set_ylabel('Frequency')
ax.set_xlabel('Plddt')
pdf.savefig()
plt.clf()
#
pdf.close()
'''
and with perfect timing, the round1 mpnn designs af2 resultds are finished !

scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_og_af2.pdf ~/desktop/ntf2_og_af2.pdf






ROUND 1 MPNN DESIGNS AF2 RESULTS ANALYSIS AND COMPARISON WITH OG!
        w/ the designs the pdb strings are different
        additionally has design number (10 submitted per input)
'''
import os
import numpy as np
import json
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd1mpnnresults/tocf'
directory_prefix='ntf2_rd1_designs'
jsonname='ntf2_mpnnrd1_af2.json'
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
for p in prediction_pdb_paths:
    print(str(p))
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
    ########################################
    cstrc_name=p.split('/')[-1].split('_')[1]
    currstrc_design_number=str(p.split('/')[-1].split('_')[3])
    currstrc_name=('_').join([cstrc_name,currstrc_design_number])
    ########################################
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    ########################################
    currstrc_model=str(p.split('/')[-1].split('_')[6])
    ########################################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open(jsonname,'w'))
##################################################################
##################################################################
##################################################################
'''
import json
rjson='ntf2_og_af2.json'
df=open(rjson,'r')
ogdata=json.load(df)
df.close()

rjson='ntf2_mpnnrd1_af2.json'
df=open(rjson,'r')
data_allstrc=json.load(df)
df.close()
'''
#plot avg plddt (across 5 models) for each structure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
##################################################################
outfilename='ntf2_og_vs_rd1_plddts.pdf'
##################################################################
pdf = PdfPages(outfilename)
#get list of avg plddt values for each strc to plot
aplddts=[]
dd=defaultdict(list)
for key in data_allstrc.keys():
    plddts2=[]
    plddts=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        plddts.append(cplddt)
        plddts2.append((cplddt,k2))
    aplddt=np.mean(plddts)
    aplddts.append(aplddt)
    plddts2=sorted(plddts2, reverse=True,key=lambda first: first[0])
    lowplddt=plddts2[0]
    dd[key.split('_')[0]].append((lowplddt[0],key.split('_')[1],lowplddt[1]))
#get the lowest plddt for each backbone
best_designs={}
for key in dd.keys():
    l=dd[key]
    l=sorted(l, reverse=True,key=lambda first: first[0])
    best_designs[key]=l[0]
aplddts3=[]
for key in best_designs.keys():
    aplddts3.append(best_designs[key][0])
#get other list of avg plddt values for each strc to plot
aplddts2=[]
for key in ogdata.keys():
    plddts=[]
    for k2 in ogdata[key]:
        cplddt=ogdata[key][k2][0]
        plddts.append(cplddt)
    aplddt=np.mean(plddts)
    aplddts2.append(aplddt)
#plot them
fig = plt.figure()
ax = fig.add_subplot()
#
ax.hist(aplddts,color='r',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts2,color='b',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts3,color='g',alpha=0.5)#bins=int(len(allscores)/20)

#
ax.axvline(np.mean(aplddts),c='r',linestyle='dashed')
ax.axvline(np.mean(aplddts2),c='b',linestyle='dashed')
ax.axvline(np.mean(aplddts3),c='g',linestyle='dashed')
ax.text(0.9,0.85,'Rd1 mean = '+str(np.mean(aplddts))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
ax.text(0.9,0.9,'OG mean = '+str(np.mean(aplddts2))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
ax.text(0.9,0.95,'Rd1 bests mean = '+str(np.mean(aplddts3))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
#
ax.set_title('Average plddt for each Sequence across 5 Models')
ax.set_ylabel('Frequency')
ax.set_xlabel('Plddt')
pdf.savefig()
plt.clf()
#
pdf.close()
'''
scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_og_vs_rd1_plddts.pdf ~/desktop/ntf2_og_vs_rd1_plddts.pdf

NOW GETTING BEST DESIGNS TO A NEW DIRECTORY TO THEN SUBMIT 4 RELAX->MPNN->CF
best_designs[lucsmodel]=(plddt,design_number,af2_model)
'''
import os
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd1mpnnresults/tocf'
directory_prefix='ntf2_rd1_designs'
outputdirectory='/wynton/home/kortemme/cgalvin/ntf2/round1_best'
os.makedirs(outputdirectory,exist_ok=True)
#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
#
l=[]
for p in prediction_pdb_paths:
    cstrc_name=p.split('/')[-1].split('_')[1]
    best_d_id=str(best_designs[cstrc_name][1])
    best_model_id=str(best_designs[cstrc_name][2])
    currstrc_design_number=str(p.split('/')[-1].split('_')[3])
    currstrc_model=str(p.split('/')[-1].split('_')[6])
    if currstrc_design_number==best_d_id:
        if currstrc_model==best_model_id:
            cmd='cp '+p+' '+os.path.join(outputdirectory,p.split('/')[-1])
            # os.system('cp '+p+' '+os.path.join(outputdirectory,p.split('/')[-1]))
            l.append(cmd)

len(list(best_designs.keys()))





'''
RELAXING AF2 MODELS OF BEST ROUND1 MPNN DESIGNS
in directory with my best strc...
'''
import os
# os.mkdir('relaxed')
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('rr1b.sh','w') ####################enter name
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
mkdir cluster_output
qsub -cwd -t 1-110 -l mem_free=4G -o cluster_output -e cluster_output rr1b.sh


I SOMEHOW ENDED UP WITH ONLY 110 INSTEAD OF THE 125 BBS I STARTED WITH...
but whatever im running with it for right now, i dont wanna get sidetracked
can investigate later after i get designs n shit

AY DONT FORGET TO COMPARE ROSETTA METRICS OF OG NTF2S AND DESIGNS












NEXT THOUGH IS SUBMITTING ROUND1 BEST TO ANOTHER ROUND OF MPNN
'''
#making the folders with the pdbs
#i will likewise submit ten per job
n_strc_per_job=10
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round1_best/relaxed'
directory_prefix='r1b'
###############################################################
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
import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round1_best/relaxed'
directory_prefix='r1b'
shf_prefix='r1bmpnn2'
#
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
            LOWERING TEMP TO 0.15 AND MAKING 1K SEQS PER STRC
                (its so fucking fast why not dawg...)
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=16G -o mpnncout -e mpnncout '+j
    os.system(cmd)
















import os
import numpy as np
import json
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/tocf'
directory_prefix='r1bd'
jsonname='ntf2_mpnnrd2_af2.json'
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
for p in prediction_pdb_paths:
    print(str(p))
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
    ########################################
    cstrc_name=p.split('/')[-1].split('_')[2]
    currstrc_design_number=str(p.split('/')[-1].split('_')[12])
    currstrc_name=('_').join([cstrc_name,currstrc_design_number])
    ########################################
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    ########################################
    currstrc_model=str(p.split('/')[-1].split('_')[-3])
    ########################################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open(jsonname,'w'))
##################################################################
##################################################################
##################################################################
'''
import json
rjson='ntf2_og_af2.json'
df=open(rjson,'r')
ogdata=json.load(df)
df.close()

rjson='ntf2_mpnnrd1_af2.json'
df=open(rjson,'r')
data_allstrc=json.load(df)
df.close()

rjson='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/ntf2_mpnnrd2_af2.json'
df=open(rjson,'r')
data_allstrc2=json.load(df)
df.close()
'''
#plot avg plddt (across 5 models) for each structure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
##################################################################
outfilename='ntf2_og_rd1_rd2_plddts.pdf'
##################################################################
pdf = PdfPages(outfilename)
#get list of avg plddt values for each strc to plot rd1
aplddts=[]
dd=defaultdict(list)
for key in data_allstrc.keys():
    plddts2=[]
    plddts=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        plddts.append(cplddt)
        plddts2.append((cplddt,k2))
    aplddt=np.mean(plddts)
    aplddts.append(aplddt)
    plddts2=sorted(plddts2, reverse=True,key=lambda first: first[0])
    lowplddt=plddts2[0]
    dd[key.split('_')[0]].append((lowplddt[0],key.split('_')[1],lowplddt[1]))
#get the lowest plddt for each backbone rd1
best_designs={}
for key in dd.keys():
    l=dd[key]
    l=sorted(l, reverse=True,key=lambda first: first[0])
    best_designs[key]=l[0]
aplddts3=[]
for key in best_designs.keys():
    aplddts3.append(best_designs[key][0])
#get other list of avg plddt values for each strc to plot og
aplddts2=[]
for key in ogdata.keys():
    plddts=[]
    for k2 in ogdata[key]:
        cplddt=ogdata[key][k2][0]
        plddts.append(cplddt)
    aplddt=np.mean(plddts)
    aplddts2.append(aplddt)
#get list of avg plddt values for each rd2 strc to plot
aplddts4=[]
dd2=defaultdict(list)
for key in data_allstrc2.keys():
    plddts=[]
    plddts2=[]
    for k2 in data_allstrc2[key]:
        cplddt=data_allstrc2[key][k2][0]
        plddts.append(cplddt)
        plddts2.append((cplddt,k2))
    aplddt=np.mean(plddts)
    aplddts4.append(aplddt)
    plddts2=sorted(plddts2, reverse=True,key=lambda first: first[0])
    # plddts2=sorted(plddts2,key=lambda first: first[0])
    lowplddt=plddts2[0]
    dd2[key.split('_')[0]].append((lowplddt[0],key.split('_')[1],lowplddt[1]))
#get the lowest plddt for each backbone
best_designs2={}
for key in dd2.keys():
    l=dd2[key]
    l=sorted(l, reverse=True,key=lambda first: first[0])
    best_designs2[key]=l[0]
aplddts5=[]
for key in best_designs2.keys():
    aplddts5.append(best_designs2[key][0])

#plot them
fig = plt.figure()
ax = fig.add_subplot()
#
ax.hist(aplddts,color='r',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts2,color='b',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts3,color='g',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts4,color='y',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts5,color='orange',alpha=0.5)#bins=int(len(allscores)/20)

#
ax.axvline(np.mean(aplddts),c='r',linestyle='dashed')
ax.axvline(np.mean(aplddts2),c='b',linestyle='dashed')
ax.axvline(np.mean(aplddts3),c='g',linestyle='dashed')
ax.axvline(np.mean(aplddts4),c='yellow',linestyle='dashed')
ax.axvline(np.mean(aplddts5),c='orange',linestyle='dashed')
ax.text(0.9,0.85,'Rd1 mean = '+str(np.mean(aplddts))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
ax.text(0.9,0.9,'OG mean = '+str(np.mean(aplddts2))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
ax.text(0.9,0.95,'Rd1 bests mean = '+str(np.mean(aplddts3))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
ax.text(0.9,0.8,'rd2 mean = '+str(np.mean(aplddts4))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='yellow', fontsize=8)
ax.text(0.9,0.75,'Rd2 bests mean = '+str(np.mean(aplddts5))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='orange', fontsize=8)
#
ax.set_title('Average plddt for each Sequence across 5 Models')
ax.set_ylabel('Frequency')
ax.set_xlabel('Plddt')
pdf.savefig()
plt.clf()
#
pdf.close()
'''
scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_og_rd1_rd2_plddts.pdf ~/desktop/ntf2_og_rd1_rd2_plddts.pdf

scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_og_rd1_plddts.pdf ~/desktop/ntf2_og_rd1_plddts.pdf
scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_rd1__rd2_plddts.pdf ~/desktop/ntf2_rd1__rd2_plddts.pdf
scp cgalvin@log2.wynton.ucsf.edu:ntf2/ntf2_og__rd2_plddts.pdf ~/desktop/ntf2_og__rd2_plddts.pdf

ntf2_og_rd1_plddts.pdf
ntf2_rd1__rd2_plddts.pdf
ntf2_og__rd2_plddts.pdf


outfilename='ntf2_rd1__rd2_plddts.pdf'
##################################################################
pdf = PdfPages(outfilename)
#plot them
fig = plt.figure()
ax = fig.add_subplot()
#
ax.hist(aplddts,color='r',alpha=0.5)#bins=int(len(allscores)/20)
# ax.hist(aplddts2,color='b',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts3,color='g',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts4,color='y',alpha=0.5)#bins=int(len(allscores)/20)
ax.hist(aplddts5,color='blue',alpha=0.5)#bins=int(len(allscores)/20)

#
ax.axvline(np.mean(aplddts),c='r',linestyle='dashed')
# ax.axvline(np.mean(aplddts2),c='b',linestyle='dashed')
ax.axvline(np.mean(aplddts3),c='g',linestyle='dashed')
ax.axvline(np.mean(aplddts4),c='yellow',linestyle='dashed')
ax.axvline(np.mean(aplddts5),c='blue',linestyle='dashed')

ax.text(0.9,0.9,'Rd1 mean = '+str(np.mean(aplddts))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
# ax.text(0.9,0.95,'OG mean = '+str(np.mean(aplddts2))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
ax.text(0.9,0.85,'Rd1 bests mean = '+str(np.mean(aplddts3))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
ax.text(0.9,0.8,'rd2 mean = '+str(np.mean(aplddts4))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='yellow', fontsize=8)
ax.text(0.9,0.75,'Rd2 bests mean = '+str(np.mean(aplddts5))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
#
ax.set_title('Average plddt for each Sequence across 5 Models')
ax.set_ylabel('Frequency')
ax.set_xlabel('Plddt')
pdf.savefig()
plt.clf()
#
pdf.close()



NOW GETTING BEST DESIGNS TO A NEW DIRECTORY TO THEN SUBMIT 4 RELAX
'''
best_designs2={}
for key in dd2.keys():
    l=dd2[key]
    l=sorted(l, reverse=True,key=lambda first: first[0])
    for x in range(10):
        for u in l:
            if float(u[0])<80.0:
                l.remove(u)
    try:
        best_designs2[key]=[l[0],l[1],l[2],l[3],l[4]]
    except:
        try:
            best_designs2[key]=[l[0],l[1],l[2],l[3]]
        except:
            try:
                best_designs2[key]=[l[0],l[1],l[2]]
            except:
                try:
                    best_designs2[key]=[l[0],l[1]]
                except:
                    try:
                        best_designs2[key]=[l[0]]
                    except:
                        pass
best_designs=best_designs2
len(list(best_designs.keys()))


import os
#
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/tocf'
directory_prefix='r1bd'
outputdirectory='/wynton/home/kortemme/cgalvin/ntf2/round2_best'
os.makedirs(outputdirectory,exist_ok=True)
#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
#
l=[]
for p in prediction_pdb_paths:
    cstrc_name=p.split('/')[-1].split('_')[2]
    if cstrc_name in list(best_designs.keys()):
        for e in best_designs[cstrc_name]:
            #
            best_d_id=str(e[1])
            best_model_id=str(e[2])
            #
            currstrc_design_number=str(p.split('/')[-1].split('_')[-6])
            currstrc_model=str(p.split('/')[-1].split('_')[-3])
            #
            if currstrc_design_number==best_d_id:
                if currstrc_model==best_model_id:
                    cmd='cp '+p+' '+os.path.join(outputdirectory,p.split('/')[-1])
                    # os.system('cp '+p+' '+os.path.join(outputdirectory,p.split('/')[-1]))
                    l.append(cmd)







'''
RELAXING AF2 MODELS OF BEST ROUND2 MPNN DESIGNS
in directory with my best strc...
'''
import os
# os.mkdir('relaxed')
pdbnames=[i for i in os.listdir() if i[-3:]=='pdb']
ofile=open('rr2b.sh','w') ####################enter name
ofile.write('#!/bin/bash')
ofile.write('\n')
ofile.write('tasks=(0\n')
for strc in pdbnames[:-1]:
    ofile.write('       '+strc+'\n')
ofile.write('       '+pdbnames[-1]+')')
ofile.write('\n~/main/source/bin/relax.default.linuxgccrelease -ignore_unrecognized_res -s ${tasks[$SGE_TASK_ID]} -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -in:file:fullatom -relax:fast -out:path:all relaxed')
ofile.write('\nqstat -j "$JOB_ID"')
ofile.close()
print(str(len(pdbnames)))

'''
419
mkdir cluster_output
qsub -cwd -t 1-419 -l mem_free=2G -o cluster_output -e cluster_output rr2b.sh


MAKE POS FILES
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
PREPPING ROSETTA SCORE SCRIPT TO RUN ON rd2
'''
import os
#first copy all the strc to a single directory
# all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed'
# directory_prefix='r2b'
# #
# toscoredir='rd2scores'
# #
# input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
# #
# submitlist=[]
# for id in input_directories:
#     for i in os.listdir(id):
#         if i[-3:]=='pdb':
#             submitlist.append(os.path.join(id,i))
# #
# os.makedirs(toscoredir,exist_ok=True)
# for f in submitlist:
#     os.system('cp '+f+' '+os.path.join(toscoredir,f.split('/')[-1]))
'''
now submit that shit
'''

import os
# os.mkdir('relaxed')
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
qsub -cwd -t 1-419 -l mem_free=2G -o cluster_output -e cluster_output r2rm.sh

main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol scaffold_metrics.xml -parser:view -run:preserve_header -ignore_unrecognized_res -s XXXX.pdb -scorefile_format json -out:file:score_only scores.json


'''



'''
analyze the scores

/wynton/home/kortemme/cgalvin/ntf2/ogscores/ogscores.json
/wynton/home/kortemme/cgalvin/ntf2/round2_best/r2scores.json
        i think maybe i scored the unrealxed structures by accident?
        resubmit .....

/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxedr2scores.json
'''
import json
import os
#analysis of scores from json file
sfname='/wynton/home/kortemme/cgalvin/ntf2/ogscores/ogscores.json'

f=open(sfname,'r')
lines=[line for line in f.readlines()]
f.close()
'''
ok some lines were corrupted cus they put 2 dicts on one line somehow
lemme see if i can fix this
'''
flines=[]
for line in lines:
    dendi=[]
    for ii,i in enumerate(line):
        if i=='}':
            dendi.append(ii)
    if len(dendi)>1:
        firstline=line[:dendi[0]+1]
        secline=line[dendi[0]+1:dendi[1]+1]
        # lines.remove(line)
        flines.append(firstline)
        flines.append(secline)
    else:
        flines.append(line)
#
scores=[]
for i,line in enumerate(flines):
    try:
        scores.append(json.loads(line))
    except:
        print(i)

#add the custom metrics of hphobesasa/totalsasa
dontnorm=['ps','rescount','oversat','acc','buns_bb_heavy','buns_sc_heavy',
          'ref','yhh_planarity','np_over_tsasa']
for d in scores:
    h=d['hphobe_sasa']
    t=d['total_sasa']
    r=float(h)/t
    d['np_over_tsasa']=r
    nres=int(d['rescount'])
    toadd=[]
    for key in d.keys():
        if key not in dontnorm:
            v=d[key]
            try:
                nv=float(v)/nres
                nvs=key+'_perres'
                toadd.append((nvs,nv))
            except:
                pass
    for a,b in toadd:
        d[a]=b




#analysis of scores from json file
sfname2='/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/r2scores.json'

f=open(sfname2,'r')
lines=[line for line in f.readlines()]
f.close()
'''
ok some lines were corrupted cus they put 2 dicts on one line somehow
lemme see if i can fix this
'''
flines=[]
for line in lines:
    dendi=[]
    for ii,i in enumerate(line):
        if i=='}':
            dendi.append(ii)
    if len(dendi)>1:
        firstline=line[:dendi[0]+1]
        secline=line[dendi[0]+1:dendi[1]+1]
        # lines.remove(line)
        flines.append(firstline)
        flines.append(secline)
    else:
        flines.append(line)
#
scores2=[]
for i,line in enumerate(flines):
    try:
        scores2.append(json.loads(line))
    except:
        print(i)

#add the custom metrics of hphobesasa/totalsasa
dontnorm=['ps','rescount','oversat','acc','buns_bb_heavy','buns_sc_heavy',
          'ref','yhh_planarity','np_over_tsasa']
for d in scores2:
    h=d['hphobe_sasa']
    t=d['total_sasa']
    r=float(h)/t
    d['np_over_tsasa']=r
    nres=int(d['rescount'])
    toadd=[]
    for key in d.keys():
        if key not in dontnorm:
            v=d[key]
            try:
                nv=float(v)/nres
                nvs=key+'_perres'
                toadd.append((nvs,nv))
            except:
                pass
    for a,b in toadd:
        d[a]=b


#repeats in here for some reason
sc2=[]
good=[i for i in os.listdir('/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed') if i[-3:]=='pdb']
seen=[]
for d in scores2:
    if d['decoy'] not in seen:
        sc2.append(d)
        seen.append(d['decoy'])
scores2=sc2

terms=[key for key in scores[0].keys()]

# #make a pdf showing all score distributions
# ######
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.pyplot as plt
# def plot_dists(terms,scores,outfilename):
#     pdf = PdfPages(outfilename)
#     for term in terms:
#         allscores=[]#score
#         if term != 'decoy':
#             for d in scores:
#                 allscores.append(float(d[term]))
#         fig,ax=plt.subplots()
#         if len(allscores)!=0:
#             ax.hist(allscores)#bins=int(len(allscores)/20)
#             ax.set_title(term)
#             ax.set_ylabel('frequency')
#             ax.set_xlabel('score')
#             pdf.savefig()
#             plt.clf()
#     pdf.close()
#
# plot_dists(terms,scores,'ntf2ogrosettametrics.pdf')


'''
/wynton/home/kortemme/cgalvin/ntf2/ogscores/ogscores.json
'''
#plot avg plddt (across 5 models) for each structure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
#gotta fix this but its for stacked hists once i get the design metrics
pdfname='stackedhists_og_rd2.pdf'
pdf = PdfPages(pdfname)
#
#gonna do just for sub micromolar vs micromolar and up binders
for term in terms:
    og=[]
    des=[]
    for d in scores:
        try:
            og.append(float(d[term]))
        except:
            pass
    for d2 in scores2:
        try:
            des.append(float(d2[term]))
        except:
            pass
    try:
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(og,color='red',alpha=0.3)
        # ax.hist(x,color='green',alpha=0.3)
        ax.hist(des,color='blue',alpha=0.3)
        ax.axvline(np.mean(og),c='r',linestyle='dashed')
        # ax.axvline(np.mean(x),c='g',linestyle='dashed')
        ax.axvline(np.mean(des),c='b',linestyle='dashed')
        ax.set_xlabel(term)
        ax.set_ylabel('Count')
        ax.set_title(term)
        ax.text(0.9,0.9,'OG = '+str(np.mean(og))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
        # ax.text(0.9,0.8,'nm mean = '+str(np.mean(x))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
        ax.text(0.9,0.7,'Design = '+str(np.mean(des))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
        pdf.savefig()
        plt.clf()
        plt.close()
    except:
        pass
pdf.close()

'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/ntf2/round2_best/relaxed/stackedhists_og_rd2.pdf ~/desktop/stackedhists_og_rd2.pdf




















alright thats chillin, next is
SUBMIT ROUND2 BEST SEQUENCES FOR COLABFOLD
then i can relax and score them to compare with these scores ^^^^
    and also i will match with those relaxed structures
        NEED TO GENERATE THE MOTIFS THAT I WILL WANT TO MATCH INTO THESE SCAFFOLDS!

'''

import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/ntf2/round1_best/relaxed'
directory_prefix='r1b'
ntopseqstochoose=10
selectfastasoutputdir='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/tocf'
os.makedirs(selectfastasoutputdir,exist_ok=True)
#
input_directories=[os.path.join(all_strc_dir,i) for i in os.listdir(all_strc_dir) if os.path.isdir(os.path.join(all_strc_dir,i))==True and i[:len(directory_prefix)]==directory_prefix]
#
allresl=[]
for id in input_directories:
    resultfspath=os.path.join(id,'_'.join([id.split('/')[-1],'output']),'seqs')
    resl=[i for i in os.listdir(resultfspath) if i[-3:]=='.fa']
    for y in resl:
        allresl.append(os.path.join(resultfspath,y))
'''
analyze the mpnn results and select sequences for af2 modeling

    are the score dists different for different input bbs?
    des score dist vs og score dist vs lowest score dist
    seq recov dist
    seq entropy of designs
take top ten scoring strc for each to submit to colabfold
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
for key in add.keys():
    selectseqs=[i[2] for i in add[key][:ntopseqstochoose]]
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
#
import os
nfastasperjob=2
directory_prefix='r1bd'
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/tocf'
shf_prefix='r1bd'
l3=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if i[-6:]=='.fasta']
#
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
import os
nfastasperjob=2
directory_prefix='r1bd'
allfastasdir='/wynton/home/kortemme/cgalvin/ntf2/rd2mpnnresults/tocf'
shf_prefix='r1bd'
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
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output r1bd_197_.sh
'''
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
jobss.remove('r1bd_197_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output '+j
    os.system(cmd)
