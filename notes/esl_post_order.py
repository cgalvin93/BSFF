'''
i wanna use original design/fiultering aprams that i used to get
the designs i ordered

-BUT this time i want to disallow methionine during binding site design
-then i want to try mpnn informed design without allowing native res (unless it
in mpnn profile)
-then i want to apply the refinement script to filtered des to fix as many hbonds as possible
-and finally i want to ue fragqual for filtering

'''
#in /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean
#getting the matches i want
#doesnt have bad res:
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

#has at least one aro:
import os
good2=[]
for i in good:
    motif=i.split('_')[2]
    chars=list(motif)
    if 'W' not in chars:
        if 'Y' not in chars:
            if 'F' not in chars:
                pass
    else:
        if i not in good2:
            good2.append(i)

'''
In [12]: len(good2)
Out[12]: 11872

In [13]: len(l)
Out[13]: 45673

In [14]: len(good)
Out[14]: 30371

In [15]: len(good2)
Out[15]: 11872
'''

os.makedirs('esl4rm_1_bestmatches')
c=1
for i in good2:
    os.system('cp '+i+' esl4rm_1_bestmatches/'+i)
    c+=1
    print(c)

'''
mv esl4rm_1_bestmatches ~/esl4rm_1_bestmatches
'''

























'''
remake resfiles but now trying to filter surface res that are close to lig
'''


#resfile generation
#################################################################################
#################################################################################
#################################################################################
#################################################################################
import os
#################################################################################
inputs_dir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
sfname='srmr.sh'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
11872
'''
count=1
idx=[j for j in range(0,len(matches),5)]  #how many matches per job
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
        cmdd='time python ~/BSFF/tools/match_resfile2.py '+matchpath+' '+prm1+' 4.5 genpot'
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
now moving cst files from gp/clean directory to my new selected matches directory
so that it will be easier to submit jobs (will use single params from og dir,
then resfiles,pdbs,csts will all be in my new selected match dir)

cd /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean
'''
import os
c=0
matches1=[i.split('.')[0] for i in os.listdir('/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches') if i[-3:]=='pdb']
for i in os.listdir():
    if i[-3:]=='cst':
        if i.split('.')[0] in matches1:
            os.system('cp '+i+' /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/'+i)
            c+=1
            print(c)

















'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches

just gonna move the params file here actually so i can delete old designs
import os
inputs_dir='/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]

os.system('mv '+prm1+' /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/'+prm1.split('/')[-1])

ls -d */


mkdir enzdes
cd enzdes
'''

import os
#get a params file from old directory
inputs_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
# ndesignseachmatch='50'
ndesignseachmatch='20'
sfname='enzdes.sh'

#
inputs_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches'
matches1=[i for i in os.listdir('/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches') if i[-3:]=='pdb']
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
output_prefix='ed'+str(c)
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
     '-enzdes:cst_design True','-enzdes:design_min_cycles 3', '-enzdes:lig_packer_weight 2.5',
     '-enzdes:cst_min True',
     # '-enzdes:detect_design_interface',
     # '-enzdes:cut1','6.0',
     # '-enzdes:cut2','8.0',
     # '-enzdes:cut3','10.0',
     # '-enzdes:cut4','12.0',
     # '-score:set_weights','hbond_bb_sc','2.0',
     # '-score:set_weights','hbond_sc','2.0',
     # '-score:set_weights','hbnet','1.0',
     # '-score:set_weights','voids_penalty','0.5',
     # '-score:set_weights','buried_unsatisfied_penalty','0.5',
     # '-score:set_weights','approximate_buried_unsat_penalty','5.0',
     '-packing:soft_rep_design True',
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
~/main/source/bin/enzyme_design.linuxgccrelease -s /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/UM_6_M14W13S70S88_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_hybrid_16_1_0001_clean.pdb -in:file::extra_res_fa /wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/UM_13_F58S87T51S110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_11_1_0001_ligH_bcc.params -resfile /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/UM_6_M14W13S70S88_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_hybrid_16_1_0001_clean.resfile -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt True -enzdes:bb_min True -enzdes:chi_min True -enzdes:cst_design True -enzdes:design_min_cycles 3 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min True -packing:soft_rep_design True -relax:script InterfaceDesign2019 -no_packstat_calculation True -out:nstruct 20 -packing:linmem_ig 10 -corrections::gen_potential -load_PDB_components False -enzdes:cstfile /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/UM_6_M14W13S70S88_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_1_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_hybrid_16_1_0001_clean.cst -out:prefix ed1
'''
import os
sfname='enzdes.sh'
cmd='qsub -cwd -t 1-11872 -l mem_free=2G -o cluster_output -e cluster_output '+sfname
os.system(cmd)

'''
beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  703.60 GiB| 1000.00 GiB||  4547442|unlimited



set jobs at 11.56 am on 12.30
'''

































'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes
'''
import os
#################################################################################
inputs_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches'
allparams=[os.path.join(inputs_dir,i) for i in os.listdir(inputs_dir) if i[-6:]=='params']
prm1=allparams[0]
sfname='qa.sh'
scorefile_name='enzdesanalysis.json'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
204759

'''
count=1
idx=[j for j in range(0,len(matches),50)]  #how many matches per job
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
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes


'''





'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE



NOW THAT I HAVE REFINEMENT SCRIPT, WHICH IS CAPABLKE OF
FIXING AT LEAST ONE HYDROGEN BOND, IM GONNA SAY I CAN ALLOW ONE
UNSAT DURING THE FILTERING STAGES
    yea cus whether its a lig or ptn unsat i can address
    if lig then maybe can find hb partner
    if ptn then maybe can mutate the polar res to a np
'''
#analysis of scores from json file
sfname='enzdesanalysis.json'
##################
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
63248
94463
94461
31639
131638


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

f1=return_filtered(scores,'buns2interface','<',1.0)
f2=return_filtered(f1,'lighphobesasa','<',25.0)
f3=return_filtered(f2,'hbtolig','>',3.0)
f4=return_filtered(f3,'ligoversat','<',0.0)
# f5=return_filtered(f4,'bsE_per_res','<',-2.0)
# plot_dists(terms,f3,'esl4rmfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
# print(len(f5))
'''
131638
57852
23261
11559
11547
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

'''
cd filtered

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/1np3hb/1np/buried/genpot/clean/ed1/esl_4rm.pdf ~/desktop/esl_4rm.pdf

ADD CORRECT DIRECTORY TO BOTTOM OF CODE HERE (EXCEPT LINE)
cd filtered
pwd
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered

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
        os.chdir('/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered')
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
import os
#################################################################################
params=[i for i in os.listdir() if i[-6:]=='params']
prm1=params[0]
sfname='extrascore.sh'
#################################################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(matches))
'''
11522
'''
count=1
idx=[j for j in range(0,len(matches),3)]  #how many per job
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
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run

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
10
11414
11414
5
11417


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
f1=return_filtered(scores,'buns2interface','<',1.0)
f2=return_filtered(f1,'boltz','<',-0.3)
f3=return_filtered(f2,'ddg','<',-22.0)
f4=return_filtered(f3,'contact_molsurf','>',165)
f5=return_filtered(f4,'hbtolig','>',2.0)
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
11417
8739
1131
833
533
523
455
455
----
291
246
123
122
109
82
82
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
455
19
167
245
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
RUNNING MPNN ON MATCHES TO GENERATE SEQUENCE PROFILE
design at all positions not in binding site first shell

EXAMPLE OF FIXED POSITIONS DICT
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}

/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2
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
'''
echo "
import os
from pyrosetta import*
import json
init('-load_PDB_components False -gen_potential')
##########
allparams=[os.path.join('/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run',i) for i in os.listdir('/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run') if i[-6:]=='params']
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
    #################
    print(first_shell_res)
    for residd in first_shell_res:
        individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        individual_res_selector.set_index(residd)
        individual_res_selector.apply(p)
        sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
        sasa_metric.set_residue_selector(individual_res_selector)
        this_res_sasa = sasa_metric.calculate(p)
        if this_res_sasa>=12.:
            first_shell_res.remove(residd)
    print(first_shell_res)
    #this is one of those stupid instances wbere i have to filter a list twice
    #for no apparent reason, pisses me off -.-
    for residd in first_shell_res:
        individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        individual_res_selector.set_index(residd)
        individual_res_selector.apply(p)
        sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
        sasa_metric.set_residue_selector(individual_res_selector)
        this_res_sasa = sasa_metric.calculate(p)
        if this_res_sasa>=12.:
            first_shell_res.remove(residd)
    print(first_shell_res)
    ###############
    for resnum in first_shell_res:
        if resnum not in tofreeze:
            tofreeze.append(resnum)
    tfdata[pdb]=tofreeze

print(len(list(tfdata.keys())))
json.dump(tfdata,open(ofjsonname,'w'))
print('json output')


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


">maked.py

echo "
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python maked.py
">maked.sh



qsub -cwd maked.sh





cd /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2


mv mpnn_params.json mpnn
cd mpnn
'''

#making the folders with the pdbs
#in this case i have to do 1 per job cus each file will
#have different fixed positions
import json
import os
########
ofjsonname='mpnn_params.json'
#########
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
all_strc_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn'
directory_prefix='mpnndesign'
shf_prefix='mpnn_'

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
# jobss.remove('mpnn__1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)

'''
~30 minutes per job only












consolidating mpnn results
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn'
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
447




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
ofjsonname='mpnn_params.json'
#called pdbsdir, using to get inputs for stability redesign (strc and params)
paramsdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run'
#
allparams=[os.path.join(paramsdir,i) for i in os.listdir(paramsdir) if i[-6:]=='params']
prm1=allparams[0]
pdbsdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2'

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
    ############was removing second shell res from redesign but idk why
    # p.update_residue_neighbors()
    # ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
    # neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
    # neighborhood_selector_bool = neighborhood_selector.apply(p)
    # neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    # first_shell_res=list(neighborhood_residues_resnums)
    # #
    # second_shell=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(neighborhood_selector,10.,False)
    # ss_bool = second_shell.apply(p)
    # ss_resnums=pyrosetta.rosetta.core.select.get_residues_from_subset(ss_bool)
    # ss_res=list(ss_resnums)
    # #
    # for i in des_res:
    #     if i in ss_res:
    #         ss_res.remove(i)
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
        # elif i in ss_res:
        #     s=str(p.pdb_info().pose2pdb(i))+'NATAA'
        #     ofile.write(s+'\n')
        # elif i==p.total_residue():
        #     s=str(p.pdb_info().pose2pdb(i))+'NATAA'
        #     ofile.write(s+'\n')
        else:
            s=str(p.pdb_info().pose2pdb(i))+'NATAA'
            ofile.write(s+'\n')
    ofile.close()

'''
FUCKED UP AND FORGOT THAT I DID THIS THE RIGHT WAY BELOW
BUT WHATEVER I THINK I ALSO CHANGED THE ABOVE CODE TO DO IT RIGHT ALSO

'''



'''
cd resfiles

copy the relevant pdbs and params to this directory for stability focused design
'''

import os
#########
paramsdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run'
pdbsdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2'
############
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
447
447
712


okay having an issue here where there are too many params file,
i guess ill remove the extraneoius ones but wtf is going on with this?
                SKETCHY???????
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
447
447
447
'''





'''
FASTDESIGN with 3bop weight 5
no intup
some filters but no ddg/boltz
'''
# FASTDESIGN with beta nov16 sf
import os
########
allparams=[i for i in os.listdir() if i[-6:]=='params']
prm1=allparams[0]
sfname='fdmpnn.sh'
ndesignseachmatch='10'
outputdirectory='fd_mpnn'
scorefile_name='esl_fd_mpnn.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(10):     #this way i can spread across more nodes if neccessary
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
         '-beta_nov16',
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
447


~/main/source/bin/rosetta_scripts.default.linuxgccrelease -s ed1UM_18_F112W57T31S107_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean__DE_10_oglig_0001.pdb -parser:protocol ~/BSFF/tools/fd3bop_mpnn.xml -parser:view -run:preserve_header -in:file::extra_res_fa ed1UM_2_M58S109S95W93_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_1_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_12_oglig_0001.params -resfile ed1UM_18_F112W57T31S107_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_hybrid_95_1_0001_clean__DE_10_oglig_0001.resfile -ex1 -ex2 -extrachi_cutoff 0 -use_input_sc -out:nstruct 10 -beta_nov16 -load_PDB_components False -out:path:all fd_mpnn -out:prefix fd1 -scorefile_format json -out:file:scorefile esl_fd_mpnn.json
'''
import os
sfname='fdmpnn.sh'
jobss=[i for i in os.listdir() if sfname in i]
for j in jobss:
    cmd='qsub -cwd -t 1-447 -l mem_free=4G -o cluster_out -e cluster_out '+j
    os.system(cmd)



'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn
'''




'''

SCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORESCORE

'''




#analysis of scores from json file
sfname='esl_fd_mpnn.json'
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
264
43693
43693
132
43849
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

# plot_dists(terms,scores,'esl_fd3bmpnn.pdf')

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
f1=return_filtered(scores,'buns2interface','<',1.0)
f2=return_filtered(f1,'contact_molsurf','>',165)
f3=return_filtered(f2,'hbtolig','>',2.0)
f4=return_filtered(f3,'shape_comp','>',0.65)
f5=return_filtered(f4,'lighphobesasa','<',40.0)
f6=return_filtered(f5,'packstat','>',0.6)
f7=return_filtered(f6,'buns_bb_heavy','<',5)
f8=return_filtered(f7,'buns_sc_heavy','<',1)
f9=return_filtered(f8,'ligoversat','<',0)
f10=return_filtered(f9,'oversat','<',0)
f11=return_filtered(f10,'exphyd','<',1300)
f12=return_filtered(f11,'cav','<',140)
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
'''
43849
38332
34032
33355
30971
30965
19374
16267
8280
8241
7351
5836
5801
'''
import os
filtered_strc=[]
for d in f12:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered2',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i+'.pdb filtered2/'+i+'.pdb')
    ########################################
# os.system('mv *.pdf filtered2/')

'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2

DIVERSITY AMONGST FILTERED DESIGNS

'''
import os
des=[i for i in os.listdir() if i[-3:]=='pdb']


motifs=[]
scaffs=[]
pairs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-10]
    scaffold=('_').join(s.split('_')[7:-10])
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
5801
18
131
173
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
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2


I think the next move should be to apply the hbond refinement script to all of these
followed by fragment quality filtering
and then finally alphafold


TESTING SOME CHANGES TO THE SCRIPT
time python3 ~/desktop/BSFF/tools/refine_filtered_designs.py ~/desktop/esl.params ed1UM_9_M98T110S58S61_1_relaxed_relaxed_5tpj_39178_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_9_unrelaxed_model_2_rank_1_0001_hybrid_39_1_0001_clean__DE_24_oglig_0001.pdb solutions
time python3 ~/desktop/BSFF/tools/refine_filtered_designs.py ~/desktop/esl.params ed1UM_63_M89Y42S13N112_1_relaxed_relaxed_5tpj_371721_design_5_unrelaxed_model_3_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_2_rank_1_0001_hybrid_5_1_0001_clean__DE_20_oglig_0001.pdb solutions
time python3 ~/desktop/BSFF/tools/refine_filtered_designs.py ~/desktop/esl.params ed1UM_35_M96Y31S65N61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid_5_1_0001_clean__DE_7_oglig_0001.pdb solutions
time python3 ~/desktop/BSFF/tools/refine_filtered_designs.py ~/desktop/esl.params ed1UM_10_L58Y87Q36T110_1_relaxed_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid_1_1_0001_clean__DE_4_oglig_0001.pdb solutions
time python3 ~/desktop/BSFF/tools/refine_filtered_designs.py ~/desktop/esl.params ed1UM_10_F95T88T17Y38_1_relaxed_relaxed_5tpj_128842_design_5_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_4_rank_1_0001_design_2_unrelaxed_model_1_rank_1_0001_hybrid_10_1_0001_clean__DE_29_oglig_0001.pdb solutions





INDIVIDUAL MUTATIONS FOR REFINEMENT
'''

import os
#get a params file from old directory (og params)
params_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run'
allparams=[os.path.join(params_dir,i) for i in os.listdir(params_dir) if i[-6:]=='params']
prm1=allparams[0]
sfname='refine.sh'
#
inputs_dir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2'
matches=[i for i in os.listdir(inputs_dir) if i[-3:]=='pdb']
print(len(matches))
#
#########
###########
try:
    os.system('rm -r cluster_output')
except:
    pass
os.makedirs('cluster_output',exist_ok=True)
###########
sf=open(sfname,'w')
sf.write('#!/bin/bash')
sf.write('\nsource ~/anaconda3/etc/profile.d/conda.sh\n')
sf.write('conda activate pyr37\n')
sf.write('\n')
sf.write('tasks=(0\n')
for match in matches[:-1]:
    sf.write('       '+str(match.strip('.pdb'))+'\n')
sf.write('       '+matches[-1].strip('.pdb')+')')
sf.write('\n')
cmd=['time python3',
     '~/BSFF/tools/refine_filtered_designs.py',
     prm1,
     '${tasks[$SGE_TASK_ID]}.pdb',
     'refined']
sf.write((' ').join(cmd))
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
print(len(matches))
'''
5801


time python3 ~/BSFF/tools/refine_filtered_designs.py /wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/ed1UM_27_M99S111S88W59_1_relaxed_relaxed_5tpj_144267_design_10_unrelaxed_model_2_rank_1_0001_design_8_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_3.params fd5ed1UM_5_L52S37T91W18_1_relaxed_relaxed_5tpj_294053_design_9_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_107_1_0001_clean__DE_4_oglig_0001_0010.pdb refined

'''
import os
sfname='refine.sh'
cmd='qsub -cwd -t 1-5801 -l mem_free=2G -o cluster_output -e cluster_output '+sfname
os.system(cmd)

'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined

In [1]: import os

In [2]: len(os.listdir())
Out[2]: 5231
'''






























'''
FRAGQUAL ANALYSIS

/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined

'''
import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
print(len(l))

'''
make fasta files
'''
from pyrosetta import *
init('-ignore_unrecognized_res -load_PDB_components False') #because of this the ligand shouldnt be there
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())
    print(sequence)
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

c=1
for a in l:
    print(c)
    c+=1
    try:
        fasta(a)
    except:
        print('failure')

os.makedirs('fastas',exist_ok=True)
l2=[i for i in os.listdir() if i[-6:]=='.fasta']
for i in l2:
    os.system('mv '+i+' fastas/'+i)

############################
##############
##############
##############
##############
'''
SUBMIT FRAGMENT PICKING JOBS
'''
##############
##############
##############
##############
############################
import os
l=[i for i in os.listdir() if i[-5:]=='fasta']

# l.remove('fd9ed1UM_9_M55W62S109S51_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_16_1_0001_clean__DE_16_oglig_0001_0002_unrefined.fasta')


c=1
for fastaf in l:
    print(c)
    c+=1
    s=fastaf.split('.')[0]
    inf=os.path.join(os.getcwd(),fastaf)
    cmd='time python ~/BSFF/tools/fragment_picker/workingfragpick.py '+inf+' '+s
    os.system(cmd)

'''
on a test file it looks like job only takes about 12.4 G memory max, so Im
gonna change the 'submit...wynton.sh' script to only request 40G, which should
be way more than enough

IF I EVER HAVE LARGER PROTEINS, I SHOULD BE SURE TO CHANGE THIS TO A HIGHER VALUE
AS NECCESSARY


ONLY ABOUT 25-30 MIN for jobs of these size ptns


output has name inpuA.200.9mers




NOW APPLYING FRAGMENT QUALITY FILTER
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined

'''
import os

###################
####################
jsonoutputdir='fragment_scores'
os.makedirs(jsonoutputdir,exist_ok=True)
###################
###################

pdbs=[i for i in os.listdir() if i[-3:]=='pdb']

inputs_dict={}
for pdb in pdbs:
    pdbn=pdb.split('.')[0]
    fragments_input_path=os.path.join('fastas',pdbn,'fragments','inpuA.200.9mers')
    if os.path.exists(fragments_input_path):
        # print('yea')
        inputs_dict[pdb]=fragments_input_path

shcount=1
for key in inputs_dict.keys():
    fraginputpath=inputs_dict[key]
    fragjsonoutname=os.path.join(jsonoutputdir,'fragment_scores_'+str(shcount)+'.json')
    of=open('fragqual_'+str(shcount)+'.sh','w')
    of.write('#!/bin/bash')
    of.write('\nsource ~/anaconda3/etc/profile.d/conda.sh\n')
    of.write('conda activate pyr37\n')
    of.write('\n')
    cmd=['time python3',
         '~/BSFF/tools/fragment_picker/fragqual_filter.py',
         fraginputpath,
         key,
         fragjsonoutname]
    of.write((' ').join(cmd))
    of.write('\nqstat -j "$JOB_ID"')
    of.close()
    shcount+=1

if not os.path.exists('cluster_output'):
    os.makedirs('cluster_output',exist_ok=True)

import os
l=[i for i in os.listdir() if 'fragqual_' in i]
for i in l:
    cmd='qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output '+i
    os.system(cmd)

'''
Your job 959088 ("fragqual_5231.sh") has been submitted


/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined


ONLY TOOK LIKE 2 HRS MAYBE



ANALYZING FRAGQUAL FILTER RESULTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEED TO CHANGE FRAGQUAL SCRIPT TO OUTPUT JSON IN A CLEANER WAY SO I CAN ACTUALLY ACCESS ALL
THE SCORES
i think the easiest way here is gonna be to just output a different json
file for each run, then i can later concatenate them all
fuck it -.-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/fragment_scores

'''
import json
import os
#analysis of scores from json file
jsonoutputdir='fragment_scores'
sfs=[os.path.join(jsonoutputdir,i) for i in os.listdir(jsonoutputdir) if i[-4:]=='json']

scores=[]
for sf in sfs:
    f=open(sf,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        scores.append(json.loads(line))
##################
# import json
# starts=[]
# ends=[]
# f=open(sfname,'r')
# lines=[line for line in f.readlines()]
# f.close()
# for line in lines:
#     for ind,char in enumerate(line):
#         if char=='{':
#             starts.append(ind)
#         elif char=='}':
#             ends.append(ind)
#
# l=[len(starts),len(ends)]
# scores=[]
# for x in range(min(l)):
#     try:
#         td=lines[0][starts[x]:ends[x]+1]
#         try:
#             scores.append(json.loads(td))
#         except:
#             if len(td)>1:
#                 print('\n\n\n\n\n\n\n')
#                 print(td)
#                 print('\n\n\n\n\n\n\n')
#     except:
#         pass
terms=list(scores[0].keys())
print(len(scores))
'''
5231


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

plot_dists(terms,scores,'esl_fragscores.pdf')

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

f1=return_filtered(scores,'max_min_rmsd','<',1.3)
f2=return_filtered(f1,'design_avg_rmsd','<',2.0)
#these are the values for 6w9o (xingjies solved structure)
# plot_dists(terms,f3,'esl4rmfilt.pdf')
print(len(scores))
print(len(f1))
print(len(f2))
# print(len(f5))
'''
5231
337
337

'''
filtered_strc=[]
for d in f2:
    nn=d['description']
    filtered_strc.append(nn)
os.makedirs('filtered',exist_ok=True)
print(len(filtered_strc))
c=0
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
    c+=1
    print(str(c))
os.system('mv *.pdf filtered/')


'''
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered


'''
import os
des=[i for i in os.listdir() if i[-3:]=='pdb']


motifs=[]
scaffs=[]
pairs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-11]
    scaffold=('_').join(s.split('_')[7:-12])
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
337
6
15
15



OKAY SO THE FINAL THING TO DO IS
AF2 AND
EXTRASCORE DOUBLE CHECKING
and then finally
to somehow select like a best design from each of the unique 15
(assuming they all make it thru af2/extrascores)
I COULD ALSO CONSIDER FORWARD FOLDING? PROBABY AF2 IS JUST AS USEFUL EH
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


cd cf/fastas
mkdir cluster_output

'''
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


/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas

'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='filtdesigns'
allfastasdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas'
nonaf2desdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered'

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

plddt_threshold=85.0
carmsd_threshold=2.0
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
'''
144
'''


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
os.system('cp *.pdf af2filtered2')
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
print(len(des))
print(len(set(motifs)))
print(len(set(scaffs)))
print(len(set(pairs)))
'''
144
5
10
10
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
/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas/af2filtered2
'''


#run the extra scoring
import os
paramsdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run'
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
sf.write('\n/wynton/home/kortemme/cgalvin/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s=${tasks[$SGE_TASK_ID]} -parser:protocol /wynton/home/kortemme/cgalvin/BSFF/tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+paramspath+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only scores.json -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -beta_nov16 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
#
print(len(pdbs))
'''





qsub -cwd -t 1-144 -l mem_free=2G extrascore.sh


/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas/af2filtered2
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
0
144
144
0
144


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
f2=return_filtered(f1,'boltz','<',-0.27)
f3=return_filtered(f2,'ddg','<',-20.0)
f4=return_filtered(f3,'contact_molsurf','>',160)
f5=return_filtered(f4,'hbtolig','>',3.0)
f6=return_filtered(f5,'shape_comp','>',0.65)
f7=return_filtered(f6,'lighphobesasa','<',20.0)
f8=return_filtered(f7,'packstat','>',0.55)
f9=return_filtered(f8,'buns_bb_heavy','<',5)
f10=return_filtered(f9,'buns_sc_heavy','<',1)
f11=return_filtered(f10,'ligoversat','<',0)
f12=return_filtered(f11,'oversat','<',0)
f13=return_filtered(f12,'exphyd','<',1200)
f14=return_filtered(f13,'cav','<',140)
#
# plot_dists(terms,f14,'final_filtered_binders.pdf')
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
144
131
131
131
129
121
121
121
121
116
114
114
111
108
108
'''
import os
filtered_strc=[]
for d in f14:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i[:-5]+'.pdb filtered/'+i+'.pdb')
    ########################################
os.system('mv *.pdf filtered/')
'''




DESIGN DIVERSITY


'''



#
import os
des=[i for i in os.listdir() if i[-3:]=='pdb']

motifs=[]
scaffs=[]
pairs=[]
for d in des:
    s=d.split('.')[0]
    motif=s.split('_')[-12]
    scaffold=('_').join(s.split('_')[3:-12])
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
2
5
5
'''
# same_scaff_diff_match={}
# for scaff in set(scaffs):
#     l=[]
#     for d in des:
#         s=d.split('.')[0]
#         scaffold=('_').join(s.split('_')[3:-11])
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
d={}
for motif,scaff in pairs:
    l=[]
    for i in des:
        if motif in i:
            if scaff in i:
                l.append(i)
    d[motif+'_'+scaff]=l

'''
{'52_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid': ['fd3ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0009_refined_0001.pdb',
  'fd6ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0001_refined_0001.pdb',
  'fd3ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0005_refined_0001.pdb',
  'fd2ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0008_refined_0001.pdb',
  'fd10ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0003_refined_0001.pdb'],
 '52_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid': ['fd4ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0007_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd2ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0005_refined_0001.pdb',
  'fd4ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd10ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd2ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0002_refined_0001.pdb',
  'fd10ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0001_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd3ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd2ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0005_refined_0001.pdb',
  'fd7ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0002_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0006_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd3ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0006_refined_0001.pdb',
  'fd1ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0002_refined_0001.pdb',
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0006_refined_0001.pdb',
  'fd10ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0006_refined_0001.pdb',
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd1ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd3ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0008_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0007_refined_0001.pdb',
  'fd4ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0002_refined_0001.pdb',
  'fd2ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0007_refined_0001.pdb',
  'fd2ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd4ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd1ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0005_refined_0001.pdb',
  'fd7ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0001_refined_0001.pdb',
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd7ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0008_refined_0001.pdb',
  'fd8ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0008_refined_0001.pdb',
  'fd8ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0009_refined_0001.pdb',
  'fd10ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd3ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0005_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0002_refined_0001.pdb',
  'fd7ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0008_refined_0001.pdb',
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0003_refined_0001.pdb',
  'fd4ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0004_refined_0001.pdb',
  'fd8ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb',
  'fd6ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0007_refined_0001.pdb',
  'fd9ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0004_refined_0001.pdb'],
 '52_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid': ['fd5ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0010_refined_0001.pdb',
  'fd9ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0005_refined_0001.pdb',
  'fd5ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0008_refined_0001.pdb',
  'fd6ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0001_refined_0001.pdb',
  'fd10ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0002_refined_0001.pdb'],
 '52_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid': ['fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd2ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0005_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0005_refined_0001.pdb',
  'fd8ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0003_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0001_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd2ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0003_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0008_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd8ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0001_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd2ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0003_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0008_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0003_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0003_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0001_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0007_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0008_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0001_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd2ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd8ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd8ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0008_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0005_refined_0001.pdb',
  'fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd8ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd10ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0001_refined_0001.pdb',
  'fd7ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0010_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd4ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0009_refined_0001.pdb',
  'fd1ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0004_refined_0001.pdb',
  'fd3ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd5ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0002_refined_0001.pdb',
  'fd6ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0005_refined_0001.pdb'],
 '62_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid': ['fd6ed1UM_2_M54W89W100T38_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid_62_1_0001_clean__DE_19_oglig_0001_0007_refined_0001.pdb',
  'fd7ed1UM_2_M54W89W100T38_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid_62_1_0001_clean__DE_19_oglig_0001_0003_refined_0001.pdb'],
 '52_1_relaxed_relaxed_5tpj_293640_design_3_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid': ['fd7ed1UM_2_M42S99S95W90_1_relaxed_relaxed_5tpj_293640_design_3_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_2_oglig_0001_0010_refined_0001.pdb',
  'fd10ed1UM_2_M42S99S95W90_1_relaxed_relaxed_5tpj_293640_design_3_unrelaxed_model_1_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_2_oglig_0001_0003_refined_0001.pdb']}

'''


'''
NOW WHAT IS GOING TO BE MY FINAL FILTERING CRITERIA TO SELECT ONE DESIGN FROM EACH
UNIQUE MATCH????????

it should be a stability thing, either alphafold result or fragment quality result
or some combination of both

im gonna say fragment quality is actually more trustworthy

/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/fastas/fragment_scores
'''

import json
import os
#analysis of scores from json file
jsonoutputdir='/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/fragment_scores'
sfs=[os.path.join(jsonoutputdir,i) for i in os.listdir(jsonoutputdir) if i[-4:]=='json']

scores=[]
for sf in sfs:
    f=open(sf,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        scores.append(json.loads(line))
##################
# import json
# starts=[]
# ends=[]
# f=open(sfname,'r')
# lines=[line for line in f.readlines()]
# f.close()
# for line in lines:
#     for ind,char in enumerate(line):
#         if char=='{':
#             starts.append(ind)
#         elif char=='}':
#             ends.append(ind)
#
# l=[len(starts),len(ends)]
# scores=[]
# for x in range(min(l)):
#     try:
#         td=lines[0][starts[x]:ends[x]+1]
#         try:
#             scores.append(json.loads(td))
#         except:
#             if len(td)>1:
#                 print('\n\n\n\n\n\n\n')
#                 print(td)
#                 print('\n\n\n\n\n\n\n')
#     except:
#         pass
terms=list(scores[0].keys())
print(len(scores))
'''
5231


d=.....^^^^^^^
'''

bests=[]

for key in d.keys():
    scores_this_match=[]
    for design in d[key]:
        fragscore_name=('_').join(design.split('_')[:-1])#####################
        for entry in scores:
            if entry['description']==fragscore_name:
                scores_this_match.append((entry['max_min_rmsd'],design))
    sorted_scores=sorted(scores_this_match, key=lambda first: first[0])
    bests.append(sorted_scores[0])
'''
[(1.0123683214187622,
  'fd2ed1UM_2_M41S95S29W99_1_relaxed_relaxed_5tpj_363097_design_7_unrelaxed_model_3_rank_1_0001_design_2_unrelaxed_model_4_rank_1_0001_design_5_unrelaxed_model_5_rank_1_0001_hybrid_52_1_0001_clean__DE_6_oglig_0001_0008_refined_0001.pdb'),
 (1.1678568124771118,
  'fd5ed1UM_4_M57S53S94W61_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_design_5_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_8_oglig_0001_0010_refined_0001.pdb'),
 (1.1583961248397827,
  'fd6ed1UM_4_M84S99S58W54_1_relaxed_relaxed_5tpj_150791_design_8_unrelaxed_model_1_rank_1_0001_design_1_unrelaxed_model_3_rank_1_0001_design_4_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_15_oglig_0001_0001_refined_0001.pdb'),
 (1.0951082706451416,
  'fd9ed1UM_14_M96S98S57W53_1_relaxed_relaxed_5tpj_140033_design_10_unrelaxed_model_5_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_design_3_unrelaxed_model_2_rank_1_0001_hybrid_52_1_0001_clean__DE_4_oglig_0001_0006_refined_0001.pdb'),
 (1.2135244607925415,
  'fd7ed1UM_2_M54W89W100T38_1_relaxed_relaxed_5tpj_96333_design_10_unrelaxed_model_1_rank_1_0001_design_9_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid_62_1_0001_clean__DE_19_oglig_0001_0003_refined_0001.pdb')]


 '''

filtered_strc=[b for a,b in bests]
os.makedirs('bests',exist_ok=True)
for i in filtered_strc:
    ########################################
    os.system('cp '+i[:-4]+'.pdb bests/'+i+'.pdb')
    ########################################

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl4rm_1_bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas/af2filtered2/filtered ~/desktop/esl_post_order_filtered_designs


OKAY YEAH THIS BATCH DOES ALL LOOK PRETTY DAMN GREAT BY EYE TEST
ALTHOUGH THERE IS A WEIRD SITUATION WITH 3 HBONDS TO A GLUTAMATE CORBXYLATE O
    (i guess possible cus a resonance strc has 3 lone pairs on one O)
AND THEN THERE IS INDEED AN UNPAIRED LIG HYDROXYL IN ONE OF THEM, BUT
THE HBOND NETWORK IS OTHERWISE SO GOOD THAT I THINK IM CONTENT TO ACCEPT IT IN THIS CASE

I think I do actually wanna go ahead and order this 5, and that gives me
some leeway to make one more batch of designs from which to try to order some more
'''











'''


GENE ORDERING


'''




'''
make fastas for designs to order
we dont include methionine if has one at beginning,
because in the constructs (pet28+a) there is a start codon before the n terminal
his tag/thrombin cut site
we also ensure we are not including the ligand code
in addition, we add a * at the end of the sequence to specify a stop codon


/Users/student/desktop/esl_post_order_filtered_designs/bests
'''

import os
l=[i for i in os.listdir() if i[-3:]=='pdb']

d={}
c=1
for pdb in l:
    id='design_'+str(c)
    d[pdb]=id
    c+=1

from pyrosetta import *
init('-ignore_unrecognized_res')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())
    if sequence[0]=='M':
        seq2=sequence[1:]
        sequence=seq2
    if sequence[-1]=='Z' or sequence[-1]=='X':
        seq2=sequence[:-1]
        sequence=seq2
    ofile=open(pdb[:-4]+'.fasta','w')
    ofile.write('>'+d[pdb]+ '\n')
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
    ofile.write(str(last_s) +'*\n')
    ofile.close()

for a in l:
    try:
        fasta(a)
    except:
        print('failure')

os.makedirs('fastas',exist_ok=True)

ll=[i for i in os.listdir() if i[-6:]=='.fasta']
for i in ll:
    os.system('mv '+i+' fastas/'+i)

os.chdir('fastas')
ll=[i for i in os.listdir() if i[-6:]=='.fasta']
of=open('allfastas.fasta','w')
for i in ll:
    f=open(i,'r')
    lines=[line for line in f.readlines()]
    for line in lines:
        of.write(line)
    f.close()
of.close()
'''
now on twist website, we order CLONAL GENES
upload all fasta, codon usage table e.coli
select all, then 'change vector', --> pet28a+ w/ insertion site ndel-xhol
download genbank files of sequences and check a couple with snapgene
    to make sure correct aa sequence
all default parameters for amount of dna n shit
order through B002235895


'''

 '''
 Plasmids (pET-28a(+)) encoding the designed proteins were ordered from Twist
Bioscience. The DNA sequences of the designed proteins were inserted between the NdeI and
XhoI restriction sites, which added the DNA coding sequence for an N-terminal
MGSSHHHHHHSSGLVPRGSHM tag to the designed proteins. Escherichia coli BL21(DE3)
cells were transformed with the plasmid.
'''
