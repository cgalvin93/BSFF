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
