'''
OKAY RATHER THAN JUST USING RNG HERE LEMME JUST MAKE A SET NUMBER OF PROBLEMATIC
AND >3 RES MOTIFS....
    lets see, and ima add np contacts to some subset of them, so how many do i
    wanna make fr
        i make some number of 3p,
            from which some number 3p1np, 3p2np,
        i make some number of 4p
            from which some number of 4p1np
    matching small number 3p1np
    very small number 4p
    bigger number 4p1np
    about the same bigger number 3p2np

i only wanna match maybe 150 motifs total for each molecule
    okya no i did 100 for 3p2np/4p1np b4 and let it run thru and it wasnt that bad,
    i just want a small amount of smaller motifs cus they make so many matches
    and i dont wanna have to kill the jobs
    100 3p2np
    100 4p1np
    10 4p
    10 3p1np

AND THEN I WILL EXPERIMENT WITH DIFFERENT DESIGN METHODS TO
PROMOTE BUILDING HBOND NETWORKS AROUND BINDING SITE
    3bop
    interface/hbond upweight
    cm
ALTERNATIVELY MAYBE THERE IS SOME WAY TO PROMOTE INTER RESIDUE HBONDING IN MY
MOTIFS UP FRONT - eg modeling them explicitly again

make 100 3rp
100 4rp
    to each of these i will add the np contacts to get those

OKAY THEN, I WILL MAKE THESE NUMBERED MOTIFS SEPARATELY THEN AND FOR EACH ALLOW
20% PROBLEMATIC MOTIFS



'''

















'''
first put the desired motifs all in the same directory, i will design and
analyze the matches as a batch
                                100 3p2np
                                100 4p1np
                                10 4p
                                10 3p1np
'''
###########################################################################
###########################################################################
###########################################################################
import os
# alldir='/wynton/home/kortemme/cgalvin/round6/38e/rand_motifs'
# motifsdir='/wynton/home/kortemme/cgalvin/round6/38e/38e_polar_csts_rand'
alldir='/wynton/home/kortemme/cgalvin/round6/a8s/rand_motifs'
motifsdir='/wynton/home/kortemme/cgalvin/round6/a8s/a8s_polar_csts_rand'
###########################################################################
###########################################################################
###########################################################################
#
os.makedirs(alldir,exist_ok=True)
#
m_3p2np=os.path.join(motifsdir,'three_res','2np')
m_4p1np=os.path.join(motifsdir,'four_res','1np')
#
m_3p1np=os.path.join(motifsdir,'three_res','1np')
m_4p=os.path.join(motifsdir,'four_res')
#
for i in os.listdir(m_3p2np):
    if i[-3:]=='cst':
        f=open(os.path.join(m_3p2np,i),'r')
        lines=[line for line in f.readlines()]
        f.close()
        if len(lines)>20:
            os.system('cp '+os.path.join(m_3p2np,i)+' '+alldir+'/m1_'+i)
        else:
            print(i)
#
for i in os.listdir(m_4p1np):
    if i[-3:]=='cst':
        f=open(os.path.join(m_4p1np,i),'r')
        lines=[line for line in f.readlines()]
        f.close()
        if len(lines)>20:
            os.system('cp '+os.path.join(m_4p1np,i)+' '+alldir+'/m2_'+i)
        else:
            print(i)
#
for i in os.listdir(m_3p1np)[:10]:
    if i[-3:]=='cst':
        f=open(os.path.join(m_3p1np,i),'r')
        lines=[line for line in f.readlines()]
        f.close()
        if len(lines)>20:
            os.system('cp '+os.path.join(m_3p1np,i)+' '+alldir+'/m3_'+i)
        else:
            print(i)
#
for i in os.listdir(m_4p)[:10]:
    if i[-3:]=='cst':
        f=open(os.path.join(m_4p,i),'r')
        lines=[line for line in f.readlines()]
        f.close()
        if len(lines)>20:
            os.system('cp '+os.path.join(m_4p,i)+' '+alldir+'/m4_'+i)
        else:
            print(i)

#IT ACTUALLY MIGHT BE SMART TO RENAME ALL OF THESE MOTIFS
#SO THAT THEY HAVE SAME NUMBER OF UNDERSCORES, THAT WAY IF I AM SPLITTING
#NAMES LATER THE INDICES WILL BE THE SAME FOR ALL OF THEM
import os
l=[i for i in os.listdir() if i[-3:]=='cst']
c=1
for i in l:
    os.system('mv '+i+' m_'+str(c)+'.cst')
    c+=1





'''
                            MATCHING
'''
###########################################################################
###########################################################################
###########################################################################
###########################################################################
import os
#now submit for matching with lucs strc
#
lig_name='a8s'
# lig_name='38e'
#
paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# paramspath='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
#
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p_1np'
allmotifsdir='/wynton/home/kortemme/cgalvin/round6/a8s/new_t3_polar_csts_rand/three_res/1np'
#
shellfile_suffix='_a8s_rpm_match.sh'#
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
                           '-match::lig_name', lig_name]
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
404
n scaffolds
418

qsub -cwd -t 1-418 -l mem_free=60G -o cluster_out -e cluster_out _1_a8s_3p1np.sh

qsub -cwd -t 1-418 -l mem_free=60G -o cluster_out -e cluster_out _1_a8s_3p1np.sh

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/round6/a8s/rand_motifs/UM_2_T69S17T41N54F13_1_relaxed_5tpj_403371_design_6_unrelaxed_model_2_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_m_108_1.pdb ~/desktop/ex.pdb

WEIRD UNSAT THR HBONDS N SHIT
    wouldnt be too terrible of a problem if i had a lotta matches, maybe a good number would
    be good, but i only got a very small handful of matches and now they all suck
    and are worthless, huge L
okay ya know what, im gonna try to model hbonds explicitly again
    spatially cluster them as finely as i can
        2body energies all terms ptn lig,
        ptn ptn only faatr/rep and hbscsc(or some other way to reward network density?)
AND SOME FUCKING BULLSHIT 1 RESIDUE MOTIF SNUCK IN HERE SO NOW ITS GIVING A MILLION
MATCHES, SO ANNOYING, WHAT THE FUCK EVEN HAPPENED HERE HOW ON EARTH DID I CREATE
THIS THING UUURRRGHGHHHGHHHHHH
import os
l=[i for i in os.listdir()]
for i in l:
    if '_m_200_' in i:
        os.system('rm '+i)

With respect to hydrogen bonds, the protein environment can better position hydrogen bond donors and acceptors, so that less rearrangement is needed relative to hydrogen bonds in solution, with positioning “paid for” by protein folding energy and ligand binding energy.
Additionally, there may be more hydrogen bonds formed within an active site than on average in water, and the active site hydrogen bond donors and acceptors can have charge densities higher than those of water and, as a result, undergo a greater stabilization going from the ground to transition state


OKAY NOW I AM THINKING RATHER THAN WANTING MORE PTN-LIG HB,
    i just wanna make sure the ones im matchinga re good, and will be counted
    some kinda analysis and filtering of hbond statistics
AND I WANNA PROMOTE SECONDARY HB AMONGST FIRST SHELL RES THAT HB WITH LIG, THIS SEEMS
LIKE A REALLY GOOD LOOK DAWG
    maybe just do this in design stage or something though, i dont think i can
    really promote p front without modeling motif sidechains explicitly
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
PROBABLY WHAT I WANT TO DO IS NOT USE MIDDLE BIN VALUES, BUT ACTUAL IDEAL
VALUE OF DIST AND ANGLE FOR CSTS
    i shouldnt use match csts, just set matched residues to repack
    when i design



    IN GOOD/OTHER NAT BINDERS
        avg number of hbonds for each kind of polar atom
        res preferences for each kind of polar atom
            res pair pref for double hbonded
        geometric preferences stats for each res-ligatom
        nbb vs nsc - lig
        network density
            but breaking down ss hbs vs extra hbs

i only get 8 matches, they dont look great
im gonna try explicitly modeling and see what i can get
to do this, i need the input fuzzball and frequency json file

i is resnumber
frequency=frequency_dict[str(i)][0]
ref=frequency_dict[str(i)][1]

OKAY FIRST I CAN JUST USE CLUSTNBP SCRIPT ON POLAR FRAGS AND ITLL DO THE SPATIAL
CLUSTERING FOR ME
    i make a new version of it that will have proper 'polar residues' list
    spatial_cluster_polar.py
time python3 ~/BSFF/spatial_cluster_polar.py a8s /wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb Fragment_1 Fragment_2 Fragment_3

compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/spatial_cluster_polar.py a8s /wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb Fragment_1 Fragment_2 Fragment_3
qstat -j "$JOB_ID"
'>cnp_$compound_dir.sh
chmod ugo+x cnp_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output cnp_$compound_dir.sh

now i make a script based on fuzzball_for_assembly.py
called ffa2.py which will
produce the fuzzball and the json file that i want to send to mc assembly



STILL SPATIALLY CLUSTERING THE POLAR RES, BUT MAYBE THERE IS AN EASY WAY TO DO?

            i will simply try to get a lotta matches w/ 4p motifs, see
            if i can then
'''


import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
len(l)

'''
Out[1]: 206

Out[1]: 5
'''

#clean matches and put in own directory
import os
os.makedirs('design',exist_ok=True)
pdbs=[i for i in os.listdir() if i[:2]=='UM']
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
'''
In [3]: len(bad)
Out[3]: 0
In [3]: len(bad)
Out[3]: 0

'''
pdbs2=[i for i in pdbs if i not in bad]
firsts=pdbs2
#17033
#58645
# scaffs=[]
# motifs=[]
# firsts=[]
# for pdb in pdbs2:
#     s=pdb.split('_')
#     seq=s[2]
#     scaff=s[4:-2]
#     motif=s[-2]
#     if scaff not in scaffs:
#         if pdb not in firsts:
#             firsts.append(pdb)
#             scaffs.append(scaff)
#     if motif not in motifs:
#         if pdb not in firsts:
#             firsts.append(pdb)
#             motifs.append(motif)
# import random
# for i in range(2000):
#     pdb=random.choice(pdbs2)
#     if pdb not in firsts:
#         firsts.append(pdb)
'''
'''
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
#
# for pdb in firsts:
#     f=open(pdb,'r')
#     l=[li for li in f.readlines()]
#     cl=[]
#     for i in l:
#         if i[:4]=='ATOM':
#             cl.append(i)
#         if i[:6]=='HETATM':
#             cl.append(i)
#         if i[:10]=='REMARK 666':
#             cl.append(i)
#     ofp=os.path.join('design',pdb.split('/')[-1])
#     of=open(ofp,'w')
#     for line2 in cl:
#         of.write(line2)
#     of.close()
'''
consolidated all serthr to 1 dir
make resfiles
'''
#make a resfile for each match structure
import os
#################################################################################
prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'

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
219
qsub -cwd -t 1-219 -l mem_free=1G -o cluster_out -e cluster_out _1srmr.sh

'''
# import os
# jobss=[i for i in os.listdir() if i[-3:]=='.sh']
# #
# jobss.remove('_23r6a8ssrmr.sh')
# #
# for j in jobss:
#     cmd='qsub -cwd -t 1-100 -l mem_free=1G -o cluster_out -e cluster_out '+j
#     os.system(cmd)

# # changing cst files to match the pdb of matches names
# import os
# ################################################
# matches=[i for i in os.listdir() if i[-3:]=='pdb']
# # cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np'
# # cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np'
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p_1np'
#
#
# ################################################
# for i in matches:
#     ################################################
#     cstid='hybrid_'+i.split('_')[-2]+'.cst'
#     ################################################
#     ogcst=open(os.path.join(cstdir,cstid),'r')
#     ################################################
#     lines=[line for line in ogcst.readlines()]
#     ogcst.close()
#     newcstname=i.split('.')[0]+'.cst'
#     of=open(newcstname,'w')
#     for line in lines:
#         of.write(line)
#     of.close()
#FASTDESIGN ON MATCHES
#setup so that 1 params for all
#resfile and cst must have same name as pdb
# import os
# ########
# # prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
#
# sfname='r638efd2.sh'
# ndesignseachmatch='1'
# outputdirectory='38er6fd2'
# scorefile_name='38er6fd2.json'
# #########
# matches=[i for i in os.listdir() if i[-3:]=='pdb']
# os.makedirs('cluster_output',exist_ok=True)
# os.makedirs(outputdirectory,exist_ok=True)
# c=1
# for xx in range(100):
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
#          '-parser:protocol','~/BSFF/tools/FD_jump1_nofilters.xml',
#          '-parser:view','-run:preserve_header',
#          '-in:file::extra_res_fa',prm1,
#          '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
#          '-ex1','-ex2','-extrachi_cutoff','0',
#          '-out:nstruct',ndesignseachmatch,
#          '-score::weights','ligand',
#          # '-score:set_weights','hbond_bb_sc','2.0',
#          # '-score:set_weights','hbond_sc','2.0',
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
# 2214
# qsub -cwd -t 1-2214 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh
#
#
# qsub -cwd -t 1-18 -l mem_free=16G -o cluster_output -e cluster_output r638efd.sh
#
# '''
# import os
# jobss=[i for i in os.listdir() if 'r638efd2.sh' in i]
# for j in jobss:
#     cmd='qsub -cwd -t 1-10 -l mem_free=16G -o cluster_output -e cluster_output '+j
#     os.system(cmd)
'''
fastdesign without match csts

scp FD_jump1_nomatchcsts.xml cgalvin@log2.wynton.ucsf.edu:BSFF/tools/
'''
#FASTDESIGN ON MATCHES
#setup so that 1 params for all
#resfile and cst must have same name as pdb
import os
########
prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'

sfname='fd.sh'
ndesignseachmatch='10'
outputdirectory='a8s_st_fd'
scorefile_name='a8s_st_fd.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(5):
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
         '-parser:protocol','~/BSFF/tools/FD_jump1_nomatchcsts.xml',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-out:nstruct',ndesignseachmatch,
         '-score::weights','ligand',
         # '-score:set_weights','hbond_bb_sc','2.0',
         '-score:set_weights','hbond_sc','2.0',
         '-load_PDB_components','False',
         # '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
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
219
'''
import os
jobss=[i for i in os.listdir() if 'fd.sh' in i]
for j in jobss:
    cmd='qsub -cwd -t 1-219 -l mem_free=16G -o cluster_output -e cluster_output '+j
    os.system(cmd)























#analysis of scores from json file
sfname='a8s_st_fd.json'
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

plot_dists(terms,scores,'a8s_st.pdf')

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
f4=return_filtered(f3,'contact_molsurf','>',150)
plot_dists(terms,f4,'a8s_st_filt.pdf')

print(len(f1))
print(len(f2))
print(len(f3))
print(len(f4))
'''
2169
346
141
133
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
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/a8sserthr/a8s_st_fd/filtered ~/desktop/a8s_st_filt
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
directory_prefix='filtdesigns'
allfastasdir=os.getcwd()
shf_prefix='_38eogdcf'
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
allfastasdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/cf/fastas'
# allfastasdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas'
# nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered'
nonaf2desdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered'

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
    currstrc_name='_'.join(p.split('/')[-1].split('_')[:-5])
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
    currstrc_model=str(p.split('/')[-1].split('_')[-3])
    ##################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#output the data so its easy to load later and compare with designs n shit
##################################################################
##################################################################
##################################################################
# json.dump(data_allstrc,open('a8s3p1np_filt_af2.json','w'))
json.dump(data_allstrc,open('38e3p2np_filt_af2.json','w'))
##################################################################
##################################################################
##################################################################

import numpy as np
#id designs with plddt over threshold
plddt_threshold=85.0
carmsd_threshold=1.0
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

#move filtered to own directory to run scoring script on
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

#run the scoring in af2filtered dir
os.chdir('af2filtered')
# paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
paramspath='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
import os
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
sf=open('score_af2_filt.sh','w')
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
os.system('qsub -cwd -t 1-8 -l mem_free=8G score_af2_filt.sh')
os.chdir('..')
#
