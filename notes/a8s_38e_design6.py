'''
first ima do the 3p1np and 4p dirs, depending on how many matches i get
will determine if its worth screening the 3p2np and 4p1np motifs
    yea actually
'''
import os
#now submit for matching with lucs strc
# lig_name='a8s'
lig_name='38e'
# paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
paramspath='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/a8s_3p1np'
# shellfile_suffix='_a8s_3p1np.sh'
# #
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/a8s_3p2np'
# shellfile_suffix='_a8s_3p2np.sh'
# #
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/a8s_4p'
# shellfile_suffix='_a8s_4p.sh'
# #
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/a8s_4p1np'
# shellfile_suffix='_a8s_4p1np.sh'
#
#
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_3p1np'
# shellfile_suffix='_38e_3p1np.sh'
# #
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_3p2np'
# shellfile_suffix='_38e_3p2np.sh'
# #
# allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p'
# shellfile_suffix='_38e_4p.sh'
# #
allmotifsdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p_1np'
shellfile_suffix='_38e_4p1np.sh'#
#
#
#
#
#
#
#
#
#
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
'''
# import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('_1_a8s_3p1np.sh')
#
for j in jobss:
    cmd='qsub -cwd -t 1-418 -l mem_free=40G -o cluster_out -e cluster_out '+j
    os.system(cmd)



'''
AFTER OVERNIGHT MATCHING THE A8S 3P1NP, I AM DELETING THE JOBS BECAUSE THERE ARE
SO MANY MATCHES, AND the very first has not yet finished

    I WILL ELT THE A8S AND 38E 4P JOBS CONTINUE AS THEY ARE HERE, BUT IM
    GONNA GO AHEAD AND DELETE THE REST , AND RESUBMIT IN A WAY THAT WILL MAYBE BE FASTER?
    I GUESS LESS screens PER job
        and ya know what fuck it i dont think i rly need
        60 gb per job


'''
import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
len(l)
'''
a8s 3p1np
Out[1]: 58645

a8s 4p
Out[1]: 3190
    and this is afetr not running for very long
        Out[1]: 16641
            after leaving it for like another 10 minutes lmao ...
            stopping these jobs here

OK, NOW I WILL SUBMIT 3P2NP AND 4P1NP AND SEE IF THEY CAN GET A GOOD NUMBER OF MATCHES
    submitted for 38e

ITS CLEAR THAT I STRAIGHT UP MADE TOO MANY MOTIFS, I CAN BE PICKIER
AND MAKE A MUCH SMALLER AMOUNT AND SUBMIT FOR SMALLER,FASTER MATCHING JOBS
    im gonna keep the motifs/jobs i submitted, but delete all of them after
    the tenth, thats 100 motifs right

# 3p2np 11 = 449937
# 4p1np 1 = 450060
# 4p1np 10 = 450069
# 4p1np 48 =  450107
for i in range(449937,450060):
    print('qdel -j '+str(i))
for i in range(450069,450107):
    print('qdel -j '+str(i))

                OKAY AND I AM STARTING TO SEE SOME 3P2NP MATCHES COMING IN
scp cgalvin@log2.wynton.ucsf.edu:r6matching/38e_3p2np/UM_1_N38T96K67M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_4_1.pdb ~/desktop/ex1.pdb


getting matches -
design on matches -
starting presentation


HOW DO NATURALLY OCCURRING PROTEINS COMPENSATE THE BUIRIED
UNSAT SC ATOMS IN APO STATE THAT WOULD HBOND WITH LIG IN
HOLO STATE
    i think for a reason like this its worth doing mpnn on designs
    becuase it isnt even considering the ligand so its gonna
    do what it can to stabilize the rest of the protein
        yeah and freeze all the first shell residues


ALRIGHT, NOW I DO WANNA DESIGN ON AT LEAST SOME OF THE
3P1NP AND 4P MATCHES FOR A8S
     and then i can compare design quality with
     those of larger motifs
a8s 3p1np
Out[1]: 58645
a8s 4p
Out[1]: 3190
    i can fd on all the 4p, and take some random selection of the 3p2np

IDK IT REALLY LOOKS LIKE MATCHING JOBS ARENT USING A LOT OF MEMORY NOW,
MAYBE BECAUSE FINDING SO MUCH LESS MATCHES? MAYBE ITD BE OKAY TO SUBMIT WITH
LESS MEMFREE REQUESTED
    i need to submit these fd jobs, but its gonna take a real long time to get to them
    with all these matching jobs in the queue
        NEED TO GO AHEAD AN D DELETE JOBS 5-10 RN TOO TO ADDRESS THIS I THINK,
        IDK MAN HOPEFULLY I get like at least 100 matches or something still

aight i think thats fine, now lemme do some design on these smaller motif a8s matches
OKAY AND COOL FIRST MATCHING JOBS ARE MOSTLY FINISHED AFTER ABOUT 2 HOURS
THATS NOT SO BAD
ESPECIALLY IF I DONT NEED TO ReQUEST SO MUCH MEMORY

'''



'''
a8s3p1np filtdesigncf and filtdesign mpnn running
    while i wait on it i try to get 38e 3p2np and 4p1np to the same point

gonna go ahead and do 38e 3p2np and 4p1np
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
        if len(l)<5:
            bad.append(pdb)
'''
In [2]: len(bad)
Out[2]: 46670
In [4]: len(pdbs2)
Out[4]:

38e3p2np
In [6]: len(pdbs2)
Out[6]: 18

4p1n
In [3]: len(pdbs)
Out[3]: 10
In [4]: len(pdbs2)
Out[4]: 10
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
In [4]: len(firsts)
Out[4]: 361
In [7]: len(firsts)
Out[7]: 2218
    okay thatll do
honestly thats probably like, close to the number of bigger motif matches ill get
eh but i guess thats bad actually, the strength of doing less is getting more matches,
so i should use significantly more, how can i expand on this set? i guess just add randomly
honestly
In [3]: len(firsts)
Out[3]: 651
In [5]: len(firsts)
Out[5]: 2605
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
make resfiles
        well shit, turns out 4p matches are corrupted somehow, the pdb files are
        just empty
DELETING BAD FILES FIRST FOR 3P1NP WORKED
'''
#make a resfile for each match structure
import os
#################################################################################
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'

sfname='r638esrmr.sh'
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
'''
print(len(matches))
2214
qsub -cwd -t 1-2214 -l mem_free=1G -o cluster_out -e cluster_out _1r6a8ssrmr.sh

print(len(matches))
18
qsub -cwd -t 1-18 -l mem_free=1G -o cluster_out -e cluster_out _1r638esrmr.sh

print(len(matches))
10
qsub -cwd -t 1-10 -l mem_free=1G -o cluster_out -e cluster_out _1r638esrmr.sh


'''
# import os
# jobss=[i for i in os.listdir() if i[-3:]=='.sh']
# #
# jobss.remove('_23r6a8ssrmr.sh')
# #
# for j in jobss:
#     cmd='qsub -cwd -t 1-100 -l mem_free=1G -o cluster_out -e cluster_out '+j
#     os.system(cmd)

# changing cst files to match the pdb of matches names
import os
################################################
matches=[i for i in os.listdir() if i[-3:]=='pdb']
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np'
# cstdir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np'
cstdir='/wynton/home/kortemme/cgalvin/r6matching/38e_4p_1np'


################################################
for i in matches:
    ################################################
    cstid='hybrid_'+i.split('_')[-2]+'.cst'
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


#FASTDESIGN ON MATCHES
#setup so that 1 params for all
#resfile and cst must have same name as pdb
import os
########
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'

sfname='r638efd2.sh'
ndesignseachmatch='1'
outputdirectory='38er6fd2'
scorefile_name='38er6fd2.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
for xx in range(100):
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
2214
qsub -cwd -t 1-2214 -l mem_free=16G -o cluster_output -e cluster_output r6a8sfd.sh


qsub -cwd -t 1-18 -l mem_free=16G -o cluster_output -e cluster_output r638efd.sh

'''
import os
jobss=[i for i in os.listdir() if 'r638efd2.sh' in i]
for j in jobss:
    cmd='qsub -cwd -t 1-10 -l mem_free=16G -o cluster_output -e cluster_output '+j
    os.system(cmd)

'''
analyze

a8sr6fd1.json
38er6fd2.json
38er6fd1.json

'''
#analysis of scores from json file
sfname='38er6fd1.json'
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
                allscores.append(float(d[term]))
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

plot_dists(terms,scores,'38e_r6_3p2np.pdf')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
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
    return filtered_scores

f1=return_filtered(scores,'buns2interface','<',3.0)
f2=return_filtered(scores,'contact_molsurf','>',170)
#
f3=return_filtered(scores,'buns2interface','<',1.0)
f4=return_filtered(f3,'lighphobesasa','<',10.0)
f5=return_filtered(f4,'contact_molsurf','>',160)
plot_dists(terms,f5,'38e_r6_3p2np_filt.pdf')

'''
In [23]: len(f5)
Out[23]: 44
'''
import os
filtered_strc=[]
for d in f5:
    filtered_strc.append(d['decoy'])
os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')
os.system('mv *.pdf filtered/')

'''
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered ~/desktop/38er63p2npfiltered
'''


























'''
OKAY NOW NEXT I WANT TO
get colabfold scores of these new sequences.
    and then DO MPNN DESIGN ON THESE WHILE FREEZING
    ALL RES THAT HB WITH LiG OR HAVE OVER A CERTAIN INTERACTION ENERGY
        then af2 those and see if the scoeres improve i guess?
            then rosetta relax with match csts, full score script,
            get filtered of those

I WILL MOVE FORWARD WITH DESIGNTING THE MOTIFS I CURRENTLY HAVE
    but i will discuss in rpesentation whether i shopuld make denser
    hbond motifs and do a deeper study of contact stats for hb fragments
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
NOW RUN MPNN WITH FIRST SHELL RESIDUES FROZEN
maybe the simplest is just to do all residues within a certain distance cutoff
of ligand, that for sure has to include the hbond residues,
hmmm and then maybe check residue energy? shit but then i need to explicitly check hbonds


/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered
/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered


'''
# prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'

ofjsonname='mpnn_params.json'
import os
from pyrosetta import*
import json
init('-load_PDB_components False')

sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(fa_sol, 1)
sf.set_weight(hbond_sc, 1)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,1)
# sf=get_fa_scorefxn()


pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
tfdata={}
for pdb in pdbs:
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
    #######################################################
    ######################################################
    ###ADD MOTIF RESNUMBERS FROM STRING TO FREEZE LIST###############################
    ###THEN I CAN JUST USE NORMAL MATCH CSTS###############################
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    #now analysis of first shell residues
    ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
    neighborhood_selector_bool = neighborhood_selector.apply(p)
    neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    first_shell_res=list(neighborhood_residues_resnums)
    sf(p)
    rosetta.core.pack.optimizeH(p, sf)
    for resnum in first_shell_res:
        # if resnum not in tofreeze:
        #     tofreeze.append((resnum))
        #residue ligand energy
        w=p.energies().energy_graph().find_energy_edge(resnum,p.total_residue())
        w.fill_energy_map()
        hbsc=w[rosetta.core.scoring.hbond_sc]
        hbbbsc=w[rosetta.core.scoring.hbond_bb_sc]
        faatr=w[rosetta.core.scoring.fa_atr]
        farep=w[rosetta.core.scoring.fa_rep]
        fasol=w[rosetta.core.scoring.fa_sol]
        faelec=w[rosetta.core.scoring.fa_elec]
        res_lig_score=faelec+fasol+farep+faatr+hbbbsc+hbsc
        hbtot=hbbbsc+hbsc
        if hbtot<=-0.5:
            if resnum not in tofreeze:
                tofreeze.append((resnum,p.residue(resnum).name()))
        if res_lig_score<=0.:
            if resnum not in tofreeze:
                tofreeze.append(resnum)
    tfdata[pdb]=tofreeze

json.dump(tfdata,open(ofjsonname,'w'))

'''
MPNN WITH CONSTRAINTS
/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json
'''

os.makedirs('mpnn',exist_ok=True)
l=[i for i in os.listdir() if i[-3:]=='pdb']
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
import json
#########
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'
ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn_params.json'
#########
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn'
all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn'
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
import json
import os
#########
# ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn_params.json'

ofjsonname='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn_params.json'
#########
with open(ofjsonname,'r') as f:
    tfdata=json.load(f)
#########
n_strc_per_job=1
# all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn'
all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/mpnn'
directory_prefix='mpnndesign'
###############
shf_prefix='f38eiltdesmpnn'
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
    '--sampling_temp "0.1"',
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
import os
jobss=[i for i in os.listdir() if i[-3:]=='.sh']
#
# jobss.remove('filtdesmpnn_1_.sh')
#
for j in jobss:
    cmd='qsub -cwd -l mem_free=10G -o mpnncout -e mpnncout '+j
    os.system(cmd)














'''
consolidating mpnn results

                MAYBE I SHOULD CONSIDER LOWERING TEMP
                WHEN DESIGNING ON ROSETTA DESIGNS?
'''


import os
#
all_strc_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/mpnn'
directory_prefix='mpnndesign'
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
########################################
selectfastasoutputdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s3p1np_mpnn_best/tocf'
########################################
########################################
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
cd allfastasdir
'''
import os
nfastasperjob=1
directory_prefix='cfmpnnrd1designs'
allfastasdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s3p1np_mpnn_best/tocf'
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
directory_prefix='cfmpnnrd1designs'
allfastasdir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s3p1np_mpnn_best/tocf'
shf_prefix='cf_mpnnrd1best'
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
# jobss.remove('cf_ntf2_r1d_1_.sh')
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


'''
print(len(pdbs))
3
qsub -cwd -t 1-3 -l mem_free=8G score_af2_filt.sh
[(89.23860869565218,
  'UM_25_K88Q97Q99F53_1_relaxed_5tpj_278534_design_2_unrelaxed_model_3_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid_151_1fd1_0001',
  '3',
  0.8045230627059936),
 (93.6191754385965,
  'UM_7_K96N98N51F55_1_relaxed_5tpj_176872_design_8_unrelaxed_model_2_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_hybrid_9_1fd1_0001',
  '4',
  0.7494069337844849),
 (91.46368965517243,
  'UM_9_K98N42N53F57_1_relaxed_5tpj_384973_design_6_unrelaxed_model_2_rank_1_0001_design_7_unrelaxed_model_2_rank_1_0001_hybrid_9_1fd1_0001',
  '5',
  0.863116979598999)]


[(91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0019',
  '5',
  0.8312956929206848),
 (91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0021',
  '5',
  0.8688668847084046),
 (91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0042',
  '5',
  0.8868943095207215),
 (91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0061',
  '5',
  0.8446571588516235),
 (91.65067857142856,
  'UM_2_N38T96K67M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_28_1fd1_0081',
  '4',
  0.8520600438117981),
 (91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0079',
  '5',
  0.8209141969680787),
 (91.6506607142857,
  'UM_2_N38T96K67M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_28_1fd1_0084',
  '4',
  0.8506036877632142),
 (91.693875,
  'UM_1_N38T96N85M94I17_1_relaxed_5tpj_258209_design_10_unrelaxed_model_3_rank_1_0001_design_10_unrelaxed_model_1_rank_1_0001_hybrid_32_1fd1_0086',
  '5',
  0.7954639554023742)]
'''
#move these to a new folder where im gonna try to relax them
ns=[]
for i in aplddts:
    n='_'.join([i[1],'unrelaxed_model',i[2]])
    ns.append(n)
pn=[]
for n in ns:
    for path in prediction_pdb_paths:
        s=path.split('/')[-1]
        ss=s.split('_')[:-2]
        sss='_'.join(ss)
        if n==sss:
            pn.append(path)

os.makedirs('high_conf',exist_ok=True)
for n in pn:
    try:
        os.system('cp '+n+' high_conf/'+n.split('/')[-1])
    except:
        print(n)

#now try to relax and get ligand back into these
designs='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas/high_conf'
predesign_dir='/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design'
ligand_pdb='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb'
import os
l=[i for i in os.listdir(designs) if i[-3:]=='pdb']
ns=[]
for i in l:
    s=i.split('_')[:26]
    ss='_'.join(s)
    sss=ss[:-3]
    ns.append((sss+'.pdb',i))
l2=[i for i in os.listdir(predesign_dir) if i[-3:]=='pdb']
os.makedirs('torelax')
for n in ns:
    templatelines=[]
    liglines=[]
    for i in l2:
        if n[0]==i:
            ofp=os.path.join(predesign_dir,i)
            f=open(ofp,'r')
            lines=[line for line in f.readlines() if line[:6]=='REMARK' or line[:6]=='HETATM']
            f.close()
            for line in lines:
                if line[:6]=='REMARK':
                    templatelines.append(line)
                elif line[:6]=='HETATM':
                    liglines.append(line)
    f=open(os.path.join(designs,n[1]),'r')
    alines=[line for line in f.readlines() if line[:4]=='ATOM']
    f.close()
    ofn=os.path.join('torelax',n[1])
    of=open(ofn,'w')
    for li in templatelines:
        of.write(li)
    for li in alines:
        of.write(li)
    for li in liglines:
        of.write(li)
    of.close()
l3=[i for i in os.listdir(predesign_dir) if i[-3:]=='cst']
for n in ns:
    for i in l3:
        if n[0].strip('.pdb')==i.strip('.cst'):
            ofp=os.path.join(predesign_dir,i)
            os.system('cp '+ofp+' torelax/'+n[1].strip('.pdb')+'.cst')

#now submit relax
import os
########
prm1='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
# prm1='/wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e.params'
sfname='reldes.sh'
ndesignseachmatch='1'
outputdirectory='relaxed'
scorefile_name='a8s3p1npaf2desrelaxed.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
c=1
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
     '-parser:protocol','~/BSFF/tools/fr_addligback.xml',
     '-parser:view','-run:preserve_header',
     '-in:file::extra_res_fa',prm1,
     '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
     '-ex1','-ex2','-extrachi_cutoff','0',
     '-out:nstruct','1',
     '-score::weights','ligand',
     '-load_PDB_components','False',
     '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
     '-relax:constrain_relax_to_start_coords', '-relax:coord_cst_stdev .5',
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

qsub -cwd -t 1-12 -l mem_free=4G -o cluster_output -e cluster_output _1reldes.sh
scp -r cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas/high_conf/torelax/relaxed ~/desktop/relex




okay relax is not working but i will run ligbinderanalysis on the 3 designs
which passed af filter
then do the same for 38e

import json
rjson='a8s3p1np_filt_af2.json'
df=open(rjson,'r')
data_allstrc=json.load(df)
df.close()
'''
#plot avg plddt (across 5 models) for each structure
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
##################################################################
outfilename='38e3p2np_filt_af2.pdf'
##################################################################
pdf = PdfPages(outfilename)
#get list of avg plddt values for each strc to plot
aplddts=[]
acarmsds=[]
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        plddts.append(cplddt)
        carmsds.append(data_allstrc[key][k2][2])
    aplddt=np.mean(plddts)
    aplddts.append(aplddt)
    acarmsd=np.mean(carmsds)
    acarmsds.append(acarmsd)
#plot them
fig,ax=plt.subplots()
ax.hist(aplddts)#bins=int(len(allscores)/20)
ax.set_title('Average plddt for each Sequence across 5 Models')
ax.set_ylabel('Frequency')
ax.set_xlabel('Plddt')
pdf.savefig()
#
fig,ax=plt.subplots()
ax.scatter(acarmsds,aplddts)#bins=int(len(allscores)/20)
ax.set_title('CaRMSD vs Plddt')
ax.set_ylabel('Plddt')
ax.set_xlabel('CA RMSD')
pdf.savefig()
plt.clf()
#
pdf.close()
'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas/a8s3p1np_filt_af2.pdf ~/desktop/a8s3p1np_filt_af2.pdf

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/cf/fastas/38e3p2np_filt_af2.pdf ~/desktop/38e3p2np_filt_af2.pdf



/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/cf/fastas/af2filtered


/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas/af2filtered

import json
with open('38e_polar_contacts_clustered_liberal.json','r') as f:
    df=json.load(f)


ANALYZE SCORES OF RELAXED CF MODELS OF FILT DESIGNS
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
                allscores.append(float(d[term]))
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

plot_dists(terms,scores,'38e_r6_3p2np_relcffilt.pdf')

'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/a8s_3p1np/design/a8sr6fd1/filtered/cf/fastas/af2filtered/a8s_r6_3p2np_relcffilt.pdf ~/desktop/a8s_r6_3p2np_relcffilt.pdf

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/r6matching/done/38e_3p2np/design/38er6fd1/filtered/cf/fastas/af2filtered/38e_r6_3p2np_relcffilt.pdf ~/desktop/38e_r6_3p2np_relcffilt.pdf

INSTEAD OF ENUMERATING MOTIFS :
    find overrep for each frag (top nth percentile pop)
    randomly draw from each and then apply filters...




'''
