'''
get a8s clustering, gen5 motifs matching rn (get those and design on them)
then do new moad stuff



'''
'''
Create a directory for the new target
'''
time python ~/desktop/BSFF/new_target.py a8s
time python ~/desktop/BSFF/new_target.py 38e

'''
download a sdf/molfile of the target chemical structure (RCSB or Pubchem)
open in avogadro
    build --> add hydrogens for ph --> 7.4
    autooptimization tool --> ff=mmff94 --> start (run until dG=0)
    save as --> target_opt.mol

NOW I AM ADDING A HYDROGEN TO THE CARBOXYLATE FOR THE OPTIMIZED DEPROTONATED CONFORMER
    load optimized w/ right protonation
    use selection tool to select whole molecule
        ->extensions->mol mech->fix selected atoms
    now add the desired hydrogen and apply geo opt again :)
'''

'''
create params for both H and proper protonated version
compounds=("a8s")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}_opt.mol
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n ${I}H -p ${I}H ${I}_optH.mol
done

compounds=("38e")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}_opt.mol
done


move params and pdb to Inputs/Rosetta_Inputs:
compounds=("a8s")
for I in ${compounds[@]}; do
    mv ${I}.params ${I}/Inputs/Rosetta_Inputs/${I}.params
    mv ${I}_0001.pdb ${I}/Inputs/Rosetta_Inputs/${I}_0001.pdb
    mv ${I}H.params ${I}/Inputs/Rosetta_Inputs/${I}H.params
    mv ${I}H_0001.pdb ${I}/Inputs/Rosetta_Inputs/${I}H_0001.pdb
done

compounds=("38e")
for I in ${compounds[@]}; do
    mv ${I}.params ${I}/Inputs/Rosetta_Inputs/${I}.params
    mv ${I}_0001.pdb ${I}/Inputs/Rosetta_Inputs/${I}_0001.pdb
done

'''

'''
fragment search
with frags.txt files in target parent directories
DELETED ALKENE CHAIN FRAGMENT
MADE DIMETHYL SUBSTIT OFF RING MORE MINIMAL

        A8S
    1-carboxylate
    2-alcohol
    3-ketone
    4-cyclohexene
    5-ethyl substituent on ring (no filtered contacts)
        DFHBI
    1-carbonyl
    2-secondary imidazole N
    3-alcohol
    4-benzene
    5-cf1
    6-cf2


NOW I WILL DO FRAG SEARCH WITH PROTONATED VERSION OF LIG PDB FILE
    changed to hopefully not add redundant ligs or pdbs to the lists
compounds=("a8s")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}H_0001.pdb
    cd ..
done

compounds=("38e")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}_frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
    cd ..
done

~10 minutes

'''


'''
NOW setup to split into multiple jobs for alignment on cluster
conda activate pyr37
compounds=("a8s")
for I in ${compounds[@]}; do
    cd ${I}
    mkdir cluster_output
    mv *.pdb Inputs/Fragment_Inputs/
    mv *.mol Inputs/Fragment_Inputs/
    time python ~/desktop/BSFF/process_fragment_search.py
    cd ..
done

compounds=("38e")
for I in ${compounds[@]}; do
    cd ${I}
    mkdir cluster_output
    mv *.pdb Inputs/Fragment_Inputs/
    mv *.mol Inputs/Fragment_Inputs/
    time python ~/desktop/BSFF/process_fragment_search.py
    cd ..
done

okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r round6 cgalvin@log2.wynton.ucsf.edu:
scp -r BSFF cgalvin@log2.wynton.ucsf.edu:
scp -r 38e cgalvin@log2.wynton.ucsf.edu:round6/


submit the alignment jobs
a8s
qsub -cwd -t 1-2292 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh
qsub -cwd -t 1-1980 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh


'''





'''

SET UP SHELLFILE FOR SUBMITTING FILTER JOB
in target parent directory:
time python ~/BSFF/setup_filter_contacts.py

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/setup_filter_contacts.py
'>run_setupfilter.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_setupfilter.sh


a8s
1082

38e
1355

chmod ugo+x run_filter_array.sh
qsub -cwd -t 1-1082 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh
qsub -cwd -t 1-1355 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            changed fillter contacts array .py
            if float(occupancy)<0.7 or float(temp_factor)>40.0:
                from >0 and >60
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
doing the actual filtering is quite fast, only the setup takes many hours




CONSOLIDATE FILTERED CONATCTS
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/consolidate_filtered_contacts.py
'>run_consoldate.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_consoldate.sh




EVIDENTLY MY MATCHING JOBS WERE USING UPWARDS OF 30GB OF MEMORY ,
ALSO THEY WANT EVERYTHING AS ARRAYS RATHER THAN IND JOBS
AND ALSO THEY SAID ITS BAD IF MANY JOBS DRAWING DATA FROM SAME DIR?
SO MAYBE SPREAD THINGS OUT BETWEEN MULTIPLE DIRECTORIES



NOW OUTSIDE OF PARENT DIRECTORY
mkdir cluster_output
yoooo i should subm,it each fragment for clustering on a different node
that should be easu enough
a8s frag 5 had nothing

compound_dir='a8s'
compound_dir='38e'


            DELETE A8S PDB FILE WITH H IN IT
            FULL LIGAND PDB IN CLUSTERING IS PULLED FROM INPUTS/ROSETTAINPUTS
            JUST FINDS PDB FILES IN THERE, takes first one
            SHOULD REALLY ONLY BE ONE IN THERE THOUGH, IE THE ONE WE ACTUALLY WANT
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_1
qstat -j "$JOB_ID"
'>clust1_$compound_dir.sh

chmod ugo+x clust1_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust1_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_2
qstat -j "$JOB_ID"
'>clust2_$compound_dir.sh

chmod ugo+x clust2_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust2_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_3
qstat -j "$JOB_ID"
'>clust3_$compound_dir.sh

chmod ugo+x clust3_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust3_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_4
qstat -j "$JOB_ID"
'>clust4_$compound_dir.sh

chmod ugo+x clust4_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust4_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_5
qstat -j "$JOB_ID"
'>clust5_$compound_dir.sh

chmod ugo+x clust5_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust5_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_fragment.py '$compound_dir' Fragment_6
qstat -j "$JOB_ID"
'>clust6_$compound_dir.sh

chmod ugo+x clust6_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust6_$compound_dir.sh


            OKAY, LETS SEE HOW LONG THIS TAKES WITH FRAG CLUSTERING SEP
            ON DIFF NODE EACH
            set at 12:30 pm wed june 22


#analysis of scores from json file
sfname='a8sr6fd1.json'


import json
scores2 = [json.loads(line) for line in open(sfname,'r')]
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
    pdf.close()      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  872.43 GiB| 1000.00 GiB|| 10485998|unlimited

match analysis (stopped early but whatever)
    and design on them if id like
develop hbond filtering to favorable params only (that one paper)
    will want to apply this to new a8s/38e motif generation
clust comparison a8s round5/6
    comparison between eq frgas too, like a8s5 benzene vs 3ng5 benzene
    then i will gen motifs from r6 clusts and do new matching in a way that makes wynton happy :)
moad


'''




































'''
OKAY, DESIGN ON THE MATCHES, LETS DO SHORT RANGE FIRST
'''

#clean matches and put in own directory
import os
os.makedirs('design',exist_ok=True)
pdbs=[i for i in os.listdir() if i[:2]=='UM']
for pdb in pdbs:
    s='egrep "^ATOM|HETATM|REMARK 666" '+pdb+' >design/'+pdb
    print(s)
    os.system(s)

'''
gonna change the resfile script to find any designable
then find any repackable by finding atoms close to designable, instead of a cut
distance from the ligand
'''
#make a resfile for each match structure
import os
prm1='/wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params'
matches=[i for i in os.listdir() if i[-3:]=='pdb']
sfname='r5a8ssrmr.sh'
os.makedirs('cluster_out',exist_ok=True)
sf=open(sfname,'w')
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
print(len(matches))
502
qsub -cwd -t 1-502 -l mem_free=1G -o cluster_output -e cluster_output r5a8ssrmr.sh





'''
# changing cst files to match the pdb of matches names
import os
matches=[i for i in os.listdir() if i[-3:]=='pdb']
cstdir='/wynton/home/kortemme/cgalvin/a8sr5matches/a8spmlt3_polar_csts'
for i in matches:
    cstid=i.split('_')[6]+'.cst'
    ogcst=open(os.path.join(cstdir,cstid),'r')
    lines=[line for line in ogcst.readlines()]
    ogcst.close()
    newcstname=i.split('.')[0]+'.cst'
    of=open(newcstname,'w')
    for line in lines:
        of.write(line)
    of.close()
'''

#FASTDESIGN ON MATCHES
#setup so that 1 params for all
#resfile and cst must have same name as pdb
import os
########
prm1='/wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params'
sfname='r5a8sfd.sh'
ndesignseachmatch='1'
outputdirectory='a8sr5fd1'
output_suffix='fd1'
scorefile_name='a8sr5shorthbfd1.json'
#########
matches=[i for i in os.listdir() if i[-3:]=='pdb']
os.makedirs('cluster_output',exist_ok=True)
os.makedirs(outputdirectory,exist_ok=True)
sf=open(sfname,'w')
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
     '-out:suffix',output_suffix,
     '-scorefile_format','json',
     '-out:file:scorefile',scorefile_name]
sf.write((' ').join(cmd))
sf.write('\nqstat -j "$JOB_ID"')
sf.close()
print(len(matches))
502
qsub -cwd -t 1-502 -l mem_free=16G -o cluster_output -e cluster_output r5a8sfd.sh

check on mem reqs after finished

'''
okay i am realizing that buns2interface filter is brute
it counts every lone pair as a potential hb
so carboxylate needs 4 hb to be considered as having 0, or 2 otherwise
have some exposed to water
        so theres a question of either trying to design this many hbonds into
        the motifs
        or looking at how often they have all satisfied in moad binders and
        getting a sense of how often i leave unsat

I THINK A GOOD PLAN IS
            1 DO THE MOAD STUFF
            2 WHILE I WAIT FOR NEW CLUSTERS TO FINISH, DEVELOP A PLAN
              FOR THE NEXT GENERATION OF MOTIFS AND DESIGN AND GET THE CODE
              READY FOR WHEN CLUSTERING IS FINISHED

I COULD GREATLY INCREASE MATCH EFFICIENCY IF I FIND FULL FAVORABLE RANGES
FOR EACH RES3ATOM-LIG3ATOM COMBO and specify those in the cst blocks
instead of separating into multiple clusters with only slightly different bin values

filter by distance
    filter to 'statistically overrep' clusters only
        identify all unique res3atom-lig3atom pairs
        for each:
            find distribution for hbond params
            identify ranges for each param independtly which are considered stat overrep relative to random dist of points
            use these ranges to define single cst block for that res3atom-lig3atom pair, store these

why not just do this upfront using all of the polar contact data?
    i could, and this actually makes more sense if this is how im going at it
    then could compare these distributions with those for ptn-ptn interactions in rosetta - are they the same?
I THINK MAYBE THE SUPER LONG TAIL ON MY CLUSTER DISTRIBUTIONS
MEANS THAT THE OVERREPRESENTED CLUSTERS ARE STILL NOT ACCOUNTING FOR A LARGE SHARE
OF THE TOTAL CONTACTS
    yea okay ya know what this is convincing me its worth the work to try this
        YEA FIND PARAM DISTRIBUTIONS FOR EACH INDIVIDUAL ATOMTYPE PAIRING

THE WAY IM DOING IT NOW HAS ADVANTYAGE OF NOT TREATING PARAMS INDEPENTLY
    however i am atm penalizing residues that have more than one hbonding atom,
        i should fix that, could for instance multiply arg cluster sizes by 2
        2 equivalent atoms can make (or is it 3? idk exactly how  to treat this)

OKAY I CHANGE MY MIND, IM GONNA MAKE THAT EASY CHANGE BUT JUST WAIT TIL I HAVE CLUST RSULTS AND
TRY THAT FIRST, ALSO I SHOULD MERGE CLUSTERS WHICH DIFFER BY ONLY 1 BIN VALUE FROM TOP
CLUSTER FOR A GIVEN RES3ATOM-LIG3ATOM PAIR (this seems like a good compromise)
    maybe like, for a given bin, if i can find another clust with value +-1, i will
    gen the cst as the whole range

THEN IM JUST GONNA FOCUS ON BINDING METRICS TO FILTER, THEN GET ALL METRICS,
THEN GET AF2 METRICS, THEN SUBMIT TO MPNN, THEN SEE HOW THESE METRICS IMPROVE


3 POLAR 1 NP
3 POLAR 2 NP
4 POLAR
4 POLAR 1 NP



'''
























'''


MATCH ANALYSIS
mkdir a8sr5matches
mv a8slonghb a8sr5matches
mv a8sshorthb a8sr5matches
mv a8sr5matches ~/

import os
matches=[i for i in os.listdir() if i[:2]=='UM']
len(matches)
# Out[1]: 502 short
#n matches for each scaffold type
    #n double matches (um2/3/4/etc...)
    #n unqiue motifs
    #n unique scaffolds
    #n unique res-seqs
scaffs=[]
c1=0
c2=0
c3=0
seqs=[]
for match in matches:
    scaff=('_').join(match.split('_')[4:6])
    motif=match.split('_')[6]
    mn=match.split('_')[1]
    seq=match.split('_')[2]
    l=[]
    for s in seq:
        if not s.isdigit():
            l.append(s)
    seqs.append(('').join(l))
    scaffs.append((scaff,motif))
    if int(mn)==1:
        c1+=1
    elif int(mn)==2:
        c2+=1
    else:
        c3+=1
print(str(c1))
print(str(c2))
print(str(c3))
print(str(len(set(seqs))))
print(set(seqs))
# 486
# 16
# 0
#12 really 11!
# {'WSI', 'IGW', 'TSS', 'WGI', 'IGF', 'IGW', 'WSF', 'WGF', 'FSF', 'IWI', 'WSS', 'IWF'}
us=[]
um=[]
for a,b in scaffs:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 178
# 73

ross=[]
ntf2=[]
nat=[]
for match in matches:
    scaff=('_').join(match.split('_')[4:6])
    motif=match.split('_')[6]
    if scaff.split('_')[0]=='model':
        ross.append((scaff,motif))
    elif scaff.split('_')[0]=='5tpj':
        ntf2.append((scaff,motif))
    else:
        nat.append((scaff,motif))
print(str(len(ross)))
print(str(len(ntf2)))
print(str(len(nat)))
# short
# 210
# 13
# 279
us=[]
um=[]
for a,b in ross:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 80
# 60

us=[]
um=[]
for a,b in ntf2:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 8
# 10

us=[]
um=[]
for a,b in nat:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 90
# 53








now long

import os
matches=[i for i in os.listdir() if i[:2]=='UM']
len(matches)
# Out[1]: 4058
scaffs=[]
c1=0
c2=0
c3=0
seqs=[]
for match in matches:
    scaff=('_').join(match.split('_')[4:6])
    motif=match.split('_')[6]
    mn=match.split('_')[1]
    seq=match.split('_')[2]
    l=[]
    for s in seq:
        if not s.isdigit():
            l.append(s)
    seqs.append(('').join(l))
    scaffs.append((scaff,motif))
    if int(mn)==1:
        c1+=1
    elif int(mn)==2:
        c2+=1
    else:
        c3+=1
print(str(c1))
print(str(c2))
print(str(c3))
print(str(len(set(seqs))))
print(set(seqs))
# 2438
# 699
# 921
# 14 really 12!
# {'WGL', 'WGF', 'WGV', 'ISI', 'ISF', 'WGS', 'WGT', 'IST', 'WGI', 'ISL', 'ISV', 'WGG'}
us=[]
um=[]
for a,b in scaffs:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 803
# 125

ross=[]
ntf2=[]
nat=[]
for match in matches:
    scaff=('_').join(match.split('_')[4:6])
    motif=match.split('_')[6]
    if scaff.split('_')[0]=='model':
        ross.append((scaff,motif))
    elif scaff.split('_')[0]=='5tpj':
        ntf2.append((scaff,motif))
    else:
        nat.append((scaff,motif))
print(str(len(ross)))
print(str(len(ntf2)))
print(str(len(nat)))
# 1507
# 422
# 2129
us=[]
um=[]
for a,b in ross:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 472
# 102

us=[]
um=[]
for a,b in ntf2:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 102
# 27

us=[]
um=[]
for a,b in nat:
    if a not in us:
        us.append(a)
    if b not in um:
        um.append(b)
print(str(len(us)))
print(str(len(um)))
# 229
# 82
#how similar are matched seqs?
l1=list(... ... ...)
l2= ...
c=0
for i in l1:
    if i in l2:
        c+=1
        print(i)
print(str(c))
# WGF
# WGI
# 2

In [23]: l1   long
Out[23]:
['WGL',
 'WGF',
 'WGV',
 'ISI',
 'ISF',
 'WGS',
 'WGT',
 'IST',
 'WGI',
 'ISL',
 'ISV',
 'WGG']

In [24]: l2
Out[24]: ['IGW', 'WSF', 'WGF', 'IGF', 'IWI', 'WSI', 'WSS', 'IWF', 'WGI', 'FSF', 'TSS']


OKAY NOW YEA I DO WANNA DESIGN ON THESE, IF NOTHING ELSE I HAVE CODE
DEVELOPED AND HANDY N ALLAT

in both a8s/38e cases it is the alcohol fragment taking
the longest to cluster, which i suppose makes sense

scp cgalvin@log2.wynton.ucsf.edu:round6/a8s/a8s_Fragment_1_cluster_statistics.pdf ~/desktop/a8s_Fragment_1_cluster_statistics.pdf
    looking at frag1 clust dists, arg is strongly preferred, makes ssense
    polar clusts should reflect this!!!!!!

'''























'''

            WRITE FILTERED HB AND NP RES TO SEPARATE PDB FILES
            SUBMIT ONE FOR POLAR CLUST
            ONE FOR SPATIAL CLUST



'''
import os
#
fog=open('/wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_2/Fragment_2_contact_fuzzball.pdb','r')
ligand_path='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb'
#
lines=[line for line in fog.readlines()]
fog.close()
all_fuzzball_lines=lines
#first gotta get indices of unique residues in all_fuzzball_lines
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(all_fuzzball_lines): #collecting start/end indices for unique res
        resnum=int(all_fuzzball_lines[index][22:26])
        resname=line[17:20]
        try:
            lastresnum=int(all_fuzzball_lines[index-1][22:26])
            lastresname=all_fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=all_fuzzball_lines[index+1][17:20]
            next_resnum=int(all_fuzzball_lines[index+1][22:26])
            if resnum!=next_resnum or resname!=nextresname:
                last=index+1
                lasts.append(last)
        except:
            last=len(all_fuzzball_lines)
            lasts.append(last)
for index,start in enumerate(starts): #put the indices together for each res
    fuzzball_residue_indices.append((start,lasts[index]))
fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
#
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY'] #ONLY TRUE hbonders AND GLY
np_residues=['TYR','TRP','ALA','LEU','ILE','VAL','PHE','MET']
#separate into polar and np lists
np_fb=[]
pol_fb=[]
for (a,b) in fuzzball_residue_indices:
    resname=all_fuzzball_lines[a][17:20]
    if resname in polar_residues:
        pol_fb.append((a,b))
    if resname in np_residues:
        np_fb.append((a,b))
    else:
        pass
#load the full ligand lines to add to fuzzball outputs
#wanna put it as last line in fuzzball files
#
f=open(ligand_path,'r')
liglines=[line for line in f.readlines()]
f.close()
#now lets build the polar fuzzball first
new_residue_number=1
clean_fuzzball_lines=[]
#
for (a,b) in pol_fb: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in all_fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:22]+' '+str(new_residue_number)+'     '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:22]+' '+str(new_residue_number)+'    '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:22]+' '+str(new_residue_number)+'   '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:22]+' '+str(new_residue_number)+'  '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number<100000:
            newline=line[0:22]+' '+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1
#now add the ligand lines and output the polar cluster
for line in liglines:
    clean_fuzzball_lines.append(line)
of=open(os.path.join('a8s_polar_contact_fuzzball.pdb'),'w')
for line in clean_fuzzball_lines:
    of.write(line)
of.close()
#okay now np fuzzball
new_residue_number=1
clean_fuzzball_lines=[]
#
for (a,b) in np_fb: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in all_fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:22]+' '+str(new_residue_number)+'     '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:22]+' '+str(new_residue_number)+'    '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:22]+' '+str(new_residue_number)+'   '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:22]+' '+str(new_residue_number)+'  '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number<100000:
            newline=line[0:22]+' '+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1
#now add the ligand lines and output the polar cluster
for line in liglines:
    clean_fuzzball_lines.append(line)
of=open(os.path.join('a8s_nonpolar_contact_fuzzball.pdb'),'w')
for line in clean_fuzzball_lines:
    of.write(line)
of.close()




'''
tryna submit new pol clust script for this pol fb
compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/bsff_cluster_polar.py a8s yes Fragment_2 /wynton/home/kortemme/cgalvin/round6/a8s/a8s_polar_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params
qstat -j "$JOB_ID"
'>run_cp_$compound_dir.sh
chmod ugo+x run_cp_$compound_dir.sh
qsub -cwd -l mem_free=32G -o cluster_output -e cluster_output run_cp_$compound_dir.sh


OKAY IT APPEARS TO BE WORKING, NOW I JUST WANNA DO THE SAME THING FOR THE OTHER
A8S POLAR CLUSTERS RIGHT,
at which point i will be able to analyze the cluster stat plots, look at which
residues are favored at the top of output jsons, and finally enumerate polar motifs
which I can then try to match
    if I am getting a lot of matches again, i can cluster the np in a new separate way
    and then try a hybrid motif thing
        ONCE I SEE WHATS UP WITH THIS, I CAN DO SAME WITH 38E->MATCH+DESIGN+AF2/MPNN




submit other polar clusts and match + design
    maybe hybrid depending on match results, but get some
    designs with great binding metrics,
    then submit for af2->mpnn->af2(->relax w match csts->rosetta metrics) pipeline
        maybe redesign mpnn w cm to allow lig to readjust?
IT WOULD KINDA BE GOOD TO SEE THAT TEST CASE RESULTS TURN OUT GOOD
FOR NEW POLAR CLUSTERING WORKFLOW
    ideally i can make just script.py target
    like other bsff scripts and it will automatically iterate fragments
    and create systematically named fuzzballs,
        then likewise with new pol clustering script
        as well as an additional new spatial clustering script i need to make


########################submitting script to process fbs into p/np
compound_dir='a8s'
echo '
#!/bin/bash
time python ~/BSFF/process_filtered_fuzzballs.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_1/Fragment_1_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_1
qstat -j "$JOB_ID"
'>pff1_$compound_dir.sh
chmod ugo+x pff1_$compound_dir.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output pff1_$compound_dir.sh

echo '
#!/bin/bash
time python ~/BSFF/process_filtered_fuzzballs.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_3/Fragment_3_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_3
qstat -j "$JOB_ID"
'>pff2_$compound_dir.sh
chmod ugo+x pff2_$compound_dir.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output pff2_$compound_dir.sh






########################cluster polar residues for frags 1/3 (2 still goin)

compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/bsff_cluster_polar.py a8s yes Fragment_1 /wynton/home/kortemme/cgalvin/round6/a8s_Fragment_1_polar_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params
qstat -j "$JOB_ID"
'>run_cp1_$compound_dir.sh
chmod ugo+x run_cp1_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output run_cp1_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/bsff_cluster_polar.py a8s yes Fragment_3 /wynton/home/kortemme/cgalvin/round6/a8s_Fragment_3_polar_contact_fuzzball.pdb  /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params
qstat -j "$JOB_ID"
'>run_cp3_$compound_dir.sh
chmod ugo+x run_cp3_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output run_cp3_$compound_dir.sh


scp cgalvin@log2.wynton.ucsf.edu:round6/a8s/a8s_Fragment_2_polar_clusters_plots_liberal.pdf ~/desktop/a8s_Fragment_2_polar_clusters_plots_liberal.pdf















'''


















'''
ALTERNATIVE CLUSTERING METHOD THAT WORKS
            much,much faster than og also
            currently only clustering polar contacts,
            still need to write script to spatial cluster
            np contacts

OKAY, PLAN WILL BE TO SPLIT POLAR FUZZBALL PDB INTO A FUZZBALL
FOR EACH INDIVIDUAL RESIDUE TYPE
    - maybe the problem is just loading really large pdb files into pyrosetta?
    idfk, everything seems to work fine with my test pdb file of 2k lines
    and the og setup worked fine when it was loading a bunch of smaller cluster pdb
    files instead of these one very large polar fuzzball files, so i will give it a shot


time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_1/Fragment_1_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_1
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_2/Fragment_2_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_2
time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_3/Fragment_3_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_3

time python3 ~/BSFF/cluster_polar.py a8s /wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb Fragment_1 Fragment_2 Fragment_3



    core.pack.pack_missing_sidechains: {0} packing residue number 9990 because of missing atom number 10 atom name  OE1
core.pack.pack_missing_sidechains: {0} packing residue number 9991 because of missing atom number 10 atom name  NE2
core.pack.pack_missing_sidechains: {0} packing residue number 11927 because of missing atom number 10 atom name  OE1
core.pack.pack_missing_sidechains: {0} packing residue number 11928 because of missing atom number 10 atom name  NE2

            I REALLY DONT WANT IT PACKING RES WITH MISSSING ATOMS
            check if zero occ in pdb file
            ignore_zero_occupancy false
            -packing pack_missing_sidechains false

i will want nanomolar binder monomer af2 scores to compare against
'''






'''
ENUMERATING POALR MOTIFS
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json n_clust_each_frag ligand_resname nmotifres resultdir_prefix ndesired_motifs dist_cutoff_index
    dist_ind[5]=2.4 (allows d up to 2.5 ang)
    no more 1 arg
    no more 1 trp
    no more than 2 instances any other res
    residues are 'banned' from being in any more motifs once >1/4 the number of desired motifs

time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 30 a8s 3 a8s6_new_way_3respol_motifs 500 5
['LYS', 'ASN']
SER
6
GLN
6
ASN
126
LYS
126
ARG
93
                hmmm, should i add a thing to not add too many charges?
                yea fuck it ima try that
                IT WOULD STILL INCREASE MY EFFICENCY HERE to merge clusters differing by 1 param n shit
            added no more than 1 arg/lys, changed cutoff for ind res to half of all motifs ,
            will up to 50 clust each frag
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 50 a8s 3 a8s6_new_way_3respol_motifs 500 5
    interestinly, now that only 1 arg or lys is allowed per motif, asn and gln are the
    residues getting banned first, must be a lot with these, maybe 2 instances of them
405
['ASN', 'GLN']
405
SER
135
GLN
250
ARG
197
TRP
172
ASN
251
LYS
207
            okay these seems pretty balanced actually

I COULD CHANGE POSITION FILES TO INCREASE MATCHING EFFICIENCY,
ONLY ALLOW POSITIONS IN BINDING POCKET? IDK SEEMS LIKE A PAIN ACTUALLY,
I WILL PROCEED NORMAL WAY
    highest maxvmem val i see in the angry wynton email is about 37
    should i specify 38? only problem is how many cpus do we have with
    that kinda power, while so many of my jobs are so much less than that too


scp cgalvin@log2.wynton.ucsf.edu:bsff5/a8sshorthb/UM_1_I59G58W89_1_clean_3en8A_223_1.pdb ~/desktop/a8s5/UM_1_I59G58W89_1_clean_3en8A_223_1.pdb
turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait
movie record ; turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait; movie encode ~/desktop/desexspin.mp4




    NEW MATCHING STRATEGY
    will be to submit a different tas array for each motif
    i will simply submit one for now, see how many matches i get,
    and based off of that decide how many more i wanna submit, or whether it is
    maybe worth doing a hybrid cst approach now that the sequence stats look
    good again
    also i can see how much memory each individual matching job will use,
    and i can base future mem requests off of that
        idk if more ram gets used when is runs consecutive cmds on same core?

NEED TO UNZIP ROSSMAN FOLDS SO THEY ARE ALL PDB
import os
l=[i for i in os.listdir() if i[-3:]=='.gz']
for i in l:
    os.system('gunzip '+i)
'''


import os
#now submit for matching with lucs strc
lig_name='a8s'
allmotifsdir='/wynton/home/kortemme/cgalvin/round6/a8s/a8s6_new_way_3respol_motifs_polar_csts'
shellfile_suffix='_a8snew3rpsr.sh'
paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
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
    for cst_path in motif_paths[key][:1]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        of.write('\n')
        of.write('tasks=(0\n')
        for scaff in scaffolds[:-1]:
            of.write('       '+os.path.join(scaffold_path,scaff.strip('.pdb'))+'\n')
        of.write('       '+os.path.join(scaffold_path,scaffolds[-1].strip('.pdb'))+')')
        of.write('\n')
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
        of.write('\nqstat -j "$JOB_ID"')
        print(count)
        count+=1
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)
print('n motifs')
print(len(paths))
print('n scaffolds')
print(str(len(scaffolds)))
'''
clean_1mpd.pdb
1
n motifs
404
n scaffolds
1505
'''
qsub -cwd -t 1-1505 -l mem_free=32G -o cluster_out -e cluster_out _1_a8snew3rpsr.sh
'''
okay looks like we runnin without error here
okay and instantly getting a gazillion matches again,
to the point that after a couple minutes i already got motifs matching
into scaffolds over 100 times lmao
            LETS DO HYBRIDS :)))
            need to do np clustering new way,
            then maybe npffa and mc if i wanna try 2 np contacts in motifs
        AND MAYBE FORGET ABOUT MATCHING WITH ROSSMAN?idk ill keep 4 next rd but moving forward may be a thing

buns no
sub uM mean = 0.3155
buns
sub uM mean = 2.2755
nhb/pot
sub uM mean = 0.3337
buns/pot ????????
rpnp
sub uM mean = 0.3318



time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/round6/a8s/Transformed_Aligned_PDBs/Fragment_4/Fragment_4_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb a8s Fragment_4


time python3 ~/BSFF/cluster_nonpolar.py target_name fragfbdir ligand_path frag1 frag2 frag3
parent_directory='a8s'
fragfbdir='/wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs '
ligand_path='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s_0001.pdb'
f1='Fragment_4'

okay i have it working where its spittin out the clusters nice,
now i gotta think about how i wanna pick np csts to add
        maybe to add 1, i will simply select from top ten most populated clusts











ADDING 1 NP CST TO POL CSTS
'''
#writes the csts to working directory
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random

fragment_fuzzballs_dir='/wynton/home/kortemme/cgalvin/round6/a8s/a8sfragfbs/a8s_Fragment_4_residue_contacts'
ligparams=['/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params']
polar_cst_paths='/wynton/home/kortemme/cgalvin/round6/a8s/a8s6_new_way_3respol_motifs_polar_csts'

clust_paths=[]
for i in os.listdir(fragment_fuzzballs_dir):
    if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
        for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
            if clust[-3:]=='pdb':
                clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))

clust_paths_pops=[]
for clustpath in clust_paths:
    pop=int(clustpath.split('/')[-1].split('.')[0].split('_')[3])
    clust_paths_pops.append((clustpath,pop))

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
                                          distance_tolerance=0.35, angle_A_tolerance=15, angle_B_tolerance=15,
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


'''
hmmmm okay now how do i wanna define np csts for top clusts
    can score and use lowest scoring
    can simply take the first residue
i think the np contacts are sufficiently similar within each cluster that taking the first res
should be fine, so ill try that
'''
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
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
for clust,freq in nps[:20]:
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
    if sccs[0][1]<-1.5:
        accepted_clusts.append((clust,freq,sccs[0][1]))
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(sccs[0][0]).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
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
                                                       distance_tolerance=0.35, angle_A_tolerance=15, angle_B_tolerance=15,
                                                       torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                       torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                       distance_constraint_sample_number=0)
        for wp in range(freq):
            np_cst_blocks.append(cstblock)


c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        f=open(os.path.join(polar_cst_paths,cst),'r')
        lines=[line for line in f.readlines()]
        f.close()
        ofname='hybrid_'+str(c)+'.cst'
        of=open(ofname,'w')
        for line in lines:
            of.write(line)
        rblock=random.choice(np_cst_blocks)
        of.write('CST::BEGIN\n')
        of.write('\n'.join(rblock))
        of.write('\n')
        of.write('CST::END\n')
        of.write('\n')
        of.close()
        print(c)
        c+=1

'''

DAMN I AM SEEING THAT THE CST BLOCK CAN BE VERY DIFFERENT FOR NO RESIDUES IN THE SAME
CLUSTER
MAYBE IT ALSO MAKES SENSE TO CLUSTER THESE BASED OFF OF CST PARAMETERS,
EH BUT I GUESS IN THAT CASE I WOULD NEED TO INCLUDE THE OTHER 2 TORSIONS,
COULD BE VERY FEW IN SAME BINS






MATCH 3P1NP A8S HYBRIDS
'''
import os
#now submit for matching with lucs strc
lig_name='a8s'
allmotifsdir='/wynton/home/kortemme/cgalvin/round6/a8s/a8s6_new_way_3respol_motifs_polar_csts/hybridcsts'
shellfile_suffix='_a8s3p1npr6.sh'
paramspath='/wynton/home/kortemme/cgalvin/round6/a8s/Inputs/Rosetta_Inputs/a8s.params'
scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
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
    for cst_path in motif_paths[key][1:]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        of.write('\n')
        of.write('tasks=(0\n')
        for scaff in scaffolds[:-1]:
            of.write('       '+os.path.join(scaffold_path,scaff.strip('.pdb'))+'\n')
        of.write('       '+os.path.join(scaffold_path,scaffolds[-1].strip('.pdb'))+')')
        of.write('\n')
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
        of.write('\nqstat -j "$JOB_ID"')
        print(count)
        count+=1
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)
print('n motifs')
print(len(paths))
print('n scaffolds')
print(str(len(scaffolds)))
'''
clean_1mpd.pdb
1
n motifs
404
n scaffolds
1505

qsub -cwd -t 1-1505 -l mem_free=32G -o cluster_out -e cluster_out _1_a8s3p1npr6.sh
qsub -cwd -t 1-1505 -l mem_free=32G -o cluster_out -e cluster_out _4_a8s3p1npr6.sh
in termsof adding more csts,
    -i can enumerate combinatorialk combos of 2 contacts for frag4 and score them
        to see if good sterically,
    -could also add 2 hbonds to the carbonyl

scp cgalvin@Log2.wynton.ucsf.edu:round6/a8s/a8s6_new_way_3respol_motifs_polar_csts/hybridcsts/UM_1_K73N70N31V33_1_model_34395_hybrid_1_1.pdb ~/desktop/UM_1_K73N70N31V33_1_model_34395_hybrid_1_1.pdb
the hbond network looks great but the valine looks like shit the sc isnt evem
pointing at the ligand wtf









july 1 2022
DO THE MPNN NTF2 THING N JUST USE THOSE
    getting 2 many matches cus i have too many scaff nd they aint all good
try matching 4 hb
OR a better fucking np contact
OR multiple np contacts
MATCHES I JUST GOT W/ 3P1MNP HYBRIDS FOR ONLY 2 MOTIFS
    In [2]: len(os.listdir())
    Out[2]: 141


            CHANGING ADD1NPCST CODE ABOVE SO THAT BB ATOMS CANNOT
            BE FIRST ATOM IN RES CST,
            CHANGING N TOP CLUSTERS CONSIDERED FROM 10 TO 20,
            AND CHANGING THE ENERGY THRESHOLD FROM -0.5 TO -1.5 REU
    i could pretty easily use the 4 contacts returned to make pdbs of 2contacts
    and score them and look if any have no clashes and i could use for 5res cst


SOME JOBS W LIKE 64 GB MAXVMEM
        maybe i will forget about the naturally occurring scaffolds
        id rather use de novo anyway
        and these can be much bigger proteins so probably use more mem

LET ME TRY TO GET A ROUND OF MPNN NTF2 DESIGNS
    read the paper and see what they did, but i think just
    ->run af2 on current filtered ntf2s
        ->fastrelax af2 structures
( carmsd ogrosetta:unrelaxed af2,   plddt of each struc,    carmsd vs plddt,
lddt per residue plots, diff rossetta metrics ogrosetta:rosettarelaxedaf2   )
    ->run MPNN on current filtered ntf2s
        (scores,    sequence diversity)
        ->run af2 on top N scoring sequences for each
        (plddt distribution of these)
            ->take highest plddt for each bb, fastrelax, submit for another
            round of MPNN
                ->again, top N scoring sequences for each bb to AF2, hiest plddt
                each bb->fastrelax->MPNN
                    CONTINUE THIS PROCESS ITERATIVELY UNTIL


4 presentation maybe some kinda analysis of hbond ckluster stats and param stats
and how they compare to rosetta enrgies for those interactions or something






'''

turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait
movie record ; turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait; movie encode ~/desktop/a8s3p1np_ex1.mp4
