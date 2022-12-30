'''
MOLECULES USED HERE
a8s - 3 polar groups to sat
nps - 2 pol, 1 ether - nice target
3ng - big, rotatable bond, hard target
jge - (r)ketoprofen - 2 entries - 2 pol nice targ but a lil bigger than nps
dex - 5 oxygens, kinda like

'''



'''
Create a directory for the new target
'''
time python ~/desktop/BSFF/new_target.py a8s
time python ~/desktop/BSFF/new_target.py 3ng
time python ~/desktop/BSFF/new_target.py jge


'''
download a sdf/molfile of the target chemical structure (RCSB or Pubchem)
open in avogadro
    build --> add hydrogens for ph --> 7.4
    autooptimization tool --> ff=mmff94 --> start (run until dG=0)
    save as --> target_opt.mol

'''

'''
create params
compounds=("a8s" "3ng" "jge")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}_opt.mol
done

move params and pdb to Inputs/Rosetta_Inputs:
compounds=("a8s" "3ng" "jge")
for I in ${compounds[@]}; do
    mv ${I}.params ${I}/Inputs/Rosetta_Inputs/${I}.params
    mv ${I}_0001.pdb ${I}/Inputs/Rosetta_Inputs/${I}_0001.pdb
done

'''


'''
fragment search
with frags.txt files in target parent directories

compounds=("a8s" "3ng" "jge")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
    cd ..
done

'''


'''
NOW setup to split into multiple jobs for alignment on cluster
conda activate pyr37
compounds=("a8s" "3ng" "jge")
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

scp -r bsff5 cgalvin@log2.wynton.ucsf.edu:
scp -r BSFF cgalvin@log2.wynton.ucsf.edu:

submit the alignment jobs
a8s
qsub -cwd -t 1-2148 -l mem_free=2G -o cluster_output -e cluster_output run_align.sh
jge
qsub -cwd -t 1-785 -l mem_free=2G -o cluster_output -e cluster_output run_align.sh
3ng
qsub -cwd -t 1-1009 -l mem_free=2G -o cluster_output -e cluster_output run_align.sh


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


compounds=("a8s" "3ng" "jge")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/BSFF/setup_filter_contacts.py
    cd ..
done

3ng
1839

jge
1694

a8s
5651

chmod ugo+x run_filter_array.sh

qsub -cwd -t 1-1694 -l mem_free=1G -o cluster_output -e cluster_output run_filter_array.sh
qsub -cwd -t 1-1839 -l mem_free=1G -o cluster_output -e cluster_output run_filter_array.sh
qsub -cwd -t 1-5651 -l mem_free=1G -o cluster_output -e cluster_output run_filter_array.sh

REDO
w occ 0.25, pol max dist 3.6
set 4g mem for setup and changed to 500 strc per job from 100
    maybe i wasnt including heavy atoms in filter script, so 2.5 ang between heavy atoms (no H)
    was killing me
        i will actually then make both dist max 4.1 ang in filter script, because im
        later filtering on hbond dist less than 2.6 more carefully anyway
3ng
368
jge
339


qsub -cwd -t 1-368 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh
qsub -cwd -t 1-339 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh
qsub -cwd -t 1-5651 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh

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


time python ~/BSFF/consolidate_filtered_contacts.py









NOW OUTSIDE OF PARENT DIRECTORY
mkdir cluster_output


compound_dir='3ng'
compound_dir='jge'
compound_dir='a8s'


echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/cluster_filtered_contacts.py '$compound_dir'
qstat -j "$JOB_ID"
'>clust_$compound_dir.sh

chmod ugo+x clust_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output clust_$compound_dir.sh

damn now they all taking hellaaaaa long, who knows how long a8s will be since it
already took dumby long last time
i guess it makes sense because there will be a lot more contacts this time,
which is good, and exactly what i want, but hopefully its not like a week or two
long runtime or anything crazy like that
    the problem is this isnt really parallelizable since i need to look at all contacts
    to accurately assign clusters
        however, i can simply change it to exclude polar residues since im gonna
        do them a different way anyway, and then i can change reclustpolar script
        to act on like, maybe a single pdb file for all polar residues or one
        for each polar res type or something
            nah not as simply as i thought, many res need both pol and
            np clustering, ie tyr, trp, maybe others. +the script is very gnarly
            to play around with, i think im just at the mercy of waiting for it to
            finish the normal way

CLUSTERING FINALLY FINISHING UP AFTER ABOUT 5 DAYS
    jge taking the longest this timne,interestingly enough
    a8s frags 5/6 have no clustered contacts somehow...kinda surprising to me ?

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  933.16 GiB| 1000.00 GiB|| 11583132|unlimited

beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
(base) [cgalvin@log2 ~]$ beegfs-ctl --getquota --storagepoolid=11 --uid cgalvin
      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
       cgalvin| 61046||  809.10 GiB| 1000.00 GiB||  8988994|unlimited








                TO DO
        SUBMIT OTHER MOAD FOR ANALYSIS + DO W/ GENPOT

        SVM CLASSIFIER FOR MOAD

        SETUP MOTIF GEN PIPELINE
            i think i wanna do np mc w/ hb frag exclusion, match, then match
            hbonds into resultant matches (be sure to exclude matched np from pos file)

TRY CLASSIFY AGAIN BUT REMOVE MIDDLING BINDERS, SEE IF CAN DISTINGUISH
REALLY BAD FROM REALLY GOOD
SHIT MAN IS THERE ANY? WAY I COULD MAKE CLUSTERING FASTER....?



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
https://assets.researchsquare.com/files/rs-1666858/v1/901b0566-ca2a-4fb5-81f4-7a2117fecda6.pdf?c=1652899389:

We found that for proteins with good multisequence alignment (MSA) depths, including
large proteins and protein complexes, the AF2-scores are highly correlated with the root mean square
fluctuations (RMSF) calculated from MD simulations.
 For proteins with little or no MSA hits (the IDP and
randomized protein), the AF2-scores do not correlate with the RMSF from MD, especially for the
intrinsically disordered proteins (IDPs). Our results indicate that the protein structures predicted by AF2
also convey information of the residue flexibility, i.e., protein dynamics.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





setup mc np motif assembly + recluster polar stuff

scp npffa.py cgalvin@log2.wynton.ucsf.edu:BSFF/
scp np_mcma.py cgalvin@log2.wynton.ucsf.edu:BSFF/

            WAIT YOU CANT DO DOUBLE MATCHING
            CUS CANT KEEP LIGAND THERE LOL

            radakas design analysis
            matching all 3 targs, hybrid motifs after all



RECLUSTER POLAR

make a different copy of bsff_recluster_polar.py for each ligand that
specifies the relevant ligand hbonding atoms in the dictionary at the top of
the script (relevant_atoms)

3ng
relevant_atoms={'Fragment_1':['O1','O2'],
                'Fragment_2':['N2','H1']}
3=PYRIDINE
4=CHLOROBENZENE
5=PYRIDINE
6=BENZENE


jge
relevant_atoms={'Fragment_1':['O1','O3'],
                'Fragment_2':['O2']}
3=BENZENE
4=BENZENE

a8s
relevant_atoms={'Fragment_1':['O4','O3'],
                'Fragment_2':['O1','H1'],
                'Fragment_3',['O2']}
FRAG4=CYCLOHEXANONE

            FOR SOME OF THESE NP FRAGS, ITS DEFINITELY POSSIBLE TO CONTACT FROM MULT
            ANGLES, IE RINGS IN JGE, (MAYBE) RING IN A8S

okay okay, so recluster polar frags, get json, use polarmotif script to generate a
bunch filtered for multiarg pure pol
then, npffa the other frags, npmca with the desired number of motif res, generate
a bunch of np cst this way
            i will have M polar csts and N np CSTS
            i could either sample random combos of np and polar
            or i could enumerate the top MxN
                4reference 100x100=10,000



                OUTSIDE PARENT DIRS
compound_dir='3ng'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rp_3ng.py 3ng Fragment_1 Fragment_2
qstat -j "$JOB_ID"
'>run_rp_$compound_dir.sh
chmod ugo+x run_rp_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rp_$compound_dir.sh


compound_dir='jge'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rp_jge.py jge Fragment_1 Fragment_2
qstat -j "$JOB_ID"
'>run_rp_$compound_dir.sh
chmod ugo+x run_rp_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rp_$compound_dir.sh


HAVE TO DELETE FOLDERS FOR FRAGS 5/6 IN A8S CUS THEY DIDNT RETURN ANYTHING
compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rp_a8s.py a8s Fragment_1 Fragment_2 Fragment_3
qstat -j "$JOB_ID"
'>run_rp_$compound_dir.sh
chmod ugo+x run_rp_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rp_$compound_dir.sh


hmmmm top hb clussters seem very weird, different, lots of trp
now instead of arg, redo with og dist?

            ON FILTER SCRIPT
            CHANGE OCC 0.25 TO 0.5
            MINIMUM POLAR DIST BACK TO 3.6 FROM 2.6
                my guess is that having it lower filtered out a ton of
                polar residues for some reason
                this is not good because ive seen majority of hb at 2 angs so why
                are they disappearing?idk but i will retry with these og settings







MAKE FAVORABLE NP MOTIFS
its only using contacts for nonpolar fragments because only those are evaluated in ffa script
#time ipython ~/desktop/prj/bsff/bsff_scripts/fuzzball_for_assembly.py a8s frag1 frag2 frag3 frag4

            SHOULD I NOT USE HBOND TERMS IN FFA/MCA CUS ITS NP ONLY?YE ILL DO THAT
            i think the easiest way to do this is just set the weights to 0 rather than
            trying to delete instanbces of loading or using the terms t/o the script
time python ~/BSFF/npffa.py a8s Fragment_4
time python ~/BSFF/npffa.py jge Fragment_3 Fragment_4
time python ~/BSFF/npffa.py 3ng Fragment_3 Fragment_4 Fragment_5 Fragment_6


compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/npffa.py a8s Fragment_4
'>npffa_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output npffa_$compound_dir.sh

compound_dir='jge'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/npffa.py jge Fragment_3 Fragment_4
'>npffa_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output npffa_$compound_dir.sh

compound_dir='3ng'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/npffa.py 3ng Fragment_3 Fragment_4 Fragment_5 Fragment_6
'>npffa_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output npffa_$compound_dir.sh





                MCMA NP
#time python ~/desktop/prj/bsff/bsff_scripts/mc_motif_assembly.py -d 'a8s' -w 1.0 -t 1000 -r 1000 -n 4000 -p test -m 2/3/4

3ng
relevant_atoms={'Fragment_1':['O1','O2'],
                'Fragment_2':['N2','H1']}
3=PYRIDINE
4=CHLOROBENZENE
5=PYRIDINE
6=BENZENE


jge
relevant_atoms={'Fragment_1':['O1','O3'],
                'Fragment_2':['O2']}
3=BENZENE
4=BENZENE

a8s
relevant_atoms={'Fragment_1':['O4','O3'],
                'Fragment_2':['O1','H1'],
                'Fragment_3',['O2']}
FRAG4=CYCLOHEXANONE

                    okay, for a8s, i will try just 1 nop res w/ the 3 hbs,
                    which i have done before and know should work,
                    but i will additionally do 2 residue motifs,
                    and shit should i try 3 also just to be ambitious? lemme
                    see what im doing with the others first
                    actually i think i fucked up and didnt make npmcma script
                    work for 1 res, ill have to do the old fashion way,
                    but still proceed with trying 2 here then

            just gonna output top 1k motifs and maybe filter for diversity,
            or scores or something, idk but ill figure it out from there

compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 100 -r 1000 -n 1000 -p a8s_try1_2resnp -m 2 -f 0 -g 0
'>npmcma_$compound_dir.sh

qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output npmcma_$compound_dir.sh


                jge only got 2 hb, ill try both 2 and 3 res np

compound_dir='jge'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 100 -r 1000 -n 1000 -p jge_try1_2resnp -m 2 -f 0 -g 0
'>npmcma_$compound_dir.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output npmcma_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 1000 -r 1000 -n 1000 -p jge_try1_3resnp -m 3 -f 0 -g 0
'>npmcma3r_$compound_dir.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output npmcma3r_$compound_dir.sh


            FUUUUUUUUUUUCK JGE DIDNT WORK
            Traceback (most recent call last):
  File "/wynton/home/kortemme/cgalvin/BSFF/np_mcma.py", line 124, in <module>
    pose_from_file(fuzzball_pose, fuzzball)
RuntimeError:

File: /home/benchmark/rosetta/source/src/core/kinematics/Stub.cc:105
Error in Stub::from_four_points():

File: /home/benchmark/rosetta/source/build/PyRosetta/linux/clang-3.4.2/python-3.6/release.serialization.thread/source/src/numeric/xyzVector.hh:670
Cannot create normalized xyzVector from vector of length() zero.  Error occurred when trying to normalize the vector between points A and B.  A=[-0.820000,2.000000,0.000000], B=[-0.820000,2.000000,0.000000].

real	489m24.094s
user	468m17.799s
sys	20m41.871s

                    okay, am i adding ter lines in the input fuzzball?
                    could be something weird with thinking atoms are bonded
                    when they are not
                            YES I AM ADDING TER LINES BETWEEN RESIDUES
                    otherwise, i need to go through the pdb file and
                    find the bonded atoms that have the same coordinates
                        OKAY, IF IT HAS TO BE BONDED ATOMS
                        THEN I CAN JUST CHECK COORDINATES OF ATOMS WITHIN
                        EACH RESIDUE, IF THERE ARE SAME COORD ATOMS I
                        DELETE THE RESIDUE

scp cgalvin@log2.wynton.ucsf.edu:bsff5/jge/clean_fuzzball_for_assembly_np.pdb ~/desktop/clean_fuzzball_for_assembly_np.pdb

scp -r cgalvin@log2.wynton.ucsf.edu:bsff5/a8s/a8s_try1_2resnp_motif_pdbs/clean_motifs ~/desktop/a8s2rmotifs
            FIRST OF ALL ONLY PRODUCED A HANDFUL OF MOTIFS HERE
                why are there so few?
                only 28
                guessing a lot were redundant and it only spit out nr set out of 1k desired?
                so top 1k only have 28 nr motifs?
            SECOND OF ALL SOME OF THEM LOOK OKAUY, BUT SOME OF THEM LOOK
            DUMB IN THAT IT IS STILL PUTTING THE OPPOSITE CHARGES OF
            THE BACKBONE GROUPS NEAR TO EACH OTHER
            I NEED TO TRY THIS AGAIN BUT WITH JUST FAATR/FAREP I THINK,
            FORGET ABOUT POLAR STUFF
                    ALTHOOOOOOOUGH, EVEN BETTER YET WOULD PROBABLY BE TO
                    SIMPLY MATCH POLAR INTERACTIONS AND THEN BUILD FROM
                    OTHER BINDING SITE POSITIONS IDEALIZED INTERACTIONS
                    IE SOME KIND OF SPECIAL ROT DESIGN METHOD
import json
df=open("a8s_try1_2resnp_motif_solutions.txt",'r')
data=json.load(df)
df.close()
l=data['Solutions']
l2=[]
for ll in l:
    l2.append((ll[0],ll[1]))
l3=set(l2)

a8s2r
In [2]: len(l)
Out[2]: 28
In [10]: len(l3)
Out[10]: 28
                        WAIT THERES ONLY 28 IN THE SOLUTIOINS LIST? WHY THE FUCK?
                        and some of them are redundant too still? this kinda sucks
                        the script is supposed to filter for redundancy,
                        maybe it is to a large extent but some are slipping by?idk
3nf3r
In [8]: len(l)
Out[8]: 992
In [7]: len(l3)
Out[7]: 839
                YEAH OKAY ADDS UP, AND ACTUALLY THESE 28 ARENT UNIQUE,
                SOME ARE PALINDROMES
l4=[]
for a,b in list(l3):
    if (a,b) not in l4:
        if (b,a) not in l4:
            l4.append((a,b))
In [14]: len(l4)
Out[14]: 14
3ng3r
In [10]: len(l4)
Out[10]: 837
                DAMN ONLY GETTING 14 UNIQUE MOTIFS OUT OF 1K IS CRAZY, THAT SHIT REALLY
                REAAAALLLY CONVERGED ON VERY FEW SOLNS (this for a8s 2resmotifs)
                for 3ng 3res motifs:
                In [2]: len(os.listdir())
                Out[2]: 992
                AHHHHHH WAIT THIS 3NG STUFF IS WRONG BECAUSE I DIDNT DO ALL 3 RES


OKAY, FOR NOW I CAN CONTINUE WITH THIS MATCHING PLAN LIKE THIS
BUT IM REALLY THINKING I MIGHT BE ABLE TO PUT TOGETHER A PYTHON SCRIPT THAT
MIRROS THE BEHAVIOR OF THE COUPLEDMOVES DESIGN PROTOCOL, WHILE ALSO
INCORPORATING CONTINUOUS STATISTICAL SCORE BONUSES FOR MY SPECIAL ROTAMERS
    it looks like


f=open('clean_fuzzball_for_assembly_np.pdb','r')
lines=[line for line in f.readlines()]
f.close()
terindices=[]
for i,e in enumerate(lines):
    if e[:3]=='TER':
        terindices.append(i)

# In [15]: len(terindices)
# Out[15]: 205959
# In [2]: len(terindices)
# Out[2]: 205880
def a(p1,p2,p3):
    x1=p1[0]
    y2=p2[1]
    y3=p3[1]
    x2=p2[0]
    y1=p1[1]
    x3=p3[0]
    calc=0.5 * ((x1*(y2 - y3)) + (x2 * (y3 - y1)) + (x3 * (y1 - y2)))
    if calc==0.:
        return True
    else:
        return False
badres=[]
c=1
for i,e in enumerate(terindices[:-1]):
    currentresatoms=lines[e+1:terindices[i+1]]
    atomcoords=[]
    print(c)
    c+=1
    for atom in currentresatoms:
        x=float(atom[31:39].strip(''))
        y=float(atom[39:47].strip(''))
        z=float(atom[47:54].strip(''))
        if (x,y,z) not in atomcoords:
            atomcoords.append((x,y,z))
        # if (x,y,z)==(-0.820,2.000000,0.000000):
        #     print('this one')
        else:
            # print('caught redundant')
            badres.append((e+1,terindices[i+1]))
        for p1 in atomcoords:
            for p2 in atomcoords:
                for p3 in atomcoords:
                    if p1!=p2:
                        if p1!=p3:
                            if p2!=p3:
                                if a(p1,p2,p3):
                                    # print('caught colinear')
                                    # print(p1)
                                    # print(p2)
                                    # print(p3)
                                    badres.append((e+1,terindices[i+1]))
                    else:
                        badres.append((e+1,terindices[i+1]))
                        else:
                            badres.append((e+1,terindices[i+1]))
                            else:
                                badres.append((e+1,terindices[i+1]))

sbr=list(set(badres))
print(str(len(sbr)))
# In [32]: len(sbr)
# Out[32]: 79
sbrr=[]
for a,b in sbr:
    for i in range(a,b+1):
        sbrr.append(i)

newlines=[]
for i,e in enumerate(lines):
    if i not in sbrr:
        newlines.append(e)

of=open('cleaned_clean_fuzzball_for_assembly_np.pdb','w')
for line in newlines:
    of.write(line)
of.close()

scp cleaned_clean_fuzzball_for_assembly_np.pdb cgalvin@log2.wynton.ucsf.edu:
        RENAMED THE OG FUZZBALL AND PUT THIS IN THE A8S DIR WITH PROPER NAME
                        OKAY SOMEHOW CHECKING ONLY WITHIN RES FOR REDUNDANT
                        COORDS IS NOT RETUYRNING ANYTHING, I WILL NOW CHECK
                        ACROSS ALL RSIDUES, BUT IN THIS CASE I WOULDNT EXPECT IT
                        TO CAUSE A PROBLEM? IDK
                        UNLESS ITS NOT REDUNDANT COORDS BUT RATHER COLINEAR ATOMS WITHIN
                        RESIDUES, BUT HOW COULD THAT BE THE CASE?
                        ALSO IN THE ERROR TEXT ITS SAYING ITS RESULTING FROM
                        TRYING TO COMPUTE SOMETHING FOR THE SAME POINTS
d['x']=line[31:39]
d['y']=line[39:47]
d['z']=line[47:54]


                3ng really needs at least 3, so i try that
compound_dir='3ng'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 1000 -r 1000 -n 1000 -p 3ng_try1_3resnp -m 3 -f 0 -g 0
'>npmcma_$compound_dir.sh

qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output npmcma_$compound_dir.sh



            off rip it appears to be running ...




                    WOULD BE INTERESTING TO DO A DEEPER ANALYSIS
                    OF THE HYDROGEN BONDS NOW


ALSO I THINK MAYBE I JUST WANT TO CONCATENATE NP AND HB MOTIFS FROM DIRECTOREIS OF BOTH
OF THEM

redoing 3ng/a8s with just faatr/farep

compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 100 -r 1000 -n 1000 -p a8s_try2_2resnp -m 2 -f 0 -g 0 -c 0 -e 0
'>mca_$compound_dir.sh

qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output mca_$compound_dir.sh

compound_dir='3ng'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 1000 -r 1000 -n 1000 -p 3ng_try2_3resnp -m 3 -f 0 -g 0 -c 0 -e 0
'>mca_$compound_dir.sh

qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output mca_$compound_dir.sh


now trying jge again without colinear atoms

compound_dir='jge'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 1000 -r 500 -n 1000 -p jge_try2_2resnp -m 2 -f 0 -g 0 -c 0 -e 0
'>mca2r_$compound_dir.sh

qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output mca2r_$compound_dir.sh

echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/np_mcma.py -d '$compound_dir' -w 1.0 -t 1000 -r 1000 -n 1000 -p jge_try2_3resnp -m 3 -f 0 -g 0 -c 0 -e 0
'>mca3r_$compound_dir.sh

qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output mca3r_$compound_dir.sh


THIS ERROR FROM FIRST JGE RUN
core.kinematics.tree.Atom_: {0} [ ERROR ] Issue getting stub for atom  atomno= 2 rsd= 169843  -- possibly due to degenerate/colinear atoms:
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 2 rsd= 1  --     -0.7920000000000000       4.096000000000000       3.390000000000000
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 1 rsd= 1  --     -0.5300000000000000       5.002000000000000       4.486000000000000
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 10 rsd= 1  --     -0.3338068592337709     -0.9426414910921784  -1.154402884809525E-16
MAKE IT SEEM LIKE MAYBE THE PROBLEMATIC ATOMS DO HAVE DIFFERENT COORDS,
BUT THEY ARE INDEED COLINEAR, WHICH FUCKING SUCKS CUS I DONT EVEN REALLY KNOW HOW TO
TEST FOR THAT
            ill just use some python code from google to do it or something
            and ill run this first before checking all atoms for degeneracy again
            cus i think in this case they definitely have to be bonded to
            one another

scp -r cgalvin@Log2.wynton.ucsf.edu:bsff5/a8s/a8s_try2_2resnp_motif_pdbs/clean_motifs ~/desktop/a8st2


STILL GETTING JGE VECTOR ERROR!!!!///??????
 Issue getting stub for atom  atomno= 2 rsd= 169783  -- possibly due to degenerate/colinear atoms:
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 2 rsd= 1  --     -0.7920000000000000       4.096000000000000       3.390000000000000
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 1 rsd= 1  --     -0.5300000000000000       5.002000000000000       4.486000000000000
core.kinematics.tree.Atom_: {0} [ ERROR ] 	  atomno= 10 rsd= 1  --     -0.3338068592337709     -0.9426414910921784  -1.154402884809525E-16

            hmmmmm how can i catch this in my clean up script? wtf
            WHAT IS THE VALUE RETURNED BY MY A FUNCTION WHEN I PASS THESE 3 POINTS?






ENUMERATE POLAR MOTIFS


time python ~/BSFF/enumerate_polar_motifs.py a8s_polar_contacts_clustered.json 20 a8s nmotifres(2/3/4)



OUTSIDE PARENT DIR
compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/enumerate_polar_motifs.py a8s_polar_contacts_clustered.json 20 a8s 3 a8s3rm
'>epm_$compound_dir.sh
qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output epm_$compound_dir.sh

compound_dir='jge'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/enumerate_polar_motifs.py jge_polar_contacts_clustered.json 20 jge 2 jge2rm
'>epm_$compound_dir.sh
qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output epm_$compound_dir.sh

compound_dir='3ng'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/enumerate_polar_motifs.py 3ng_polar_contacts_clustered.json 20 3ng 2 3ng2rm
'>epm_$compound_dir.sh
qsub -cwd -l mem_free=2G -o cluster_output -e cluster_output epm_$compound_dir.sh


OKAY I NEED TO CHANGE THE ENUMERATEPOLMOT SCRIPT TO PUT THE
MOTIFS IN THEIR OWN DIRECTORY
            SEEMS TO BE WORKING VERY WELL AND VERY VERY
            QUICKLY NOW WOOOOOOOOOOOO LESSGOOOOOO







CONCATENATING POLAR/NP MOTIFS
    just a8s/3ng rn cus jge npmca still broken (vector problem)

MAKING SURE NO REDUNDANT NP CSTS
import os
csts=[i for i in os.listdir() if i[-3:]=='pdb']
len(csts)
l=[]
for i in csts:
    s1=i.split('_')[0]
    s2=i.split('_')[1].split('.')[0]
    if (s1,s2) not in l:
        if (s2,s1) not in l:
            l.append((s1,s2))
        else:
            os.system('rm '+i)
# a8stry22resnp
# In [3]: len(l)
# Out[3]: 13
#
# In [4]: len(csts)
#    ...:
# Out[4]: 25

MOTIF2CST
conda activate pyr37

import os
paramspath='/wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params'
csts=[i for i in os.listdir() if i[-3:]=='pdb']
#####angles 15 deg
######dist 0.3 angs for np
for motif in csts:
    cstprefix=motif.split('.')[0]
    cstname=cstprefix+'.cst'
    cmds=['time',
          'python',
          '~/BSFF/tools/motif_to_cst.py',
          motif,
          paramspath,
          cstname]
    cmd=(' ').join(cmds)
    os.system(cmd)


CREATE HYBRID CSTS IN THEIR OWN DIRECTORY
import os
###
polarcsts='/wynton/home/kortemme/cgalvin/bsff5/a8s3rm_polar_csts'
npcsts='/wynton/home/kortemme/cgalvin/bsff5/a8s/a8s_try2_2resnp_motif_pdbs/clean_motifs'
###
hybridcsts='a8shybridmotifs'
####
####
####
os.makedirs(hybridcsts,exist_ok=True)
pc=[i for i in os.listdir(polarcsts) if i[-3:]=='cst']
npc=[i for i in os.listdir(npcsts) if i[-3:]=='cst']
# In [5]: len(pc)
# Out[5]: 2900
#
# In [6]: len(npc)
# Out[6]: 13
####
pcst_threshold=300
####
filtered_pc=[]
for pcst in pc:
    s=pcst.split('.')[0]
    if int(s)<=pcst_threshold:
        filtered_pc.append(pcst)
####
#its good to write the np ones first cus harder to find, could cause
#search to stop earlier if no hits
c=1
for npcst in npc:
    f=open(os.path.join(npcsts,npcst),'r')
    lines=[line for line in f.readlines()]
    f.close()
    for pcst in filtered_pc:
        f2=open(os.path.join(polarcsts,pcst),'r')
        lines2=[line2 for line2 in f2.readlines()]
        f2.close()
        newhybridname='hybrid_'+str(c)+'.cst'
        of=open(os.path.join(hybridcsts,newhybridname),'w')
        for el in lines:
            of.write(el)
        of.write('\n')
        for oel in lines2:
            of.write(oel)
        of.close()
        c+=1


'''


'''
matching
3900 a8s hybrid motifs w/ whole filt library
    3 polar + 2 np interactions (300 pm, 13 npm)
'''

import os
#now submit for matching with lucs strc
lig_name='a8s'
allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/a8shybridmotifs'
motifs=[i for i in os.listdir(allmotifsdir) if i[-4:]=='.cst']
paths=[os.path.join(allmotifsdir,i) for i in motifs]
print(len(paths))
motif_paths={}
motif_paths[lig_name]=paths
shellfile_suffix='_a8shm5.sh'
#
os.makedirs('cluster_out',exist_ok=True)
for key in list(motif_paths.keys()):
    paramspath='/wynton/home/kortemme/cgalvin/bsff5/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
    scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']
    lig_name=str(key)
    count=1
    for cst_path in motif_paths[key]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for scaffold in scaffolds:
            scaffold_input=os.path.join(scaffold_path,scaffold)
            pos_name=scaffold.split('.')[0]+'.pos'
            pos_path=os.path.join(scaffold_path,pos_name)
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_input,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', pos_path,
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
        count+=1
        print(count)
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)

import os
jobss=[i for i in os.listdir() if i[-2:]=='sh']
for sfn in jobss:
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=2G -o cluster_out -e cluster_out '+sfn)


~/main/source/bin/match.linuxgccrelease -in:file:s=/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered/model_23835.pdb.gz -in:file:extra_res_fa=/wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params -packing:extrachi_cutoff=0 -packing:use_input_sc -packing:ex1 -packing:ex2 -match:lig_name=a8s -match:geometric_constraint_file=/wynton/home/kortemme/cgalvin/bsff5/a8shybridmotifs/hybrid_1158.cst -match:scaffold_active_site_residues=/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered/model_23835.pos -match:consolidate_matches -match:output_matches_per_group=1 -match:output_format=PDB -match:match_grouper=SameSequenceGrouper -match:enumerate_ligand_rotamers=false









'''


            BENCHMARK SEED OF AN IDEA
redesigning binding site with bsff but in existing binding site where i preserve
one of the original interactions


RADHIKA XTAL ANALYSIS
    xtal strc of nat binder
    xtal strcs of suz nat binders
MATCHING 3 TARGS
AF2 SUBMISSION
OTHER MOAD STUFF
    apo scores, genpot, fragqual, allat!

'''




'''

EMAIL SALLY

CONTINUE A8S CLUSTERING N SHIT
GET JGE/3NG TO MATCHING
CONTINUE MOAD ANALYSIS - og params analysis (mcs, lig vs ref)
    do i need to relax in the sf im gonna analyze with? probably actually


REDOING A8S/JGE/3NG FILT-->Cluster
    but i still need to decide exactly what id like to do when i get my clusters
    and whatnot. how do i wanna do motif gen, matching, design
CONTINUE WITH MOAD
    analyze monomer results, submit dimers/double dimers, look at ddg vs affinity
        get all the code for analysis setup and ready to go
    setup ligand sf run - pretalaris and whatnot
    setup genpot run - add mcs last cus it may not work, get rid of sc cus i know it dont work




SALLY ARRIVES LATE, 6/12
 sboerma17@gmail.com

'''





'''
just make some motifs and try to get some matches w/ all kinds of scaff for
3 diff molecules then redesign with ifc + ex rotatmers n everything,
+CM w/ proper bb flex allowed, allat, get some good expressible designs

af2 predictions of moad binders
    + whole moad sideproject generally
        can develop workflow to process resultant data while waiting on full results



proteinMPNN design
    benchmarking?
        could rebenchmark fd/cm/gp vs ref with ifc option at the same time !



OKAY WHAT DO I HAVE AGAIN FOR THESE MOLECULES
    3ng
    a8s
    jge
ive got polar clusters for all
jge np motifs is not working because fuzzball has these messed up vectors problem








f=open('clean_fuzzball_for_assembly_np.pdb','r')
lines=[line for line in f.readlines()]
f.close()
terindices=[]
for i,e in enumerate(lines):
    if e[:3]=='TER':
        terindices.append(i)

# In [15]: len(terindices)
# Out[15]: 205959
# In [2]: len(terindices)
# Out[2]: 205880
def a(p1,p2,p3):
    x1=p1[0]
    y2=p2[1]
    y3=p3[1]
    x2=p2[0]
    y1=p1[1]
    x3=p3[0]
    calc=0.5 * ((x1*(y2 - y3)) + (x2 * (y3 - y1)) + (x3 * (y1 - y2)))
    if calc==0.:
        return True
    else:
        return False
badres=[]
c=1
for i,e in enumerate(terindices[:-1]):
    currentresatoms=lines[e+1:terindices[i+1]]
    atomcoords=[]
    print(c)
    c+=1
    for atom in currentresatoms:
        x=float(atom[31:39].strip(''))
        y=float(atom[39:47].strip(''))
        z=float(atom[47:54].strip(''))
        if (x,y,z) not in atomcoords:
            atomcoords.append((x,y,z))
        # if (x,y,z)==(-0.820,2.000000,0.000000):
        #     print('this one')
        else:
            # print('caught redundant')
            badres.append((e+1,terindices[i+1]))
        for p1 in atomcoords:
            for p2 in atomcoords:
                for p3 in atomcoords:
                    if p1!=p2:
                        if p1!=p3:
                            if p2!=p3:
                                if a(p1,p2,p3):
                                    # print('caught colinear')
                                    # print(p1)
                                    # print(p2)
                                    # print(p3)
                                    badres.append((e+1,terindices[i+1]))
                    else:
                        badres.append((e+1,terindices[i+1]))
                        else:
                            badres.append((e+1,terindices[i+1]))
                            else:
                                badres.append((e+1,terindices[i+1]))

sbr=list(set(badres))
print(str(len(sbr)))
# In [32]: len(sbr)
# Out[32]: 79
sbrr=[]
for a,b in sbr:
    for i in range(a,b+1):
        sbrr.append(i)

newlines=[]
for i,e in enumerate(lines):
    if i not in sbrr:
        newlines.append(e)

of=open('cleaned_clean_fuzzball_for_assembly_np.pdb','w')
for line in newlines:
    of.write(line)
of.close()





scp -r cgalvin@log2.wynton.ucsf.edu:bsff5/3ng/3ng_try2_3resnp_motif_pdbs/clean_motifs ~/desktop/3ng3rm



OKAY A FIRST THING TO PARSE OUT IS THE DISCREPENCY BETWEEN
ROSETTA ENERGY AND CONTACT FREQUENCY IN PDB
    AND HAVE TO BE SURE TO CONTROL FOR HOMOLOGY N SHIT SOMEHOW
        ie is an interaction observed a lot cus its good for that chemical group
        or because there 1032458394 instances of a certain molecule with that chemical group
        and in the specific context of that molecule the interaction is great
BENCHMARKING
    can i reproduce 'idealized' binding sites

MY NANOMOLAR MOAD BINDERS ARE THE ONLY BINDERS I CAN BE SURE ARE
'IDEALIZED' TO STUDY
ALL THE OTHERS IN PDB IF I DONT KNOW AFFINITY I DONT ACTUALLY KNOW HOW
GOOD THEY ARE
            this is where i should look for benchmarking
            perhaps multiple binders for same molecule amongst this set



compound_dir='jge'
echo '
#!/bin/bash
time python ~/BSFF/check_fuzzball_vectors.py
qstat -j "$JOB_ID"
'>cfv_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output cfv_$compound_dir.sh



okay why dont i try 3NG
    3rnp + 2hb


3ng
relevant_atoms={'Fragment_1':['O1','O2'],
                'Fragment_2':['N2','H1']}
3=PYRIDINE
4=CHLOROBENZENE
5=PYRIDINE
6=BENZENE


jge
relevant_atoms={'Fragment_1':['O1','O3'],
                'Fragment_2':['O2']}
3=BENZENE
4=BENZENE

a8s
relevant_atoms={'Fragment_1':['O4','O3'],
                'Fragment_2':['O1','H1'],
                'Fragment_3',['O2']}
FRAG4=CYCLOHEXANONE




CONCATENATING POLAR/NP MOTIFS for 3ng

MAKING SURE NO REDUNDANT NP CSTS OR ALL 3 RES THE SAME RESTYPE
import os
csts=[i for i in os.listdir() if i[-3:]=='pdb']
len(csts)
l=[]
# for i in csts:
#     s1=i.split('_')[0]
#     s2=i.split('_')[1].split('.')[0]
#     if (s1,s2) not in l:
#         if (s2,s1) not in l:
#             l.append((s1,s2))
#         else:
#             os.system('rm '+i)
for i in csts:
    s1=i.split('_')[0]
    s2=i.split('_')[1]
    s3=i.split('_')[2].split('.')[0]
    if (s1,s2,s3) not in l:
        if (s3,s2,s1) not in l:
            if (s1,s3,s2) not in l:
                if (s3,s1,s2) not in l:
                    if (s2,s1,s3) not in l:
                        if (s2,s3,s1) not in l:
                            l.append((s1,s2,s3))
        else:
            os.system('rm '+i)
print(str(len(csts)))
print(str(len(l)))
from pyrosetta import *
init('-ignore_unrecognized_res -load_PDB_components False')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())#[:-1] actually i think i need the whole sequence since im passing ignoreunrec/loadpdbfalse options
    if sequence[0]==sequence[1]:
        if sequence[1]==sequence[2]:
            os.system('rm '+pdb)
for i in csts:
    fasta(i)
csts=[i for i in os.listdir() if i[-3:]=='pdb']
print(str(len(csts)))
In [12]: len(csts)
Out[12]: 790

MOTIF2CST
conda activate pyr37

import os
# paramspath='/wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params'
paramspath='/wynton/home/kortemme/cgalvin/bsff5/3ng/Inputs/Rosetta_Inputs/3ng.params'

csts=[i for i in os.listdir() if i[-3:]=='pdb']
#####angles 15 deg
######dist 0.5 angs for np
for motif in csts:
    cstprefix=motif.split('.')[0]
    cstname=cstprefix+'.cst'
    if not os.path.exists(cstname):
        cmds=['time',
              'python',
              '~/BSFF/tools/motif2cst_tolerant.py',
              motif,
              paramspath,
              cstname]
        cmd=(' ').join(cmds)
        os.system(cmd)



CREATE HYBRID CSTS IN THEIR OWN DIRECTORY
import os
polarcsts='/wynton/home/kortemme/cgalvin/bsff5/3ng2rm_polar_csts'
npcsts='/wynton/home/kortemme/cgalvin/bsff5/3ng/3ng_try2_3resnp_motif_pdbs/clean_motifs'
###
hybridcsts='3nghybridmotifs3np2p'
# ###
# polarcsts='/wynton/home/kortemme/cgalvin/bsff5/a8s3rm_polar_csts'
# npcsts='/wynton/home/kortemme/cgalvin/bsff5/a8s/a8s_try2_2resnp_motif_pdbs/clean_motifs'
# ###
# hybridcsts='a8shybridmotifs'
####
####
####
os.makedirs(hybridcsts,exist_ok=True)
pc=[i for i in os.listdir(polarcsts) if i[-3:]=='cst']
npc=[i for i in os.listdir(npcsts) if i[-3:]=='cst']
print('polar')
print(str(len(pc)))
print('nonpolar')
print(str(len(npc)))
# polar
# 400
# nonpolar
# 790
#3900 was too many last time, lemme just do like half that at most
#lets say 1500
#this time np dist is more tolerant at 0.5 angstroms
#im just gonna sample a random 2k unique combinations
sampled_combos=[]
n_desired_motifs=1000
import random
while len(sampled_combos)<n_desired_motifs:
    np_cst=npc[random.randrange(len(npc))]
    p_cst=pc[random.randrange(len(pc))]
    nn=np_cst.split('.')[0]
    pn=p_cst.split('.')[0]
    combo_id=('_').join([nn,pn])
    if combo_id not in sampled_combos:
        sampled_combos.append(combo_id)
        f=open(os.path.join(npcsts,np_cst),'r')
        lines=[line for line in f.readlines()]
        f.close()
        f2=open(os.path.join(polarcsts,p_cst),'r')
        lines2=[line2 for line2 in f2.readlines()]
        f2.close()
        newhybridname=combo_id+'.cst'
        of=open(os.path.join(hybridcsts,newhybridname),'w')
        for el in lines:
            of.write(el)
        of.write('\n')
        for oel in lines2:
            of.write(oel)
        of.close()

'''
import os
#now submit for matching with lucs strc
lig_name='3ng'
allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/3nghybridmotifs3np2p'
shellfile_suffix='_3nghm.sh'
#
motifs=[i for i in os.listdir(allmotifsdir) if i[-4:]=='.cst']
paths=[os.path.join(allmotifsdir,i) for i in motifs]
print(len(paths))
motif_paths={}
motif_paths[lig_name]=paths
#
os.makedirs('cluster_out',exist_ok=True)
for key in list(motif_paths.keys()):
    paramspath='/wynton/home/kortemme/cgalvin/bsff5/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
    scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']
    lig_name=str(key)
    count=1
    for cst_path in motif_paths[key]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for scaffold in scaffolds:
            scaffold_input=os.path.join(scaffold_path,scaffold)
            pos_name=scaffold.split('.')[0]+'.pos'
            pos_path=os.path.join(scaffold_path,pos_name)
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_input,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', pos_path,
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
        count+=1
        print(count)
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)

import os
jobss=[i for i in os.listdir() if i[-2:]=='sh']
for sfn in jobss:
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)


'''
                aAHHHHHHH FUCK SOMEHOW MY NP CSTS ARE MISSING CST BLOCKS
                SOME ONLY HAVE 1 OR 2 WHEN THEY ALL SHOULD HAVE 3
                IDK HOW THIS HAPPENED BUT ILL NEED TO RUN IT BACK AND GENERATE
                MORE FROM MC :(((((

'''



'''

okay now lemme redo matching 4 a8s w only 1 np again


time python ~/BSFF/hybrid_csts_add_1_np_contact.py /wynton/home/kortemme/cgalvin/bsff5/a8s_polar_contacts_clustered.json 20 a8s 3 /wynton/home/kortemme/cgalvin/bsff5/a8s/Inputs/Rosetta_Inputs/a8s.params
        it puts all in woprking direcotyr, gotta do 1 at a time then move them all
        ie make sure no other csts in working dir before running
            ACCIDENTALLY MADE THE NP CSTS 0.25 ANG INSTEAD OF 0.35
                redoing right quick id rather get some decent matches than get nothing :(



mv *.cst a8shyb2
'''
import os
#now submit for matching with lucs strc
lig_name='a8s'
allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/a8shyb2'
shellfile_suffix='_a8s3p1.sh'
nmotifstouse=1000
#
motifs=[i for i in os.listdir(allmotifsdir) if i[-4:]=='.cst']
paths=[os.path.join(allmotifsdir,i) for i in motifs]
print(len(paths))
motif_paths={}
pathstouse=[]
import random
while len(pathstouse)<nmotifstouse:
    randmotif=paths[random.randrange(len(paths))]
    if randmotif not in pathstouse:
        pathstouse.append(randmotif)
motif_paths[lig_name]=pathstouse
#
os.makedirs('cluster_out',exist_ok=True)
for key in list(motif_paths.keys()):
    paramspath='/wynton/home/kortemme/cgalvin/bsff5/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
    scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']
    lig_name=str(key)
    count=1
    for cst_path in motif_paths[key]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for scaffold in scaffolds:
            scaffold_input=os.path.join(scaffold_path,scaffold)
            pos_name=scaffold.split('.')[0]+'.pos'
            pos_path=os.path.join(scaffold_path,pos_name)
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_input,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', pos_path,
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
        count+=1
        print(count)
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)

import os
jobss=[i for i in os.listdir() if i[-2:]=='sh']
for sfn in jobss:
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)


'''
scp -r cgalvin@log2.wynton.ucsf.edu:bsff5/a8shyb2/3ngmatchex ~/desktop/3ngmatchex

            AFTER A NIGHT OF MATCHING, ONLY 5
            BUT WHEN I PREVIOUSLY DID THIS SAME THING (3P1NP) FOR A8S
            I GOT HUNDREDS OF MATCHES STILL
            THIS IMPLIES TO ME THAT LONG RANGE HBONDS ARE EASIER TO MATCH SOMEHOW?
                OR THAT I USED 0.5 ANGS FOR NP MATCHING?
                    OR BOTH?
                            !!!!!!!OKAY NO, IT LOOKS LIKE I USED 0.25 DIST CONSTRAINT!!!!!!!!!!cd

            additionally im thinking about wy these lucs rossman folds are easier to match into
            maybe it has to do with having an open face on any pocket pretty much
                and this is why i cant achieve comparable lighphobesasas to natural binders

are binding sites inherently defects in packing? how do natural proteins compensate for this?
systematic matching thing, def ask this hb dist question, maybe look at waters in nat bs,
    also wanna look at difference in number of matches on different scaffolds and quality of resultant designs
custom moad analysis + so many other params (n charges/net charge in bs yadda yadda)
    what about hbond distances
    if lig has aro + how many aro contacts
    fasol/elec per lig polar atom?

            okay firstly, it simply just isnt enough to do look at anything on one
            single molecule, i have no choice but to use like at least 3
i hate to say it but i think i should start yet another round and be cleaner about it
this time, make clustering only np (if i didnt already?), do more analysis of interaction stats,
make it more user friendly and easier to change important options


just match short vs all polonly motifs w small set of scaffs n see what i get

hbonding depends on availiability of precise geometric placement,
but pscking should not be so sensitive to this, we should just be able to get
away with using some general rules

        IN CM VS FD BENCHMARKING, IS MUCH OF CM IMPROVEMENT
        JUST CUS LIG-PTN INTERAXNS ARE UPQEIGHTED BY A FACCTOR OF 2!?


A MOAD PLAN
AND A MATCH/DESIGN PLAN ideally allowing for benchmarking


CURENTLY UP FOR LAB MEETING JULY 7 ~2.5 weeks

wait, when i cluster hbonds, am i correcting for distance?
    NO I AM NOT
    gonna recluster including all to 3.5 angs, then i can just only
    use clusters with distance bins below a certain number if i wanna filter to short (i do)
        also may be a good idea to write csts with a bit more leeway for the hbonds,
        after all they energetically should have a little more than im giving them



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    how do i know geometric biases arent poisoning my results
        I WANT TO CLUSTER THEM MORE LIBERALLY, WIDEN BINS, THEN I WANT TO EXCLUDE CLUSTERS
        OUTSIDE OF OBSERVED FAVORABLE WINDOWS ON HISTOGRAMS, THIS WILL ENSURE I GET ENOUGH
        CONTACTS IN EA CLUSTER AND THAT BIASES ARENT RUINING, THEN I WILL WRITE
        THE CSTS TO BE MORE LIBERAL IN THE SAME MANNER AS MY NEW CLUSTERS
            i also wanna right cluster pdbs so i can see what they look like,
            ideally from og pdb files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WHEN I RUN BACK ROUND 6 i should leave hydrogen on carboxylate for search
then just at some point switch with the properly protonated one in the pdb files n shit


compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rp_a8s.py a8s Fragment_1 Fragment_2 Fragment_3
qstat -j "$JOB_ID"
'>run_rp2_$compound_dir.sh
chmod ugo+x run_rp2_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rp2_$compound_dir.sh


tried some very simple changes to polreclust to double the bin sizes, lets see if
it works

compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rcpa8slib.py a8s Fragment_1 Fragment_2 Fragment_3
qstat -j "$JOB_ID"
'>run_rp3_$compound_dir.sh
chmod ugo+x run_rp3_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rp3_$compound_dir.sh


    MAKE THIS SO FILTERS TO ONLY DESIRED BIN RANGES
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 50 a8s 3 a8spmlt1 500 5
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 60 a8s 3 a8spmlt2 500 5
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 60 a8s 3 a8spmlt3 500 5
    try 1/6
        actually in the end ill make it a certain number of the desired motifs,
        cus i think the things i actually want is that im not trying to match
        a certain res with everu soimgle m otif
          i guess ill make it like 1/2 of all motifs

359
['TRP', 'ILE', 'GLY']
360
['TRP', 'ILE', 'GLY']
361
['TRP', 'ILE', 'GLY', 'PHE']
362
['TRP', 'ILE', 'GLY', 'PHE']
363
['TRP', 'ILE', 'GLY', 'PHE']
    could only make 3 motifs without these four res while meeting all other condish
        31 minutes with all constraints including distance

            okay, it finally worked, though im not convinced i have enough seq diversity
            a mongst the whole set of motifs, eben tho ive ensured it within motifs
            a bit more (no more than one trp)
                ok i implement feature that bans res after it accounts for
                40% all contacts across all motifs

now i try allowing long too:
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 60 a8s 3 a8spmlt4 500 20
        33 minutes with no distance constraint
PHE
250
LEU
86
TRP
250
ILE
250
THR
118
SER
162
VAL
130
GLY
251
        no arg seems super weird to me, i mean just generally not enough polar res...

            I CAN DEFINITELY USE WAY MORE CLUSTERS IN ORDER TO GET TO 500 TOTAL OR WHATEVER
            SO MUCH SO (LIKE AT LEAST 200) THAT I SHOULD BE ABLE TO MAINTAIN CURR CSTS,
            MAYBE EVEN EXPAND, LIKE NOT ALLOWING TOO MANY BB HB CUS THEY HARDER TO MATCH?
            OR EVEN ALLOWING LESS OF THE SAME RESIDUES ACROSS ALL CONTACTS
                experimentally adding something which might make it faster
                by removing all elements remaining in the clust indices which
                have the index of a known bad res

experimentally deleting indices thing:
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 60 a8s 3 a8spmlt5 500 20
    ok nvm i did it way wrong, but i will actually go ahead and
            actually now i think about it maybe that wouldve worked
        do it with only 0.25 for each res and using 100 clust and dist constraints used
time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 100 a8s 3 a8spmlt5 500 5


            HMMMM OKAY MAYBE A FASTER WAY TO DO THIS IS TO LOOK AT TOP CLUSTERS FOR
            EACH FRAG INDIVIDUALLY, LOOK AT RESIDUE COUNTS, AND DELKETE CLUSTERS ENTIRELY OF RES
            THAT ARE TOO MUCH
                or record them up front and use that as i add to badres list
            OK YEA NO THE THING TO DO IS FIND LIKE TOP 10 CLUST OF FIRST 10 UNIQUE RESTYPES
            AND MAKE A NEW DF OF THAT AND JUST USE THAT

in the meantime though i am gonna go ahead and just match these two
directories, long and short filtered 3 res hb motif for a8s
 3 short 4 long
'''


import os
#now submit for matching with lucs strc
lig_name='a8s'
# allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/a8spmlt3_polar_csts'
# shellfile_suffix='_a8s3ps.sh'
# allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/a8spmlt4_polar_csts'
# shellfile_suffix='_a8s3pl.sh'
allmotifsdir='/wynton/home/kortemme/cgalvin/bsff5/a8shybridcsts'
shellfile_suffix='_a8shs.sh'
#
motifs=[i for i in os.listdir(allmotifsdir) if i[-4:]=='.cst']
paths=[os.path.join(allmotifsdir,i) for i in motifs]
print(len(paths))
motif_paths={}
motif_paths[lig_name]=paths
#
os.makedirs('cluster_out',exist_ok=True)
for key in list(motif_paths.keys()):
    paramspath='/wynton/home/kortemme/cgalvin/bsff5/'+str(key)+'/Inputs/Rosetta_Inputs/'+str(key)+'.params'
    scaffold_path='/wynton/home/kortemme/cgalvin/lucs_scaffolds/filtered'
    scaffolds=[i for i in os.listdir(scaffold_path) if i[-3:]=='.gz' or i[-3:]=='pdb']
    lig_name=str(key)
    count=1
    for cst_path in motif_paths[key]:
        sfn='_'+str(count)+shellfile_suffix
        of=open(sfn,'w')
        of.write('#!/bin/bash\n')
        for scaffold in scaffolds:
            scaffold_input=os.path.join(scaffold_path,scaffold)
            pos_name=scaffold.split('.')[0]+'.pos'
            pos_path=os.path.join(scaffold_path,pos_name)
            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', scaffold_input,
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', pos_path,
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
        count+=1
        print(count)
        of.close()
        # os.system('chmod ugo+x '+sfn)
        # os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)

import os
jobss=[i for i in os.listdir() if i[-2:]=='sh']
for sfn in jobss:
    os.system('chmod ugo+x '+sfn)
    os.system('qsub -cwd -l mem_free=4G -o cluster_out -e cluster_out '+sfn)


'''
I HATE THAT MY CLUSTERS ARE DOMINATED BY NP RES, I NEED TO DO SOME
ANALYSIS OF CLUSTER STATISTICS TO SEE WHATS GOING ON HERE, WHY ISNT ARG DOMINANT FOR CARBOXYLATE NOW?

i will compare how many matches i get between hb params and between diff scaffolds
i will do PROPER fd/cm design on filtered matches, trying to upqeight interface,
    and i will compare design metrics with nat(4 a8s and all nanobinders),
        using diff between design and input instead of raw value to filter so i keep diff scaff
i will also try MPNN while fixing crucial binding site res

i will also try hybrid motifs again with my short range and long range,
    similarly compare number of matches and do design and everything
    IF it turns out i cant match hybrids that well, maybe i can indeed look to
    a packing optimization function for help, rather than trying to match np contacts
        OKAY JUST WORKED FROM TERMINAL USING ADDNOPCIONTACT.PY SCRIPT CODE
        ADDED 1 NP TO ALL OF 262 SHORT HB CSTS
        NOW WILL MATCH above


and i will also try no bb motifs, or rather polar res only motifs
lemme start that rn, shouldnt take as long cus i exclude more res
compound_dir='a8s'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/rcp/rpa8snobb.py a8s no Fragment_1 Fragment_2 Fragment_3
qstat -j "$JOB_ID"
'>run_rpnbb_$compound_dir.sh
chmod ugo+x run_rpnbb_$compound_dir.sh
qsub -cwd -l mem_free=8G -o cluster_output -e cluster_output run_rpnbb_$compound_dir.sh
    LOOKS LIKE TRP/GLY STILL DOMINATING FRAG 1 CLUSTS, WHAT GIVES
        this my first time doing with new search function right? so
        maybe its because i have more/different contacts
            also could try adding h on search thing, but pretty sure
            i didnt add it for round4 either where i got diff clusts/hella matches
import os
l=[i for i in os.listdir()]
c=0
for i in l:
    s=i.split('_')[3]
    ss=s.split('.')[0]
    c+=int(ss)
print(str(c))
1431 og
In [2]: len(l)
Out[2]: 23974 og
52391 now
import json
clusterjson='PDB_search_results.json'
df=open(clusterjson,'r')
dff=json.load(df)
df.close()
len(dff['Fragment_1']['PDBs'])
643
62187
len(set(dff['Fragment_1']['PDBs']))
47147
len(dff['Fragment_1']['Ligands'])
len(set(dff['Fragment_1']['Ligands']))

fuck man i need to redo with opt without h then add h for search then take it away
after

            IS THIS MATCHING TROUBLE BECAUSE OF USING TOO MANY BB?
            CUS IT SEEMS EVEN LR MOTIFS ARENT GETTING THAT MANY MATCHES


time python ~/BSFF/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered_liberal.json 100 a8s 3 a8spmlt5 500 5


scp cgalvin@log2.wynton.ucsf.edu:bsff5/a8sshorthb/UM_1_I59G58W89_1_clean_3en8A_223_1.pdb ~/desktop/a8s5/UM_1_I59G58W89_1_clean_3en8A_223_1.pdb
turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait
movie record ; turn y 1 90 ; wait ; turn y -1 180; wait; turn y 1 90; wait; movie encode ~/desktop/desexspin.mp4

'''


'''


OPTIMIZED BINDERS XTAL STRUCTURES, SUBNANOMOLAR
ALANINE MUTATE + REDESIGN HOW GET BEST RECOVERY
ALANINE EXCEPT HBOND RES, SEQ RECOVERY CAN I PACK IT WELL?
    maybe perturb some torsions or something too?
AND JUST LOOK AT THEIR RULES!

how many 'hot spot' residues or like low energy interactions do
natural binders have with targets




CATALYSIS BY PROXIMITY, IE SIMPLY HOLDING TWO MOLECULES IN THE CORRECT ORIENTATION
COULD MATCH ENTIRE ACTIVE SITE INTO SCAFFOLDS PERHAPS
    also destabilization, binding a partner in a conformer which is not its
    lowest energy conformer
    preferentially stabilizing transition state
    80% of enz rxns are redox,addition/elim,hydrolysis,or decarboxylation



We found that the closer the predicted binder structure was
to the Rosetta-designed structure in Cɑ RMSD, and the higher the prediction confidence
(pLDDT; Fig 1d), the more likely a binder was to be successful (Fig S4; the pLDDT of AF2 and
RF2 were equally discriminative). These results suggest that Type 1 failures contribute to the
low success rate of binder design, and that such failures can in part be identified by
discrepancies between design models and AF2 or RF2 structure predictions.
    they also predict complexes and similarly af2 metrics are predictive of success (much more than just ddg, which is also kinda effective)
Since unlike Rosetta, ProteinMPNN keeps the protein backbone fixed, it is sensitive to the input
backbone structure quality. Inspired by the very efficient alternation between sequence
optimization and structure refinement in Rosetta flexible backbone design
13
, we evaluated
similar cycling between ProteinMPNN and Rosetta structure refinement (FastRelax), hoping to
converge on a high-quality backbone that would then allow ProteinMPNN to generate a
high-quality sequence. This hybrid ProteinMPNN/Rosetta sequence design protocol (henceforth
referred to as ProteinMPNN-FR) generated AF2 pAE_interaction<10 structures at an average
rate of ~4.3% with a throughput of 1 design per 120 CPU-s for an average efficiency of 1.4x10-6
.
The average per-target efficiency improvement of ProteinMPNN-FR over Rosetta-design is
~6-fold (Fig 2c)
    This protocol takes as input a protein complex structure. The ProteinMPNN is then fed the
complex structure with the binder sequence masked and asked to assign the binder a
sequence. The new sequence is then threaded back onto the binder structure in the complex
and the complex structure is relaxed using Rosetta FastRelax. The relaxed complex structure
can then be used as the input to the ProteinMPNN to continue the cycle. A python script to
perform this design technique is available from the authors upon reasonable request and will be
available with the published peer-reviewed manuscript.





MATCH CONTACTS INTO SEC STRC ELEMENTS FOR HALLCUINATION?
    or cluster bb coords of hi freq contacts and find consecutive pieces of sec strc?

scp cgalvin@log2.wynton.ucsf.edu:bsff5/a8s_polar_clusters_plots_liberal.pdf ~/desktop/a8s_polar_clusters_plots_liberal.pdf
scp cgalvin@log2.wynton.ucsf.edu:bsff5/a8s_polar_clusters_plots.pdf ~/desktop/a8s_polar_clusters_plots.pdf

does rosetta have any capacity to recognize the detailed energetic differences
between aromatic interactions and simple vdw forces? can it recognize a cation or anion-pi
interaction?



    REMEMBER TO ALLOW LIG FLEX IN CM RESFILE ALSO
    DESIGNING 2 HYDROGEN BONDS WITH SINGLE CAROBONYL GROUP
        how many do they have in nat subnano binders?
    CAN I LIKE, RUN OPTH ON MATCHED STRUCTURES TO GET HBONDS IN CONSENSUS BEST POSITION?
        or just check beforehand that a hydrogen bond will be considered good by chimera/rosetta?
            and only enumerate motifs using such contacts
'''
