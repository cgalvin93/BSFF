
#create a useful directory structure for target
time python ~/desktop/BSFF/new_target.py esl

'''
gotta redo with smaller pentadiol frag
'''


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

compounds=("esl")
for I in ${compounds[@]}; do
    ~/Desktop/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}.sdf
done

move params and pdb to Inputs/Rosetta_Inputs:
compounds=("esl")
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

esl
1 benzene + oh
2 pentadiol
3 benzene
4 methyl pentane
5 all carbon except benzene and methyl group

CARBONS HAVE TO BE IN ORDER LEAST TO GREATEST

compounds=("esl")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
    cd ..
done

'''



'''
NOW setup to split into multiple jobs for alignment on cluster
conda activate pyr37
compounds=("esl")
for I in ${compounds[@]}; do
    cd ${I}
    mkdir cluster_output
    mv *.pdb Inputs/Fragment_Inputs/
    mv *.mol Inputs/Fragment_Inputs/
    time python ~/desktop/BSFF/process_fragment_search.py
    cd ..
done

getting smiles...
392

getting smiles...
893

okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r esl cgalvin@log2.wynton.ucsf.edu:










CHANGED EXTRACT_ALIGN TO USE RMSD 0.5 IF LESS THAN 10 ATOMS, 1.0 IF GREATER

in target par dir on cluster:
SUBMIT ALIGNMENT JOB
qsub -cwd -t 1-893 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh


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

esl
125

322



chmod ugo+x run_filter_array.sh
qsub -cwd -t 1-322 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh






CONSOLIDATE FILTERED CONTACTS
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/consolidate_filtered_contacts.py
'>run_consoldate.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_consoldate.sh














PROCESS FILTERED CONTACTS - do for each fragment
compounds=("Fragment_1" "Fragment_2" "Fragment_3")
for I in ${compounds[@]}; do
    time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/esl/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb esl $I
done

NOTHING FOR 4 AND 5

mkdir residue_contact_fuzzballs
mv *_residue_contacts residue_contact_fuzzballs






CLUSTER ALL
time python3 ~/BSFF/spatial_cluster.py /wynton/home/kortemme/cgalvin/dog/fragfbs/dog_Fragment_1_residue_contacts/dog_Fragment_1_SER_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb

import os
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','ILE','LEU','MET','PHE','VAL','GLY']
l=[]
for i in polar_residues:
    for x in os.listdir():
        if x[-3:]=='pdb':
            if i in x:
                l.append(os.path.join(os.getcwd(),x))

f=open('sc3.sh','w')
f.write('#!/bin/bash')
f.write('\n')
f.write('source ~/anaconda3/etc/profile.d/conda.sh')
f.write('\n')
f.write('conda activate pyr37')
f.write('\n')
f.write('tasks=(0\n')
for match in l[:-1]:
    f.write('       '+str(match)+'\n')
f.write('       '+l[-1]+')')
f.write('\n')
cmd='time python3 ~/BSFF/spatial_cluster.py ${tasks[$SGE_TASK_ID]} /wynton/home/kortemme/cgalvin/esl/Inputs/Rosetta_Inputs/esl_0001.pdb'
f.write(cmd)
f.write('\nqstat -j "$JOB_ID"')
f.close()

print(len(l))

qsub -cwd -t 1-17 -l mem_free=16G sc3.sh


FUNNA MAKE SOME PLOTS ON CLUSTERING STATS
1 benzene + oh
2 pentadiol
3 benzene


OKAY RES PREFERENCES ON HBOBNDING FRAGS
'''
import os
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP'] #no his
rescounts={}
for i in os.listdir():
    if os.path.isdir(i)==True:
        rn=i.split('_')[3]
        if rn in polar_residues:
            os.chdir(i)
            l=[i for i in os.listdir() if i[-3:]=='pdb']
            count=0
            for clust in l:
                cc=clust.split('_')[-1].split('.')[0]
                ncc=int(cc)
                count+=ncc
            os.chdir('..')
            rescounts[rn]=count

'''
1
{'GLN': 56,
 'LYS': 105,
 'THR': 160,
 'ASN': 238,
 'ARG': 191,
 'SER': 174,
 'ASP': 328,
 'GLU': 384,
 'TYR': 408,
 'HIS': 350,
 'TRP': 389}

2
{'GLN': 1212,
 'LYS': 2767,
 'THR': 2786,
 'HIS': 1960,
 'SER': 4601,
 'TYR': 2122,
 'ASN': 2992,
 'ARG': 2576,
 'GLU': 4442,
 'TRP': 3449,
 'ASP': 7370}
'''

import matplotlib.pyplot as plt
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP'] #no his

# d=...
tc=0
for key in d.keys():
    tc+=d[key]
freqs=[]
for i in polar_residues:
    try:
        freqs.append((i,d[i]/float(tc)))
    except:
        pass
freqs=sorted(freqs, reverse=True, key=lambda nmem: nmem[1])
labs=[]
vals=[]
for a,b in freqs:
    labs.append(a)
    vals.append(b)

fig = plt.figure()
ax = fig.add_subplot()
ax.bar(labs,vals,color='orange',width=0.1)
ax.set_title('Hbond Residue Preferences')
ax.set_ylabel('Frequency')
ax.set_xlabel('Residue Type')
# plt.savefig('quinonerespref.pdf')
plt.savefig('secohrespref.pdf')
plt.clf()
plt.close()

'''
25_TYR_cluster_17.pdb

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_1_residue_contacts/esl_Fragment_1_TYR_contact_fuzzball_clustered/25_TYR_cluster_17.pdb ~/desktop/25_TYR_cluster_17.pdb

220_ASP_cluster_178.pdb

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_2_residue_contacts/esl_Fragment_2_ASP_contact_fuzzball_clustered/41_ASP_cluster_133.pdb ~/desktop/41_ASP_cluster_133.pdb


MAYBE IT MAKES SENSE TO CHECK WITHIN THESE CLUSTERS WHETHER OR NOT THERE ACTUALLY
IS A HYDROGEN BOND? THEN GET SAME STATS BUT FOR JUST THOSE EXAMPLES?
    yea gonna do this
AND THEN MAYBE NEXT REDO BENAENE THING BUT WITH ONLY BENZENE RING SO NO CLASHES!
'''
import os
from pyrosetta import *
init('-pack_missing_sidechains False')
#
sf=get_fa_scorefxn()

#
ligparams=['/wynton/home/kortemme/cgalvin/esl2/Inputs/Rosetta_Inputs/esl.params']
#
clusts_data={}
#
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
#
hbond_clusts=[]
for resdir in os.listdir():
    if os.path.isdir(resdir)==True:
        rn=resdir.split('_')[3]
        if rn in polar_residues:
            os.chdir(resdir)
            clusts=[i for i in os.listdir() if i[-3:]=='pdb']
            clust_data=[]
            for clust in clusts:
                nhb=0
                clust_pop=clust.split('_')[-1].split('.')[0]
                nclust_pop=int(clust_pop)
                p = Pose()
                generate_nonstandard_residue_set(p,ligparams)
                pose_from_file(p, clust)
                res_scores=[]
                #
                ligand_resnum=p.total_residue()
                for res_ind in range(1,ligand_resnum):
                    contact_pose=Pose()
                    ligand_pose=p.residue(ligand_resnum).clone()
                    res_pose=p.residue(res_ind).clone()
                    contact_pose.append_residue_by_jump(res_pose, 1)
                    contact_pose.append_residue_by_jump(ligand_pose, 1)
                    res_resname=contact_pose.residue(1).name()[0:3]
                    lig_resname=contact_pose.residue(2).name()[0:3]
                    if res_resname==lig_resname:
                        print('uh oh')
                        continue
                    else:
                        pass
                    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
                    contact_pose.update_residue_neighbors()
                    rosetta.core.scoring.hbonds.fill_hbond_set(contact_pose, False, hbond_set)
                    if hbond_set.nhbonds()>0:
                        nhb+=1
                    if nhb>=(0.1*nclust_pop):
                        hbond_clusts.append(os.path.join(resdir,clust))
                        break
            os.chdir('..')


rescounts={}
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
for i in polar_residues:
    c=0
    for j in hbond_clusts:
        rn=j.split('_')[3]
        if rn==i:
            cp=j.split('_')[-1].split('.')[0]
            cpn=int(cp)
            c+=cpn
    rescounts[i]=c

'''
1
{'SER': 41,
 'THR': 31,
 'GLN': 6,
 'ASN': 110,
 'TYR': 104,
 'TRP': 11,
 'GLY': 0,
 'HIS': 101,
 'ARG': 158,
 'LYS': 78,
 'GLU': 160,
 'ASP': 49}

 2
{'SER': 1582,
 'THR': 1012,
 'GLN': 370,
 'ASN': 1259,
 'TYR': 420,
 'TRP': 361,
 'GLY': 0,
 'HIS': 558,
 'ARG': 1305,
 'LYS': 1555,
 'GLU': 585,
 'ASP': 1172}
 'esl_Fragment_2_SER_contact_fuzzball_clustered/167_SER_cluster_47.pdb',
 'esl_Fragment_2_LYS_contact_fuzzball_clustered/258_LYS_cluster_84.pdb',
 'esl_Fragment_2_LYS_contact_fuzzball_clustered/49_LYS_cluster_82.pdb',

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl/residue_contact_fuzzballs/esl_Fragment_2_residue_contacts/esl_Fragment_2_SER_contact_fuzzball_clustered/1347_SER_cluster_21.pdb ~/desktop/1347_SER_cluster_21.pdb



AYO TRY ONE MORE FRAG CLUSTERING HBOND STATS THING USING THE 2 HYDROXYLS BUT NOT THE WHOLE
PENTANE RING AND MAYBE I CAN GET INTERESTING RESULTS?
'''
import matplotlib.pyplot as plt
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','HIS','ARG','LYS','GLU','ASP'] #no his

# d=...
tc=0
for key in d.keys():
    tc+=d[key]
freqs=[]
for i in polar_residues:
    try:
        freqs.append((i,d[i]/float(tc)))
    except:
        pass
freqs=sorted(freqs, reverse=True, key=lambda nmem: nmem[1])
labs=[]
vals=[]
for a,b in freqs:
    labs.append(a)
    vals.append(b)

fig = plt.figure()
ax = fig.add_subplot()
ax.bar(labs,vals,color='orange',width=0.1)
ax.set_title('Hbond Residue Preferences')
ax.set_ylabel('Frequency')
ax.set_xlabel('Residue Type')
# plt.savefig('quinonerespref.pdf')
plt.savefig('secohrespref.pdf')
plt.clf()
plt.close()











'''
lemme look at og cluster 3 (benzene) and see if i can plot score vs clust freq
REDO WITH JUST BENZEE FRAG SO NO CLASEHES
'''
#
echo "
import os
dirs=[i for i in os.listdir() if os.path.isdir(i)==True]
#
from pyrosetta import *
init('-pack_missing_sidechains False')
#
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep,fa_elec#, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
#
ligparams=['/wynton/home/kortemme/cgalvin/esl2/Inputs/Rosetta_Inputs/esl.params']
#
clusts_data={}
#
for resdir in dirs:
    os.chdir(resdir)
    clusts=[i for i in os.listdir() if i[-3:]=='pdb']
    clust_data=[]
    for clust in clusts:
        clust_pop=clust.split('_')[-1].split('.')[0]
        nclust_pop=int(clust_pop)
        p = Pose()
        generate_nonstandard_residue_set(p,ligparams)
        pose_from_file(p, clust)
        res_scores=[]
        #
        ligand_resnum=p.total_residue()
        for res_ind in range(1,ligand_resnum):
            contact_pose=Pose()
            ligand_pose=p.residue(ligand_resnum).clone()
            res_pose=p.residue(res_ind).clone()
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
            res_scores.append(score_this_one)
        bestscore=min(res_scores)
        clust_data.append((nclust_pop,bestscore))
    clusts_data[resdir.split('_')[3]]=clust_data
    os.chdir('..')

import json
json.dump(clusts_data,open('clustsdata.json','w'))
">clustsdata.py

'''
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python clustsdata.py
qstat -j "$JOB_ID"
'>cd.sh
qsub -cwd -l mem_free=10G cd.sh
'''
import json
with open('clustsdata.json','r') as f:
    clusts_data=json.load(f)

#plot score vs freq across all
x=[]
y=[]
for key in clusts_data.keys():
    for i in clusts_data[key]:
        if -5<i[1]<5:
            x.append(i[1])
            y.append(i[0])

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot()
#
ax.scatter(x,y)
#
ax.set_title('Rosetta Energy Score vs. Cluster Population')
ax.set_ylabel('Population')
ax.set_xlabel('Rosetta Energy Score (REU)')
plt.savefig('scorevspop.pdf')
plt.clf()

import scipy
from scipy.stats import pearsonr
scipy.stats.pearsonr(x,y)

'''
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl2/residue_contact_fuzzballs/esl_Fragment_3_residue_contacts/scorevspop.pdf ~/desktop/scorevsreubenzene.pdf
'''
import os
l=[i for i in os.listdir() if i[-3:]=='pdb']
for i in l:
    pop=i.split('_')[-1].split('.')[0]
    n=int(pop)
    if n>30:
        print(i)
'''
333_TRP_cluster_115.pdb
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl2/residue_contact_fuzzballs/esl_Fragment_3_residue_contacts/esl_Fragment_3_TRP_contact_fuzzball_clustered/333_TRP_cluster_115.pdb ~/desktop/benzclust.pdb

207_TYR_cluster_76.pdb
13_TYR_cluster_40.pdb

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl2/residue_contact_fuzzballs/esl_Fragment_3_residue_contacts/esl_Fragment_3_TYR_contact_fuzzball_clustered/207_TYR_cluster_76.pdb ~/desktop/207_TYR_cluster_76.pdb
scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl2/residue_contact_fuzzballs/esl_Fragment_3_residue_contacts/esl_Fragment_3_TYR_contact_fuzzball_clustered/13_TYR_cluster_40.pdb ~/desktop/13_TYR_cluster_40.pdb



I WANNA SHOW STACKED HISTS OF MOAD SUBNANOBINDERS AND ESL DESIGNS

scp cgalvin@log2.wynton.ucsf.edu:/wynton/home/kortemme/cgalvin/esl1/1np3hb/1np/buried/genpot/clean/ed1/filtered/analysis/run/scores.json ~/desktop/enzdesanalysis.json

ddg
lighphobesasa

ddg_nligheavy
lighphobesasa_nheavynp

'''
import json
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

df=open("alldata_custom.json",'r')
d=json.load(df)
df.close()

# df2=open("enzdesanalysis.json",'r')
# d2=json.load(df2)
# df2.close()
#analysis of scores from json file
sfname='enzdesanalysis.json'
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

newscores=[]
for dd in scores:
    try:
        ddg=dd['ddg']
        nddg=ddg/21
        dd['ddg_nligheavy']=nddg
        #
        lhps=dd['lighphobesasa']
        n=lhps/18
        dd['lighphobesasa_nheavynp']=n
        newscores.append(dd)
    except:
        pass

scores=newscores
terms=list(scores[0].keys())
print(len(bad))
d2=scores

relterms1=[]
for a in d:
    for key in list(a.keys()):
        relterms1.append(key)
relterms1=list(set(relterms1))


relterms=[i for i in relterms1 if i in terms]
#


pdfname='stackedhistssss.pdf'
pdf = PdfPages(pdfname)
#
# relterms=['lighphobesasa_nheavynp']
#gonna do just for sub micromolar vs micromolar and up binders
for term in relterms:
    subx=[];x=[];morex=[]
    suby=[];y=[];morey=[]
    try:
        for dict in d:
            if float(dict['KD_nm'])<1000:
                suby.append(float(dict['KD_nm']))
                subx.append(float(dict[term]))
        for dict2 in d2:
            morex.append(float(dict2[term]))
    except:
        pass
    try:
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.hist(subx,color='red',alpha=0.3,density=True)
        # ax.hist(x,color='green',alpha=0.3)
        ax.hist(morex,color='blue',alpha=0.3,density=True)
        ax.axvline(np.mean(subx),c='r',linestyle='dashed')
        # ax.axvline(np.mean(x),c='g',linestyle='dashed')
        ax.axvline(np.mean(morex),c='b',linestyle='dashed')
        ax.set_xlabel(term)
        ax.set_ylabel('Count')
        ax.set_title(term)
        ax.text(0.9,0.9,'High Affinity Binders = '+str(np.mean(subx))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
        # ax.text(0.9,0.8,'nm mean = '+str(np.mean(x))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
        ax.text(0.9,0.7,'Designs = '+str(np.mean(morex))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
        pdf.savefig()
        plt.clf()
        plt.close()
    except:
        pass
pdf.close()



fig = plt.figure()
ax = fig.add_subplot()
ax.hist(subx,color='red',alpha=0.3,density=True)
# ax.hist(x,color='green',alpha=0.3)
ax.hist(morex,color='blue',alpha=0.3,density=True)
ax.axvline(np.mean(subx),c='r',linestyle='dashed')
# ax.axvline(np.mean(x),c='g',linestyle='dashed')
ax.axvline(np.mean(morex),c='b',linestyle='dashed')
ax.set_xlabel(term)
ax.set_ylabel('Count')
ax.set_title(term)
ax.text(0.9,0.9,'High Affinity Binders = '+str(np.mean(subx))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='red', fontsize=8)
# ax.text(0.9,0.8,'nm mean = '+str(np.mean(x))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='green', fontsize=8)
ax.text(0.9,0.7,'Designs = '+str(np.mean(morex))[0:6],verticalalignment='top', horizontalalignment='right',transform=ax.transAxes,color='blue', fontsize=8)
plt.savefig('lhps.pdf')
plt.clf()
plt.close()
