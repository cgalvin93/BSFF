time python ~/desktop/BSFF/new_target.py dog

'''
AVOGADRO
    hydrogens for ph 7.4
    optimize energy with mm94 ff
    Extensions --> Molecular Mechanics --> Setup Force Field...
        mm94, default
    Extensions --> Molecular Mechanics --> Conformer Search...
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
    1-carbonyl (no ether o)
    2-alcohol
    3-alcohol
    4-alcohol

compounds=("dog")
for I in ${compounds[@]}; do
    cd ${I}
    time python ~/desktop/BSFF/Fragment_Search.py ${I}_frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
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


okay now it time to move the target folders over to the cluster
and we will process the search results so that we can split across many nodes

scp -r dog cgalvin@log2.wynton.ucsf.edu:

submit the alignment jobs
dog
2888
qsub -cwd -t 1-2888 -l mem_free=1G -o cluster_output -e cluster_output run_align.sh


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

dog
2704
chmod ugo+x run_filter_array.sh
qsub -cwd -t 1-2704 -l mem_free=2G -o cluster_output -e cluster_output run_filter_array.sh



CONSOLIDATE FILTERED CONTACTS
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python ~/BSFF/consolidate_filtered_contacts.py
'>run_consoldate.sh
qsub -cwd -l mem_free=4G -o cluster_output -e cluster_output run_consoldate.sh





PROCESS FILTERED CONTACTS - do for each fragment
compounds=("Fragment_1" "Fragment_2" "Fragment_3" "Fragment_4")
for I in ${compounds[@]}; do
    time python3 ~/BSFF/process_filtered_contacts2.py /wynton/home/kortemme/cgalvin/dog/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb dog $I
done

mkdir fragfbs
mv *_residue_contacts fragfbs





CLUSTER POLAR
gotta add dictionary to top of script specifying relevant hbond atoms
for each fragment

scp cgalvin@log2.wynton.ucsf.edu:dog/Inputs/Rosetta_Inputs/dog_0001.pdb ~/desktop/dog_0001.pdb
relevant_atoms={'Fragment_1':['O5'],
                'Fragment_2':['O1','H1'],
                'Fragment_3':['O2','H2'],
                'Fragment_4':['O3','H3']}
scp cpdog.py cgalvin@log2.wynton.ucsf.edu:BSFF/

time python3 ~/BSFF/cpdog.py dog /wynton/home/kortemme/cgalvin/dog/fragfbs /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb Fragment_1 Fragment_2 Fragment_3 Fragment_4

compound_dir='dog'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/cpdog.py dog /wynton/home/kortemme/cgalvin/dog/fragfbs /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb Fragment_1 Fragment_2 Fragment_3 Fragment_4
qstat -j "$JOB_ID"
'>cp_$compound_dir.sh
chmod ugo+x cp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cp_$compound_dir.sh


TO GET STATS
compound_dir='dog'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/new_cpdog.py dog /wynton/home/kortemme/cgalvin/dog/fragfbs /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog.params /wynton/home/kortemme/cgalvin/dog/Inputs/Rosetta_Inputs/dog_0001.pdb Fragment_1 Fragment_2 Fragment_3 Fragment_4
qstat -j "$JOB_ID"
'>nscp_$compound_dir.sh
chmod ugo+x nscp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output nscp_$compound_dir.sh











CLUSTER NONPOLAR - dont have any np frags for this rn but ill leave code here since
its usually a part of workflow
compound_dir='38e'
echo '
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pyr37
time python3 ~/BSFF/cluster_nonpolar.py 38e /wynton/home/kortemme/cgalvin/round6/38e/38efragfbs /wynton/home/kortemme/cgalvin/round6/38e/Inputs/Rosetta_Inputs/38e_0001.pdb Fragment_4 Fragment_5 Fragment_6
qstat -j "$JOB_ID"
'>cnp_$compound_dir.sh
chmod ugo+x cnp_$compound_dir.sh
qsub -cwd -l mem_free=16G -o cluster_output -e cluster_output cnp_$compound_dir.sh



GENERATE POLAR MOTIF CSTS - two on carb, 1 on all the OHs = 5 residue motif
polar_motifs_rand.py
'''
