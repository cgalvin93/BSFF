#time ipython ~/desktop/prj/bsff/bsff_scripts/fuzzball_for_assembly.py BSFF_COMPOUND_PARENT_DIRECTORY(EG DEX)
#parent directory should have 'fixed_transformed_aligned_pdbs' directory in it
#and each fragment directory should contain another directory with the cluster results, called
#clustered_contacts (default from cluster_contacts.py)
#time ipython ~/desktop/prj/bsff/bsff_scripts/fuzzball_for_assembly.py a8s frag1 frag2 frag3 frag4

import os
import sys
import json
import math
from pyrosetta import *
init('-load_PDB_components false')

#LOAD SCOREFUNCTION
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc#,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(fa_sol, 1)
sf.set_weight(hbond_sc, 0.)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,0.)
# sf.set_weight(lk_ball_iso,-0.38)
scoreterms=['total_score',
            'fa_elec',
            'fa_sol',
            'fa_rep',
            'fa_atr',
            'hbond_bb_sc',
            'hbond_sc',
            # 'lk_ball_iso',
            'score_weighted_freq']

stat_weight=1

#FRAGMENT CLUSTER DIRECTORIES AND PARAMS
# TAPDBS_path='a8s'+'/Transformed_Aligned_PDBs/'
# params_directory_path='a8s'+'/Inputs/Rosetta_Inputs/'
TAPDBS_path=sys.argv[1]+'/Transformed_Aligned_PDBs/'
params_directory_path=sys.argv[1]+'/Inputs/Rosetta_Inputs/'
params_file_path=[]
for i in os.listdir(params_directory_path):
    if i[-7:]=='.params':
        params_file_path.append(os.path.join(params_directory_path,i))

#
relfrags=[]
try:
    f1=sys.argv[2]
    relfrags.append(f1)
except:
    pass
try:
    f2=sys.argv[3]
    relfrags.append(f2)
except:
    pass
try:
    f3=sys.argv[4]
    relfrags.append(f3)
except:
    pass
try:
    f4=sys.argv[4]
    relfrags.append(f4)
except:
    pass
#
cluster_dir_paths=[]
for frag in relfrags:
    cluster_path=os.path.join(TAPDBS_path,frag,'clustered_contacts')
    cluster_dir_paths.append(cluster_path)


###################
fuzzball_pose=Pose()
n_all_residues=0
n_fuzzball_residues=0
fuzzball_pose_frequencies={}
#go through clusters and score
npres=['LEU','ILE','VAL','PHE','TYR','TRP','MET']
for cluster_dir in cluster_dir_paths: ############################
    try:
        rawpdbfiles=[i for i in os.listdir(cluster_dir) if i[-3:]=='pdb']
        pdbfiles=[]         #filter to only np res
        for f in rawpdbfiles:
            s=f.split('_')[1]
            if s in npres:
                pdbfiles.append(f)
    except:
        continue
    total_n_contacts_for_fragment=[]
    for file in pdbfiles:########################################################################## #first get the total n contacts so i can norm for each cluster and use to weight scores
        try:
            s=file.split('_')[3]
            s2=s.split('.')[0]
            n_contacts_current_cluster=int(s2)-1
            if n_contacts_current_cluster>=1:
                total_n_contacts_for_fragment.append(n_contacts_current_cluster)
        except:
            print('error with: '+str(os.path.join(cluster_dir,file)))
    n_all_contacts=sum(total_n_contacts_for_fragment)
    n_clusters_for_fragment=len(total_n_contacts_for_fragment)
    ref_state=(1./n_clusters_for_fragment)
    for file in pdbfiles: #########################################################################now go through and actually score
        try:
            p = Pose()
            generate_nonstandard_residue_set(p,params_file_path)
            pose_from_file(p, os.path.join(cluster_dir,file))
            ligand_resnum=p.total_residue()
            if p.residue(ligand_resnum).name()!=sys.argv[1]:##########################################################
                continue
            else:
                n_contacts_current_cluster=p.total_residue()-1
                cluster_frequency=n_contacts_current_cluster/n_all_contacts
                for i in range(n_contacts_current_cluster): #HERE ITERATING THROUGH RES IN CLUSTER
                    resnum=i+1
                    n_all_residues+=1
                    contact_pose=Pose()
                    ligand_pose=p.residue(p.total_residue()).clone()
                    res_pose=p.residue(resnum).clone()
                    contact_pose.append_residue_by_jump(res_pose, 1)
                    contact_pose.append_residue_by_jump(ligand_pose, 1)
                    rosetta.core.pack.optimizeH(contact_pose, sf)
                    sf(contact_pose)
                    w=contact_pose.energies().energy_graph().find_energy_edge(1,2)
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    hbbbsc=w[rosetta.core.scoring.hbond_bb_sc]
                    faatr=w[rosetta.core.scoring.fa_atr]
                    farep=w[rosetta.core.scoring.fa_rep]
                    fasol=w[rosetta.core.scoring.fa_sol]
                    faelec=w[rosetta.core.scoring.fa_elec]
                    # lkballiso=w[rosetta.core.scoring.lk_ball_iso]
                    total_score=faelec+fasol+farep+faatr+hbbbsc+hbsc#+lkballiso
                    weighted_score=total_score-(stat_weight*(math.log(cluster_frequency/ref_state)))
                    if weighted_score<=0.0:
                        fuzzball_pose.append_residue_by_jump(res_pose, 1)
                        n_fuzzball_residues+=1
                        fuzzball_pose_frequencies[n_fuzzball_residues]=(cluster_frequency,ref_state)
        except:
            print('error with: '+str(os.path.join(cluster_dir,file)))


fuzzball_pose.append_residue_by_jump(ligand_pose, 1)
os.chdir(sys.argv[1])
# os.chdir('a8s')

fuzzball_pose.dump_pdb('fuzzball_for_assembly.pdb')
json.dump(fuzzball_pose_frequencies,open('np_fuzzball_for_assembly_frequencies.json','w'))
print('\n')
print('\n')
print('\n')
print('\n')

print('THE NUMBER OF NP FUZZBALL CONTACTS IS: '+str(n_fuzzball_residues))
print('OUT OF A TOTAL OF '+str(n_all_residues)+' CONTACTS PARSED')
qstats=open('Fuzzball_filter_stats.txt','w')
qstats.write('THE NUMBER OF FUZZBALL CONTACTS IS: '+str(n_fuzzball_residues)+'\n')
qstats.write('OUT OF A TOTAL OF '+str(n_all_residues)+' CONTACTS PARSED')
qstats.close()

# f=open('fuzzball_for_assembly.pdb','r')
# lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[0:3]=='TER' or line[0:6]=='HETNAM' or line[0:6]=='HETATM']
# f.close()
# ofile=open('clean_fuzzball_for_assembly.pdb','w')
# cleanlines=[]
# for line in lines:
#     if line[0:4]=='ATOM':
#         cleanlines.append(line)
#     if line[0:6]=='HETATM':
#         if line[17:20]=='HOH':
#             pass
#         else:
#             newline='ATOM  '+line[6:]
#             cleanlines.append(line)
#     if line[:3]=='TER':
#         cleanlines.append(line)
# for i,line in enumerate(cleanlines):
#     try:
#         resnum=int(line[22:29].strip())
#         nextresnum=int(cleanlines[i+1][22:29].strip())
#         if resnum==nextresnum:
#             ofile.write(line)
#         else:
#             ofile.write(line)
#             ofile.write('TER\n')
#     except:
#         ofile.write(line)
# ofile.close()
#
# os.system('rm fuzzball_for_assembly.pdb')
#
# os.chdir('..')


'''
d=json.load(open('fuzzball_for_assembly_frequencies.json','r'))
'''

#now renumber fuzzball res
f=open('fuzzball_for_assembly.pdb','r')
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
og_n_res=str(len(fuzzball_residue_indices))

clean_fuzzball_lines=[]
#reorder the residue numbers of the filtered contact residues because
new_residue_number=1
for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in fuzzball_lines[a:b]]
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
        elif 9999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+''+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1

ofile=open('clean_fuzzball_for_assembly_np.pdb','w')
cleanlines=[]
for line in clean_fuzzball_lines:
    if line[0:4]=='ATOM':
        cleanlines.append(line)
    if line[0:6]=='HETATM':
        if line[17:20]=='HOH':
            pass
        else:
            cleanlines.append(line)
    if line[:3]=='TER':
        cleanlines.append(line)
for i,line in enumerate(cleanlines):
    try:
        resnum=int(line[22:29].strip())
        nextresnum=int(cleanlines[i+1][22:29].strip())
        if resnum==nextresnum:
            ofile.write(line)
        else:
            ofile.write(line)
            ofile.write('TER\n')
    except:
        ofile.write(line)
ofile.close()

os.system('rm fuzzball_for_assembly.pdb')

os.chdir('..')
