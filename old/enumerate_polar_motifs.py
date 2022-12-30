#time python enumerate_polar_motifs.py clusterresults.json n_clust_each_frag lig_resname nmotifres prefix
#time python ~/bsff_scripts/enumerate_polar_motifs.py a8s_polar_contacts_clustered.json 20 a8s 3


import json
import sys
import os

#first load the clustering results
clusterjson=sys.argv[1]
df=open(clusterjson,'r')
clusterdf=json.load(df)
df.close()

n_clust_each_frag=int(sys.argv[2])

ligand_resname=sys.argv[3]

nmotifres=int(sys.argv[4])

resultdir_prefix=sys.argv[5]
'''
#dev
df=open('a8s_polar_contacts_clustered.json','r')
clusterdf=json.load(df)
df.close()
n_clust_each_frag=5
ligand_resname='a8s'
'''
#load the list of bin values for each parameter
#
dist_ind=[i/10 for i in range(14,36,1)]
ang_ind=[i for i in range(0,190,10)]
tor_ind=[i for i in range(-180,190,10)]
#
dist_ind_pairs=[]
for i,e in enumerate(dist_ind[:-1]):
    upperlim=dist_ind[i+1]
    dist_ind_pairs.append((e,upperlim))
ang_ind_pairs=[]
for i,e in enumerate(ang_ind[:-1]):
    upperlim=ang_ind[i+1]
    ang_ind_pairs.append((e,upperlim))
tor_ind_pairs=[]
for i,e in enumerate(tor_ind[:-1]):
    upperlim=tor_ind[i+1]
    tor_ind_pairs.append((e,upperlim))


def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.05, angle_A_tolerance=5, angle_B_tolerance=5,
                                          torsion_A_tolerance=5, torsion_AB_tolerance=5, torsion_B_tolerance=5,
                                          torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                          distance_constraint_sample_number=0, greasy_sampling=True):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue3'
    # Increase tolerance/sampling by 5/1 for greasy residues (ACFILMVWY)
    if greasy_sampling and residue_resname in ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']:
        angle_A_tolerance += 10
        angle_B_tolerance += 10
        torsion_A_tolerance += 10
        torsion_AB_tolerance += 10
        torsion_B_tolerance += 10
        torsion_constraint_sample_number += 1
        angle_constraint_sample_number += 1

    cstblock = [
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


motif_clust_indices=[]
#now how can I properly enumerate the indices of the motifs that I want
if nmotifres==2:
    for x in range(n_clust_each_frag):
        for y in range(n_clust_each_frag):
            motif_clust_indices.append((x,y))
if nmotifres==3:
    for x in range(n_clust_each_frag):
        for y in range(n_clust_each_frag):
            for z in range(n_clust_each_frag):
                motif_clust_indices.append((x,y,z))
elif nmotifres==4:
    for x in range(n_clust_each_frag):
        for y in range(n_clust_each_frag):
            for z in range(n_clust_each_frag):
                for j in range(n_clust_each_frag):
                    motif_clust_indices.append((x,y,z,j))

#lists of middle values of bins for each parameter
#will use this as input cst value with tolerance aligning with bin range
#ie ang bins = 10 degrees, ex val = 175, tol = +/-5
dist_ind=[i/100 for i in range(145,350,10)]
ang_ind=[i for i in range(5,180,10)]
tor_ind=[i for i in range(-175,180,10)]
residue_tag = 'residue3'
torsion_constraint_sample_number=8  #IDK WHAT IDEAL VALUE IS FOR THIS, TRY 8 AND SEE WHAT HAPPENS BUT DONT FORGET TO PLAY AROUND WITH THIS
#okay now writing motif csts
count=1
os.makedirs(resultdir_prefix+'_polar_csts',exist_ok=True)
for indset in motif_clust_indices:
    motifresnames=[]
    cstname=os.path.join(resultdir_prefix+'_polar_csts',str(count)+'.cst')
    of=open(cstname,'w')
    for frag_ind in range(len(indset)):
        contact=clusterdf[list(clusterdf.keys())[frag_ind]][indset[frag_ind]][1]
        contact_data=contact.split('_')
        res_atoms=[contact_data[0],contact_data[1],contact_data[2]]
        if 'H' in res_atoms[0]:
            resdonor=True
        else:
            resdonor=False
        lig_atoms=[contact_data[3],contact_data[4],contact_data[5]]
        res_resname=contact_data[6]
        motifresnames.append(res_resname)
        dbin=int(contact_data[7])-1
        phibin=int(contact_data[8])-1
        psibin=int(contact_data[9])-1
        try:
            chibin=int(contact_data[10])-1
        except:
            pass
        distanceAB=dist_ind[dbin]
        distance_tolerance=0.1
        distance_constraint_sample_number=0
        if resdonor==True:
            angle_A=ang_ind[psibin]
            angle_B=ang_ind[phibin]
            try:
                torsion_A=tor_ind[int(contact_data[10])-1]
                torsion_A_tolerance=5
                torsion_A_constraint_sample_number=1
            except:
                torsion_A=0
                torsion_A_tolerance=180
                torsion_A_constraint_sample_number=torsion_constraint_sample_number
            torsion_AB=0
            torsion_B=0
            angle_A_tolerance=5
            angle_B_tolerance=5
            torsion_AB_tolerance=180
            torsion_B_tolerance=180
            cstblock = [
                '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(lig_atoms)),
                '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
                '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(res_atoms)),
                '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, res_resname),
                '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                    distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
                '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_A), angle_A_tolerance, 100, 1),
                '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_B), angle_B_tolerance, 100, 1),
                '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_A), torsion_A_tolerance, 100,
                    torsion_A_constraint_sample_number),
                '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_B), torsion_B_tolerance, 100,
                    torsion_constraint_sample_number),
                '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_AB), torsion_AB_tolerance, 100,
                    torsion_constraint_sample_number)]
        elif resdonor==False:
            angle_A=ang_ind[phibin]
            angle_B=ang_ind[psibin]
            try:
                torsion_B=tor_ind[int(contact_data[10])-1]
                torsion_B_tolerance=5
                torsion_B_constraint_sample_number=1
            except:
                torsion_B=0
                torsion_B_tolerance=180
                torsion_B_constraint_sample_number=torsion_constraint_sample_number
            torsion_AB=0
            torsion_A=0
            angle_A_tolerance=5
            angle_B_tolerance=5
            torsion_AB_tolerance=180
            torsion_A_tolerance=180
            cstblock = [
                '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(lig_atoms)),
                '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
                '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(res_atoms)),
                '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, res_resname),
                '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                    distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
                '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_A), angle_A_tolerance, 100, 1),
                '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(angle_B), angle_B_tolerance, 100, 1),
                '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_A), torsion_A_tolerance, 100,
                    torsion_constraint_sample_number),
                '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_B), torsion_B_tolerance, 100,
                    torsion_B_constraint_sample_number),
                '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                    float(torsion_AB), torsion_AB_tolerance, 100,
                    torsion_constraint_sample_number)]
        of.write('CST::BEGIN\n')
        of.write('\n'.join(cstblock))
        of.write('\n')
        of.write('CST::END\n')
        of.write('\n')
    narg=0
    for tn in motifresnames:
        if tn=='ARG':
            narg+=1
    if narg>1:
        of.close()
        os.system('rm '+cstname)
        continue
    elif len(list(set(motifresnames)))<2:
        of.close()
        os.system('rm '+cstname)
        continue
    else:
        of.close()
        count+=1
