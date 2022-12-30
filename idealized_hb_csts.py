#time python ~/bsff_scripts/enumerate_polar_motifs_liberal.py a8s_polar_contacts_clustered.json 20 a8s 3 ndesiredmotifs dist_cutoff_index
'''
dist_ind[5]=2.4 (allows d up to 2.5 ang)

what about bins that fall within what rosetta considers a hbond?

sp2 acc
chi 0-360
psi 70-180
phi 100-180

sp3
no chi
psi 80-180
phi 100-180

psi/psi 80 = ang_ind[4]

for ser/thr with sp3 acc:
    dist bin 3

ya you know what, write a certain number of idealized cst blocks for each
restype:ligatomtype interaction and fucking match those dawg, based off of the
statistics of these top clusters

HOW ABOUT IF I TRY FOR SER/THR TO SPECIFY THEIR ENTIRE RANGE OF FAVORABLE
VALS AT LEAST FOR SC MEDIATED HBONDS

scp -r cgalvin@log2.wynton.ucsf.edu:round6/a8s/rand_motifs/ex ~/desktop/ex
    ya know its not all actually so bad,
    its more like i just got too few matches,
    and of those i did get the ser/thr hbonds have too low a success rate for being
    considered an hbonbd in chimera viewer (which looks like correct assessment though)

a lotta backbone contacts have a bb n accepting an hbond, which aint happenin g,
delete all of those - have N as first contact atom
additionally, perhaps i will exclusively allow glycine to have bb contacts,
for all other res erase if first res contact atom is either N or O
        I ALSO JUST DIDNT MATCH THE 3P1NP/4P - simply do this
'''


import json
import sys
import os
import random
'''
frags_one_contact_only=['Fragment_1','Fragment_2']
clusterjson='a8s_polar_contacts_clustered_liberal.json'
ligand_resname='a8s'
resultdir_prefix='newpm_test'
ndesired_motifs=100
dist_cutoff_index=5
'''
'''
frags_one_contact_only=['Fragment_2','Fragment_3']
clusterjson='38e_polar_contacts_clustered_liberal.json'
ligand_resname='38e'
resultdir_prefix='38e'
ndesired_motifs=100
dist_cutoff_index=5
'''
#first load the clustering results
# clusterjson=sys.argv[1]
df=open(clusterjson,'r')
clusterdf=json.load(df)
df.close()
#
# n_clust_each_frag=int(sys.argv[2])
# ligand_resname=sys.argv[3]
# nmotifres=int(sys.argv[4])
# resultdir_prefix=sys.argv[5]
# ndesired_motifs=sys.argv[6]
# dist_cutoff_index=int(sys.argv[7])




#load the list of bin values for each parameter
#
dist_ind=[i/10 for i in range(14,36,2)]
ang_ind=[i for i in range(0,181,20)]
tor_ind=[i for i in range(-180,181,20)]
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

#
# motif_clust_indices=[]
# #now how can I properly enumerate the indices of the motifs that I want
# if nmotifres==2:
#     for x in range(n_clust_each_frag):
#         for y in range(n_clust_each_frag):
#             motif_clust_indices.append((x,y))
# if nmotifres==3:
#     for x in range(n_clust_each_frag):
#         for y in range(n_clust_each_frag):
#             for z in range(n_clust_each_frag):
#                 motif_clust_indices.append((x,y,z))
# elif nmotifres==4:
#     for x in range(n_clust_each_frag):
#         for y in range(n_clust_each_frag):
#             for z in range(n_clust_each_frag):
#                 for j in range(n_clust_each_frag):
#                     motif_clust_indices.append((x,y,z,j))

#lists of middle values of bins for each parameter
#will use this as input cst value with tolerance aligning with bin range
#ie ang bins = 10 degrees, ex val = 175, tol = +/-5
dist_ind=[i/100 for i in range(140,350,20)]
ang_ind=[i for i in range(0,180,20)]
tor_ind=[i for i in range(-180,180,20)]
residue_tag = 'residue3'
torsion_constraint_sample_number=5  #IDK WHAT IDEAL VALUE IS FOR THIS, TRY 8 AND SEE WHAT HAPPENS BUT DONT FORGET TO PLAY AROUND WITH THIS

#move all overrep clusters to a new df
'''
hi->low freq in moad bs:
    ser,tyr,thr,asp,asn,trp,arg,his,glu,lys,gln

would be nice if i could also filter contacts based on whether rosetta consider it
a hb ... 4 which i just need bins rosetta consider good
'''
#only ser,thr,gly spared from classification as problematic
problematic_residues=['ARG','LYS','ASN','GLN','TRP','HIS','GLU','ASP','TYR']
charged_residues=['ARG','LYS','GLU','ASP']
overrep_df={}
for frag in clusterdf.keys():
    frag_contacts=clusterdf[frag]
    #
    top_5th_percentile_index=int(0.05*(len(frag_contacts)))
    top_5th_percentile_index_val=frag_contacts[top_5th_percentile_index][0]
    #
    top_3th_percentile_index=int(0.01*(len(frag_contacts)))
    top_3th_percentile_index_val=frag_contacts[top_3th_percentile_index][0]
    #
    top_5th_percentile_contacts=[]
    for r in problematic_residues:
        rc=0
        for j in frag_contacts:
            resname=j[1].split('_')[6]
            if resname==r:
                if rc<5:
                    if j[0]>=top_3th_percentile_index_val:
                        top_5th_percentile_contacts.append(j)
                        rc+=1
    for j in frag_contacts:
        resname=j[1].split('_')[6]
        if resname not in problematic_residues:
            if j[0]>=top_5th_percentile_index_val:
                top_5th_percentile_contacts.append(j)
        else:
            pass
    filt_contacts=[]
    for fc in top_5th_percentile_contacts:
        first_atom=fc[1].split('_')[0]
        dist_i=int(fc[1].split('_')[7])
        resname=fc[1].split('_')[6]
        if first_atom!='N':
            if resname!='GLY':
                if first_atom!='O':
                    if first_atom!='H':
                        if dist_i<=dist_cutoff_index:
                            filt_contacts.append(fc)
            else:
                if dist_i<=dist_cutoff_index:
                    filt_contacts.append(fc)
    print(frag)
    print(len(filt_contacts))
    overrep_df[frag]=filt_contacts

'''
okay now for each fragment:
    for each restype:
        for each param:
            get mean, std of param
'''
restypes=['ARG','LYS','ASN','GLN','TRP','HIS','GLU','ASP','TYR','GLY','SER','THR']

for fragid in overrep_df.keys():
    for res in restypes:
        resdists=[]
        resphis=[]
        respsis=[]
        reschis=[]
        resdistsa=[]
        resphisa=[]
        respsisa=[]
        reschisa=[]
        for contact in overrep_df[fragid]:
            cdata=contact[1].split('_')
            contact_restype=cdata[6]
            if contact_restype==res:
                dbin=int(cdata[7])-1
                res_atoms=[cdata[0],cdata[1],cdata[2]]
                if 'H' in res_atoms[0]:
                    resdonor=True
                else:
                    resdonor=False
                if resdonor==True:
                    # lig_atoms=[contact_data[3],contact_data[4],contact_data[5]]
                    phibin=int(cdata[8])-1
                    psibin=int(cdata[9])-1
                    try:
                        chibin=int(contact_data[10])-1
                    except:
                        pass
                    resdists.append(dbin)
                    resphis.append(phibin)
                    respsis.append(psibin)
                    try:
                        reschis.append(chibin)
                    except:
                        pass
                elif resdonor==False:
                    # lig_atoms=[contact_data[3],contact_data[4],contact_data[5]]
                    phibin=int(cdata[8])-1
                    psibin=int(cdata[9])-1
                    try:
                        chibin=int(cdata[10])-1
                    except:
                        pass
                    resdistsa.append(dbin)
                    resphisa.append(phibin)
                    respsisa.append(psibin)
                    try:
                        reschisa.append(chibin)
                    except:
                        pass
        dm=np.mean(resdists)
        phim=np.mean(resphis)
        psim=np.mean(respsis)




distanceAB=dist_ind[dbin]
distance_tolerance=0.1
distance_constraint_sample_number=0
if resdonor==True:
    angle_A=ang_ind[psibin]
    angle_B=ang_ind[phibin]
    try:
        torsion_A=tor_ind[int(contact_data[10])-1]
        torsion_A_tolerance=10
        torsion_A_constraint_sample_number=1
    except:
        torsion_A=0
        torsion_A_tolerance=180
        torsion_A_constraint_sample_number=torsion_constraint_sample_number
    torsion_AB=0
    torsion_B=0
    angle_A_tolerance=10
    angle_B_tolerance=10
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
        torsion_B_tolerance=10
        torsion_B_constraint_sample_number=1
    except:
        torsion_B=0
        torsion_B_tolerance=180
        torsion_B_constraint_sample_number=torsion_constraint_sample_number
    torsion_AB=0
    torsion_A=0
    angle_A_tolerance=10
    angle_B_tolerance=10
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
m3rf.append(cstblock)
# of.write('CST::BEGIN\n')
# of.write('\n'.join(cstblock))
# of.write('\n')
# of.write('CST::END\n')
# of.write('\n')
ofs.append(m3rf)

#4 res
