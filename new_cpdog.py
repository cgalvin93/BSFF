'''
how about for each fragment
    for each residue type
        i get the corrected dists for each param
            i find the maxima of the distributions and use these vals for
            csts with some tolerance that's like 1 std out or something
'''

#guess ill have to add this manually for every lig....
relevant_atoms={'Fragment_1':['O5'],
                'Fragment_2':['O1','H1'],
                'Fragment_3':['O2','H2'],
                'Fragment_4':['O3','H3']}

#
#
import os
import sys
from scipy import stats
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pyrosetta import *
init('-ignore_unrecognized_res -load_PDB_components False -pack_missing_sidechains False -ignore_zero_occupancy False')
import math
import json
import random
from collections import defaultdict

#
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

#
#
parent_directory=sys.argv[1]
fragfbdir=sys.argv[2]
paramspath=sys.argv[3]
ligand_path=sys.argv[4]
f1=sys.argv[5]

try:
    f2=sys.argv[6]
except:
    pass
try:
    f3=sys.argv[7]
except:
    pass
try:
    f4=sys.argv[8]
except:
    pass

# pdfname=str(parent_directory)+'_polar_clusters_plots_liberal.pdf'
# pdf = PdfPages(pdfname)
'''
DEV
parent_directory='a8s'
f1='Fragment_1' #carboxylate
f2='Fragment_2' #cyclohexanone
f3='Fragment_4' #alcohol w some carbons
'''

#
fraglist=[f1]
try:
    fraglist.append(f2)
except:
    pass
try:
    fraglist.append(f3)
except:
    pass
try:
    fraglist.append(f4)
except:
    pass


#
polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY'] #ONLY TRUE POLAR AND GLY
# polar_residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY','ALA','LEU','ILE','VAL','PHE'] #ALL EXCEPT PRO,CYS,MET
fragdata={}
for frag in fraglist:
    all={}
    pol={}
    polfiles=[]
    nall=0
    npol=0
    nppops=[]
    polpops=[]
    clusters_path=os.path.join(fragfbdir,('_').join([parent_directory,frag,'residue_contacts']))
    clusts=[i for i in os.listdir(clusters_path) if i[-3:]=='pdb']
    for clust in clusts:
        s=clust.split('.')[0]
        ss=s.split('_')
        restype=ss[3]
        if restype in polar_residues:
            polfiles.append(clust)
        else:
            pass
    fragdata[frag]=polfiles


params = [paramspath] #ligand params
ligp = Pose()
generate_nonstandard_residue_set(ligp,params)
pose_from_file(ligp, ligand_path)
lig_resname=ligp.residue(1).name()[:3]
print('\n\n'+lig_resname)

#fragdata[frag]=polfiles (path to polar cluster)
contacts_allfrags={}
for frag in fragdata.keys():
    contacts_thisfrag={}
    for clust in fragdata[frag]: ############################################################################################################################
        params = [paramspath] #ligand params
        p = Pose()
        generate_nonstandard_residue_set(p,params)
        pose_from_file(p, os.path.join(fragfbdir,('_').join([parent_directory,frag,'residue_contacts']),clust))
        # chainX = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
        # neighborhood_selector_bool = chainX.apply(p)
        # resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
        # ligand_resnum=list(resnums)[0]
        # # ligand_resnum=p.total_residue()
        # ligand_pose=p.residue(ligand_resnum).clone()
        # lig_resname=ligand_pose.name()[0:3]
        # print('\n\n'+lig_resname)
        cluster_resdon={}
        cluster_resacc={}
        for i in range(p.total_residue()-1): #iterating over every residue in cluster
            resnum=i+1
            #make new pose for ligand and each res
            contact_pose=Pose()
            # ligand_pose=p.residue(ligand_resnum).clone()
            res_pose=p.residue(resnum).clone()
            contact_pose.append_residue_by_jump(res_pose, 1)
            contact_pose.append_residue_by_jump(ligp.residue(1).clone(), 1)
            res_resname=contact_pose.residue(1).name()[0:3]
            lig_resname=contact_pose.residue(2).name()[0:3]
            if res_resname==lig_resname:
                print('lig and res names the same...')
                # print(res_resname)
                # print(lig_resname)
                continue
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
            #remove contacts which do not include hydrogen
            for r in range(10):#not clear to me why but need to repeat to be thorough
                for ix in sd:
                    if 'H' not in ix[0] and 'H' not in ix[1]:
                        sd.remove(ix)
            #here i am trying to make sure the contact is with the right atoms
            sd2=[]
            for ix in sd:
                for atom in relevant_atoms[frag]:
                    if atom in ix[0]:
                        if ix not in sd2:
                            sd2.append(ix)
            sd=sd2
            #
            sd_only_polar=[]
            if len(sd)<1:
                continue
            else:
                pass
            for sdelementidx in range(len(sd)):
                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[sdelementidx]
                if 1.3<=float(distance_AB)<=3.5: #only check further if distance will count for hb
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
                        print('CANT ID LIGATOMSEQUENCE')
                        print(res_resname)
                        print(lig_resname)
                        print('\n\n\n')
                        continue
                    if len(resatomsequence)>0:
                        resbase1=resatomsequence[0][1]
                        resbase2=resatomsequence[0][2]
                        residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
                        residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
                    else:
                        print('CANT ID RESATOMSEQUENCE')
                        print(res_resname)
                        print(lig_resname)
                        print('\n\n\n')
                        continue
                    #fix oxt to O in res if present
                    if residue_atom1=='OXT':
                        residue_atom1='O'
                    if residue_atom2=='OXT':
                        residue_atom2='O'
                    if residue_atom3=='OXT':
                        residue_atom3='O'
                    #make sure first base atom is polar for relevant hydrogen
                    if 'H' in residue_atom1 and ('O' in residue_atom2 or 'N' in residue_atom2) or 'H' in ligand_atom1 and ('O' in ligand_atom2 or 'N' in ligand_atom2):
                        sd_only_polar.append((distance_AB,residue_atom1,residue_atom2,residue_atom3,ligand_atom1,ligand_atom2,ligand_atom3,indexres,resbase1,resbase2,indexlig,ligbase1,ligbase2))
                else:
                    continue
            #here ID closest contact between pol-h and other
            sdop_sorted=sorted(sd_only_polar, key=lambda x: x[0])
            #remove contacts which have hydrogen as first atom for both species
            # print(sdop_sorted)
            for r in range(10):#not clear to me why but need to repeat to be thorough
                for il in sdop_sorted:
                    if 'H' in il[1]:
                        if 'H' in il[4]:
                            sdop_sorted.remove(il)
            #remove contacts where "acceptor" is not polar
            for r in range(10):#not clear to me why but need to repeat to be thorough
                for il in sdop_sorted:
                    if 'H' in il[1]:
                        if 'N' not in il[4]:
                            if 'O' not in il[4]:
                                sdop_sorted.remove(il)
                    elif 'H' in il[4]:
                        if 'N' not in il[1]:
                            if 'O' not in il[1]:
                                sdop_sorted.remove(il)
            # print(sdop_sorted)
            #now if anything is in this list, use the first element to define interaxn geometry
            if len(sdop_sorted)>0:
                distance_AB=sdop_sorted[0][0]
                lig_atoms=[sdop_sorted[0][4],sdop_sorted[0][5],sdop_sorted[0][6]]
                res_atoms=[sdop_sorted[0][1],sdop_sorted[0][2],sdop_sorted[0][3]]
                sresidue_atom1=res_atoms[0]
                sresidue_atom2=res_atoms[1]
                sresidue_atom3=res_atoms[2]
                sligand_atom1=lig_atoms[0]
                sligand_atom2=lig_atoms[1]
                sligand_atom3=lig_atoms[2]
                sindexres=sdop_sorted[0][7]
                sresbase1=sdop_sorted[0][8]
                sresbase2=sdop_sorted[0][9]
                sindexlig=sdop_sorted[0][10]
                sligbase1=sdop_sorted[0][11]
                sligbase2=sdop_sorted[0][12]
                res_atom_coords=[]
                lig_atom_coords=[]
                x,y,z=contact_pose.residue(1).atom(sindexres).xyz()
                res_atom_coords.append((x,y,z))
                x,y,z=contact_pose.residue(1).atom(sresbase1).xyz()
                res_atom_coords.append((x,y,z))
                x,y,z=contact_pose.residue(1).atom(sresbase2).xyz()
                res_atom_coords.append((x,y,z))
                x,y,z=contact_pose.residue(2).atom(sindexlig).xyz()
                lig_atom_coords.append((x,y,z))
                x,y,z=contact_pose.residue(2).atom(sligbase1).xyz()
                lig_atom_coords.append((x,y,z))
                x,y,z=contact_pose.residue(2).atom(sligbase2).xyz()
                lig_atom_coords.append((x,y,z))
                #okay getting angles
                #RES1 IN CONSTRAINT FILE IS ALWAYS LIGAND
                angle_A=getAngle(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
                angle_B=getAngle(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
                #finally, getting dihedrals
                torsion_A=dihedral(lig_atom_coords[2],lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
                torsion_AB=dihedral(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
                torsion_B=dihedral(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1],res_atom_coords[2])
                #
                if 'H' in sresidue_atom1 and ('O' in sresidue_atom2 or 'N' in sresidue_atom2):
                    tcrhs=ligp.residue(1).atom_type(sindexlig)
                    # print(tcrhs)
                    hybridization_state=str(pyrosetta.rosetta.core.chemical.AtomType.hybridization(tcrhs)).split('.')[1].split('_')[0]
                    cluster_resdon[resnum]=[res_resname,lig_resname,distance_AB,angle_A,angle_B,torsion_A,torsion_AB,torsion_B,sresidue_atom1,sresidue_atom2,sresidue_atom3,sligand_atom1,sligand_atom2,sligand_atom3,hybridization_state]#resname,d,psi,phi,chi
                    # print('\n\n\nRESDON')
                    # print(res_resname)
                    # print(residue_atom1)
                    # print(residue_atom2)
                    # print(ligand_atom1)
                    # print(ligand_atom2)
                    # print('\n\n\n')
                elif 'H' in sligand_atom1 and ('O' in sligand_atom2 or 'N' in sligand_atom2):
                    tcrhs=p.residue(resnum).atom_type(sindexres)
                    # print(tcrhs)
                    hybridization_state=str(pyrosetta.rosetta.core.chemical.AtomType.hybridization(tcrhs)).split('.')[1].split('_')[0]
                    cluster_resacc[resnum]=[res_resname,lig_resname,distance_AB,angle_B,angle_A,torsion_B,torsion_A,torsion_AB,sresidue_atom1,sresidue_atom2,sresidue_atom3,sligand_atom1,sligand_atom2,sligand_atom3,hybridization_state]
                    # print('\n\n\nRESACC')
                    # print(res_resname)
                    # print(residue_atom1)
                    # print(residue_atom2)
                    # print(ligand_atom1)
                    # print(ligand_atom2)
                    # print('\n\n\n')
                else:
                    print('\n\n\nNO POTENTIAL HBONDS FOUND FOR THIS CONTACT, CANT ID AS DON OR ACC')
                    print(res_resname)
                    print(lig_resname)
                    print('\n\n\n')
                    continue
            else:
                print('\n\n\nNO POTENTIAL HBONDS FOUND FOR THIS CONTACT:SORTED POL CONTACTS EMPTY')
                print(res_resname)
                print(lig_resname)
                print('\n\n\n')
                continue
        tfd={}
        tfd['resdonors']=cluster_resdon
        tfd['resacc']=cluster_resacc
        contacts_thisfrag[clust]=tfd
    contacts_allfrags[frag]=contacts_thisfrag


#
out_filejson = open(parent_directory+'_polar_contacts_raw.json','w')
json.dump(contacts_allfrags, out_filejson, indent = 6)
out_filejson.close()

'''


jsjsjsjsjsj
jsjsjsjsjsjs
jsjsjsjsjs






sdfhsfgh



sdfsdf

sdfsdf

'''
import json
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import unique
from numpy import where
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

with open('dog_polar_contacts_raw.json','r') as f:
    df=json.load(f)

# df['Fragment_1']['dog_Fragment_1_LYS_contact_fuzzball.pdb']['resdonors']['11948']

#fragment:source pdb:resdonors/resacc:residue number in source file
testvals=[]
for resnum in df['Fragment_1']['dog_Fragment_1_LYS_contact_fuzzball.pdb']['resdonors'].keys():
    contactdata=df['Fragment_1']['dog_Fragment_1_LYS_contact_fuzzball.pdb']['resdonors'][resnum]
    contact_resname=contactdata[0]
    contact_atom=contactdata[8]
    if contact_atom!='N':
        if contact_atom!='O':
            if contact_atom!='CA':
                if contact_atom!='H':
                    dist_this_contact=float(contactdata[2])
                    hyb_this_contact=contactdata[-1]
                    psi=float(contactdata[3])
                    phi=float(contactdata[4])
                    et1=contactdata[6]
                    et2=contactdata[7]
                    #########
                    #########
                    #########
                    chi=float(abs(contactdata[5])) ######use absolute value and add 180 periodicity in csts?
                    # chi=float(contactdata[5]) ######use absolute value and add 180 periodicity in csts?
                    #########
                    #########
                    #########
                    if dist_this_contact<=2.5:      #limit distance to short range
                        # testvals.append((dist_this_contact,psi,phi,chi,et1,et2,resnum))
                        testvals.append((dist_this_contact,psi,phi,chi,resnum))
                    ##########
                    ##########

a=np.array(testvals)
a = a.astype('float64')
# a2=np.delete(a,6,1)
a2=np.delete(a,4,1)

# agglomerative clustering
scale = StandardScaler()
a2 = scale.fit_transform(a2)
# define the model
# model = AgglomerativeClustering(n_clusters=int(len(testvals)/4),linkage='single',affinity='cosine')
model = AgglomerativeClustering(n_clusters=int(len(testvals)/4))
# from sklearn.cluster import DBSCAN
# model =DBSCAN(eps=0.1)
# fit model and predict clusters
yhat = model.fit_predict(a2) #returns for each sample the label of the cluster its assigned to
# retrieve unique clusters
clusters = unique(yhat) #unique cliuster indices (labels)
#get points in each cluster
clust_resnums={}
clust_data={}
for cluster in clusters:
    clust_samples=[]
    clust_resnumss=[]
    row_ix = where(yhat == cluster)
    for ix in row_ix[0]:
        # clust_samples.append((a[ix, 0], a[ix, 1],a[ix, 2],a[ix, 3],a[ix, 4],a[ix, 5]))
        clust_samples.append((a[ix, 0], a[ix, 1],a[ix, 2],a[ix, 3]))
        # clust_resnumss.append(a[ix,6])
        clust_resnumss.append(a[ix,4])
    clust_data[cluster]=clust_samples
    clust_resnums[cluster]=clust_resnumss
#cluster populations
l=[]
for k in clust_data.keys():
    l.append(len(clust_data[k]))
l=sorted(l,reverse=True)
fig = plt.figure()
plt.bar([i+1 for i in range(len(l))],l)
plt.savefig('t3.pdf')
plt.clf()
#get overrep cluster pop value
fifth_perc_ind=int(0.03*len(l))
fifth_perc_val=l[fifth_perc_ind]
'''
scp cgalvin@log2.wynton.ucsf.edu:dog/t3.pdf ~/desktop/t3.pdf
'''
#plot of all points in data set colored by distance heatmap
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for c in clust_data.keys():
#     points=clust_data[c]
#     x = np.array([i[1] for i in points])
#     y = np.array([i[2] for i in points])
#     z = np.array([i[3] for i in points])
#     co = np.array([i[0] for i in points])
#     img = ax.scatter(x, y, z, c=co, cmap=plt.hot())
# fig.colorbar(img)
# plt.savefig('t.pdf')
# plt.clf()
'''
scp cgalvin@log2.wynton.ucsf.edu:dog/t.pdf ~/desktop/t.pdf
'''
#plot of all points colored by cluster index
# from matplotlib.pyplot import cm
# color = cm.rainbow(np.linspace(0, 1, len(l)))
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for cix,c in enumerate(list(clust_data.keys())):
#     points=clust_data[c]
#     x = np.array([i[1] for i in points])
#     y = np.array([i[2] for i in points])
#     z = np.array([i[3] for i in points])
#     img = ax.scatter(x, y, z, c=color[cix])
# ax.legend()
# plt.savefig('t2.pdf')
# plt.clf()
'''
scp cgalvin@log2.wynton.ucsf.edu:dog/t2.pdf ~/desktop/t2.pdf
'''
#find the value ranges that define each cluster by the max and min in that clust for each
#
clust_ranges={}
for c in clust_data.keys():
    points=clust_data[c]
    x = [i[1] for i in points]
    y = [i[2] for i in points]
    z = [i[3] for i in points]
    d = [i[0] for i in points]
    # j = [i[4] for i in points]
    # k = [i[5] for i in points]
    pop=len(points)
    dmin=min(d)
    dmax=max(d)
    psimin=min(x)
    psimax=max(x)
    phimin=min(y)
    phimax=max(y)
    chimin=min(z)
    chimax=max(z)
    # et1min=min(j)
    # et1max=max(j)
    # et2min=min(k)
    # et2max=max(k)
    if pop>=fifth_perc_val:
        # clust_ranges[c]=[pop,dmin,dmax,psimin,psimax,phimin,phimax,chimin,chimax,et1min,et1max,et2min,et2max]
        clust_ranges[c]=[pop,dmin,dmax,psimin,psimax,phimin,phimax,chimin,chimax]

#next thing is i need to actually look at these overrep clusters
#retrive file name and res number for the res in each cluster and spit them
#out to a new cluster fuzzball
overrep_clust_resnums={}
for k in clust_resnums.keys():
    pop=len(clust_resnums[k])
    if pop>=fifth_perc_val:
        overrep_clust_resnums[k]=clust_resnums[k]
# len(overrep_clust_resnums.keys())
# len(clust_ranges.keys())
import pandas as pd
import os

def pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
    f.close()
    df=pd.DataFrame(l)
    return df


#convert dataframe to pdb file
def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    saveresnum=1
    for i in range(dataframe.shape[0]):
        curr_resnum=int(dataframe.iloc[i]['resnum'].strip(' '))
        if curr_resnum!=saveresnum:
            of.write('TER\n')
            # print('ter')
            # print(str(curr_resnum))
            # print(str(saveresnum))
            line=''.join(dataframe.iloc[i])
            of.write(line)
            saveresnum=curr_resnum
        else:
            line=''.join(dataframe.iloc[i])
            of.write(line)
    of.close()

#convert ligand pdb to dataframe
def lig_pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:6]=='HETATM':
            d={}
            d['recordname']=line[:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=' X'
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
        else:
            pass
    f.close()
    df=pd.DataFrame(l)
    return df
#
og_fuzzball_path='/wynton/home/kortemme/cgalvin/dog/fragfbs/dog_Fragment_1_residue_contacts/dog_Fragment_1_LYS_contact_fuzzball.pdb'
df=pdb_to_df(og_fuzzball_path)
fullligdf=lig_pdb_to_df(og_fuzzball_path)
outpath='/wynton/home/kortemme/cgalvin/dog/fragfbs/dog_Fragment_1_residue_contacts/testclusts/'
if os.path.exists(outpath):
    os.system('rm -r '+outpath)
    os.makedirs(outpath,exist_ok=True)
os.makedirs(outpath,exist_ok=True)
#
rns=[]
rn=1
ll=list(df['resnum'])
for i,e in enumerate(ll):
    if i==0:
        rns.append(rn)
    else:
        if e==ll[i-1]:
            rns.append(rn)
        else:
            rn+=1
            rns.append(rn)

#
'''
okay resnum is not the listed residue number in file
but actually the index of the residue number, (pyrosetta resnum)
'''
for clustid in overrep_clust_resnums:
    cluster_df=pd.DataFrame()
    resnumbers=[int(i) for i in overrep_clust_resnums[clustid]]
    b=str(len(resnumbers))
    for resnum in resnumbers:
        try:
            rows=[i for i,e in enumerate(rns) if e == resnum]
            res_df=df.iloc[min(rows):max(rows)+1]
            cluster_df=cluster_df.append(res_df,ignore_index=True)
        except:
            print('cant find residue '+str(resnum))
    cluster_df=cluster_df.append(fullligdf,ignore_index=True)
    resnamee=str(res_df.iloc[0][4])
    ofilename=str(str(clustid)+'_'+resnamee+'_cluster_'+str(b)+'.pdb')
    df_to_pdb(cluster_df,os.path.join(outpath,ofilename))

'''
scp cgalvin@log2.wynton.ucsf.edu:dog/t2.pdf ~/desktop/t2.pdf
scp cgalvin@log2.wynton.ucsf.edu:dog/t.pdf ~/desktop/t.pdf
scp cgalvin@log2.wynton.ucsf.edu:dog/t3.pdf ~/desktop/t3.pdf

scp -r cgalvin@log2.wynton.ucsf.edu:dog/fragfbs/dog_Fragment_1_residue_contacts/testclusts ~/desktop/testclusts




CAN LIG HBONDING RESIDUES IN NAT BINDERS HBOND WITH EACH OTHER IN APO STATE?

max at 2.0 angs
The angle at the hydrogen atom shows
a clear preference for linearity for short sidechain–side-chain hydrogen bond
 Differences
between side-chain–side-chain and backbone–
backbone hydrogen bonds are also observed for
the angle at the acceptor, with maxima at around
1208 for side-chain–side-chain hydrogen bonds
and shifted to higher angles for the backbone–
backbone distributions.
The differences between
the different hybridization states at the acceptor
for side-chain–side-chain hydrogen bonds are
small, with a slightly sharper distribution (but not
shifted to significantly lower angles) for the sp3
case. A larger difference is seen comparing the geometries in short and long hydrogen bonds, with
broader distributions found in the latter case. The
dihedral angle shows fairly flat distributions
except for a well-defined peak in the case of helical
hydrogen bonds and other non-b backbone–backbone hydrogen bonds (mainly involving turns)


d 1.8-2.1
psi 110-130
phi 160-180
chi 180 +/20, periodicity 2


ASK DR SHOICHET FOR TRUE INACTIVES

do aromatic ligands prefer hbonding with aromatic residues?

MATCH NATIVE BINDING MOTIF INTO MY NTF2 MPNN LIB, TRY TO DEMONSTRATE
THE SCAFFOLD QUALITY FOR DOING THIS SHIT
using mpnn well scoring sequences as a kind of msa to inform sequence design choices
    does this improve af conf/rmsd of resultant designs
using special rot type approach during design
    is there evidence that


'''




#
#
# #im gonna make note of ambiguous H res and th
# ambiguous_res=['SER','THR','LYS','TYR']
# #########MAKE RANDOM DISTRIBUTIONS FOR NORMALIZATION##############
# ###########angle############
# lr=[]
# for i in range(10000):
#     x=random.uniform(0,3)
#     y=random.uniform(0,3)
#     z=random.uniform(0,3)
#     x2=random.uniform(0,3)
#     y2=random.uniform(0,3)
#     z2=random.uniform(0,3)
#     x3=random.uniform(0,3)
#     y3=random.uniform(0,3)
#     z3=random.uniform(0,3)
#     p1=[x,y,z]
#     p2=[x2,y2,z2]
#     p3=[x3,y3,z3]
#     a=getAngle(p1,p2,p3)
#     lr.append(a)
#
# # fig = plt.figure()
# # ax = fig.add_subplot()
# # ax.hist(lr,density=1,alpha=0.5,label='Raw Histogram')
# # plt.show()
# ###########torsion############
# lrt=[]
# for i in range(10000):
#     x=random.uniform(0,3)
#     y=random.uniform(0,3)
#     z=random.uniform(0,3)
#     x2=random.uniform(0,3)
#     y2=random.uniform(0,3)
#     z2=random.uniform(0,3)
#     x3=random.uniform(0,3)
#     y3=random.uniform(0,3)
#     z3=random.uniform(0,3)
#     x4=random.uniform(0,3)
#     y4=random.uniform(0,3)
#     z4=random.uniform(0,3)
#     p1=[x,y,z]
#     p2=[x2,y2,z2]
#     p3=[x3,y3,z3]
#     p4=[x4,y4,z4]
#     a=dihedral(p1,p2,p3,p4)
#     lrt.append(a)
#
# def qbin(l1,l2,lower_lim,upper_lim):
#     norm=1/((upper_lim**3)-(lower_lim**3))
#     qofbin=((l2**3)-(l1**3))*norm
#     return qofbin
# def qdist(d,lower_lim,upper_lim):
#     q_of_d=((d**2)/((1/3)*((upper_lim**3)-(lower_lim**3))))
#     return q_of_d
# #############NOW PLOTTING THINGS USING THESE REFERENCE DISTS TO NORMALIZE#################
# #############DISTANCE####################
# def distance_plots(alld,title_spec):
#     if len(alld)>0:
#         l=[]
#         for i in range(int(min(alld)*100),int(max(alld)*100)+5,5):
#             l.append(i/100)
#         rawcounts=[]
#         for i,e in enumerate(l[:-1]):
#             upperlim=l[i+1]
#             bincount=0
#             for d in alld:
#                 if e<=d<upperlim:
#                     bincount+=1
#             rawcounts.append((bincount,((e+upperlim)/2.)))
#         corrcounts=[]
#         for a,b in rawcounts:
#             c=a/(b**2)
#             corrcounts.append(c)
#         probs=[float((i/len(alld))) for i,b in rawcounts]
#         qx2=[]
#         qy=[]
#         for i,e in enumerate(l[:-1]):
#             lowerlim=e
#             upperlim=l[i+1]
#             qy.append(qbin(lowerlim,upperlim,1.4,3.5))
#             qx2.append(np.mean([lowerlim,upperlim]))
#         corrprobs=[]
#         for i,e in enumerate(probs):
#             z=(e/qy[i])*e
#             corrprobs.append(z)
#         # corrcounts=[i*(len(alld)) for i in corrprobs]
#         # fig = plt.figure()
#         # ax = fig.add_subplot()
#         # ax.scatter(qx2,qy,color='blue')
#         # ax.scatter(qx2,probs,color='red')
#         p=np.array(alld)
#         gkde=stats.gaussian_kde(p)
#         gkde.set_bandwidth(0.1)
#         ind = np.linspace(min(p),max(p),len(p))
#         kdepdf = gkde.evaluate(ind)
#         #
#         pofds = gkde.evaluate(p)
#         qofds=[qdist(i,min(p),max(p)) for i in p]
#         statens=[]
#         for i,e in enumerate(pofds):
#             qofd=qofds[i]
#             statens.append(-1*(math.log(e/qofd)))
#         #og hist normed automatically by matplotlib, kernel density estimate,
#         #reference distribution, stat potential
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.scatter(p,statens,color='red',label='Statistical Potential')
#         ax.scatter(p,qofds,color='yellow',label='Reference Distribution')
#         ax.hist(alld, bins=20,density=1,alpha=0.5,label='Raw Histogram')
#         ax.plot(ind, kdepdf, label='Kernel Density Estimate', color='g')
#         ax.set_title('Distribution of Donor-H Distances:'+str(title_spec))
#         ax.set_xlabel('Distance (Angstroms)')
#         ax.legend()
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#         #og histogram and weighted histogram
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.hist(alld, bins=l,alpha=0.5,color='blue',label='raw')
#         ax.bar(qx2,corrcounts,color='orange',alpha=0.5,width=0.05,label='corrected')
#         ax.legend()
#         ax.set_title('Distribution of Donor-H Distances:'+str(title_spec))
#         ax.set_xlabel('Distance (Angstroms)')
#         ax.set_ylabel('Count')
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#
#
# #####################angles############################
# def angle_plots(allphi,parametername,title_spec):
#     if len(allphi)>0:
#         phil=[]
#         for i in range(int(min(allphi)),int(max(allphi))+1,5):
#             phil.append(i)
#         rawcountsphi=[]
#         for i,e in enumerate(phil[:-1]):
#             upperlim=phil[i+1]
#             bincount=0
#             for phi in allphi:
#                 if e<=phi<upperlim:
#                     bincount+=1
#             rawcountsphi.append((bincount,((e+upperlim)/2.)))
#         probsphi=[float((i/len(allphi))) for i,b in rawcountsphi]
#         #
#         p=np.array(allphi)
#         gkde=stats.gaussian_kde(p)
#         gkde.set_bandwidth(0.1)
#         ind = np.linspace(0,180,len(p))
#         kdepdf = gkde.evaluate(ind)
#         #
#         pr=np.array(lr)
#         rkde=stats.gaussian_kde(pr)
#         rkdepdf = rkde.evaluate(ind)
#         #
#         pofangs = gkde.evaluate(p)
#         qofangs=rkde.evaluate(p)
#         statens=[]
#         for i,e in enumerate(pofangs):
#             qofang=qofangs[i]
#             statens.append(-1*(math.log(e/qofang)))
#         #og hist normed automatically by matplotlib, kernel density estimate,
#         #reference distribution
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.hist(lr,density=1,alpha=0.5,label='Random Histogram')
#         ax.hist(allphi,bins=phil,density=1,alpha=0.5,label='Raw Histogram')
#         ax.plot(ind, kdepdf, label='Kernel Density Estimate', color='g')
#         ax.plot(ind, rkdepdf, label='Kernel Density Estimate Random', color='yellow')
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_title(str(parametername)+' Distribution:'+str(title_spec))
#         ax.legend()
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#         #statistical potential
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.scatter(p,statens,color='red',label='Statistical Potential')
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_ylabel('E')
#         ax.set_title('Statistical Potential'+str(parametername)+str(title_spec))
#         ax.legend()
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#         qxa=[]
#         rawcountsphiq=[]
#         for i,e in enumerate(phil[:-1]):
#             upperlim=phil[i+1]
#             bincount=0
#             for phi in lr:
#                 if e<=phi<upperlim:
#                     bincount+=1
#             rawcountsphiq.append((bincount,((e+upperlim)/2.)))
#             qxa.append((e+upperlim)/2.)
#         probsphiq=[float((i/len(lr))) for i,b in rawcountsphiq]
#         corrprobsphi=[]
#         for i,e in enumerate(probsphi):
#             z=(e/probsphiq[i])*e
#             corrprobsphi.append(z)
#         # corrcountsphi=[i*(len(allphi)) for i in corrprobsphi]
#         corrcountsphi=[]
#         for a,b in rawcountsphi:
#             c=a/(math.sin(math.radians(b)))
#             corrcountsphi.append(c)
#         #og histogram and weighted histogram
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.hist(allphi,bins=phil,alpha=0.5,color='blue',label='raw')
#         ax.bar(qxa,corrcountsphi,color='orange',alpha=0.5,width=5,label='corrected')
#         ax.legend()
#         ax.set_title(str(parametername)+' Distribution:'+str(title_spec))
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_ylabel('Count')
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#
#
#
# #############CHI STUFF####################
# def torsion_plots(allchi,parametername,title_spec):
#     if len(allchi)>0:
#         chil=[]
#         for i in range(-180,181,10):
#             chil.append(i)
#         rawcountschi=[]
#         for i,e in enumerate(chil[:-1]):
#             upperlim=chil[i+1]
#             bincount=0
#             for chi in allchi:
#                 if e<=chi<upperlim:
#                     bincount+=1
#             rawcountschi.append((bincount,((e+upperlim)/2.)))
#         probschi=[float((i/len(allchi))) for i,b in rawcountschi]
#         #
#         p=np.array(allchi)
#         gkde=stats.gaussian_kde(p)
#         gkde.set_bandwidth(0.1)
#         ind = np.linspace(-180,180,len(p))
#         kdepdf = gkde.evaluate(ind)
#         #
#         pr=np.array(lrt)
#         rkde=stats.gaussian_kde(pr)
#         rkdepdf = rkde.evaluate(ind)
#         #
#         pofangs = gkde.evaluate(p)
#         qofangs=rkde.evaluate(p)
#         statens=[]
#         for i,e in enumerate(pofangs):
#             qofang=qofangs[i]
#             statens.append(-1*(math.log(e/qofang)))
#         #og hist normed automatically by matplotlib, kernel density estimate,
#         #reference distribution
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.hist(lrt,density=1,alpha=0.5,label='Random Histogram')
#         ax.hist(allchi,bins=50,density=1,alpha=0.5,label='Raw Histogram')
#         ax.plot(ind, kdepdf, label='Kernel Density Estimate', color='g')
#         ax.plot(ind, rkdepdf, label='Kernel Density Estimate Random', color='yellow')
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_title(str(parametername)+' Distribution:'+str(title_spec))
#         ax.legend()
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#         #statistical potential
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.scatter(p,statens,color='red',label='Statistical Potential')
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_ylabel('E')
#         ax.set_title('Statistical Potential '+str(parametername)+str(title_spec))
#         ax.legend()
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#         qxt=[]
#         rawcountschiq=[]
#         for i,e in enumerate(chil[:-1]):
#             upperlim=chil[i+1]
#             bincount=0
#             for chi in lrt:
#                 if e<=chi<upperlim:
#                     bincount+=1
#             rawcountschiq.append((bincount,((e+upperlim)/2.)))
#             qxt.append((e+upperlim)/2.)
#         probschiq=[float((i/len(lrt))) for i,b in rawcountschiq]
#         corrprobschi=[]
#         for i,e in enumerate(probschi):
#             z=(e/probschiq[i])*e
#             corrprobschi.append(z)
#         # corrcountschi=[i*(len(allchi)) for i in corrprobschi]
#         # corrcountschi=[]
#         # for a,b in rawcountschi:
#         #     c=a/(math.sin(math.radians(b)))
#         #     corrcountschi.append(c)
#         #og histogram and weighted histogram
#         fig = plt.figure()
#         ax = fig.add_subplot()
#         ax.hist(allchi,bins=chil,alpha=0.5,color='blue',label='raw')
#         # ax.bar(qxt,corrcountschi,color='orange',alpha=0.5,width=10,label='corrected')
#         ax.legend()
#         ax.set_title(str(parametername)+' Distribution: '+str(title_spec))
#         ax.set_xlabel(str(parametername)+' (Degrees)')
#         ax.set_ylabel('Count')
#         plt.tight_layout()
#         pdf.savefig()
#         plt.clf()
#         plt.close()
#
# ########now plot the stuff
# for frag in fraglist:
#     #lists for plotting distance data
#     sp2d_unambig=[]
#     sp3d_unambig=[]
#     sp2d_ambig=[]
#     sp3d_ambig=[]
#     #lists for plotting angle data
#     sp2_short_unambig_phi=[]
#     sp3_short_unambig_phi=[]
#     sp2_long_unambig_phi=[]
#     sp3_long_unambig_phi=[]
#     sp2_short_ambig_phi=[]
#     sp3_short_ambig_phi=[]
#     sp2_long_ambig_phi=[]
#     sp3_long_ambig_phi=[]
#     #
#     sp2_short_unambig_psi=[]
#     sp3_short_unambig_psi=[]
#     sp2_long_unambig_psi=[]
#     sp3_long_unambig_psi=[]
#     sp2_short_ambig_psi=[]
#     sp3_short_ambig_psi=[]
#     sp2_long_ambig_psi=[]
#     sp3_long_ambig_psi=[]
#     #
#     sp2_short_unambig_chi=[]
#     sp3_short_unambig_chi=[]
#     sp2_long_unambig_chi=[]
#     sp3_long_unambig_chi=[]
#     sp2_short_ambig_chi=[]
#     sp3_short_ambig_chi=[]
#     sp2_long_ambig_chi=[]
#     sp3_long_ambig_chi=[]
#     #THESE RESIDUES WILL BE CONSIDERED SEPARATELY FROM OTHERS BECAUSE THEIR SC
#     #HB GROUPS ARE SP3 AND H CANT BE UNAMBIGUOUSLY PLACED
#     contacts_thisfrag=contacts_allfrags[frag]
#     for clust in contacts_thisfrag.keys():
#         contacts_thisclust=contacts_thisfrag[clust]
#         for resdonor_contact in contacts_thisclust['resdonors'].keys():
#             contactdata=contacts_thisclust['resdonors'][resdonor_contact]
#             contact_resname=contactdata[0]
#             dist_this_contact=contactdata[2]
#             hyb_this_contact=contactdata[-1]
#             psi=contactdata[3]
#             phi=contactdata[4]
#             chi=contactdata[5]
#             if contact_resname in ambiguous_res:
#                 if hyb_this_contact=='SP2':
#                     sp2d_ambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp2_short_ambig_phi.append(phi)
#                         sp2_short_ambig_psi.append(psi)
#                         sp2_short_ambig_chi.append(chi)
#                     else:
#                         sp2_long_ambig_phi.append(phi)
#                         sp2_long_ambig_psi.append(psi)
#                         sp2_long_ambig_chi.append(chi)
#                 elif hyb_this_contact=='SP3':
#                     sp3d_ambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp3_short_ambig_phi.append(phi)
#                         sp3_short_ambig_psi.append(psi)
#                         sp3_short_ambig_chi.append(chi)
#                     else:
#                         sp3_long_ambig_phi.append(phi)
#                         sp3_long_ambig_psi.append(psi)
#                         sp3_long_ambig_chi.append(chi)
#                 else:
#                     print('unknown hybridization state for contact in '+str(clust))
#             else:
#                 if hyb_this_contact=='SP2':
#                     sp2d_unambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp2_short_unambig_phi.append(phi)
#                         sp2_short_unambig_psi.append(psi)
#                         sp2_short_unambig_chi.append(chi)
#                     else:
#                         sp2_long_unambig_phi.append(phi)
#                         sp2_long_unambig_psi.append(psi)
#                         sp2_long_unambig_chi.append(chi)
#                 elif hyb_this_contact=='SP3':
#                     sp3d_unambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp3_short_unambig_phi.append(phi)
#                         sp3_short_unambig_psi.append(psi)
#                         sp3_short_unambig_chi.append(chi)
#                     else:
#                         sp3_long_unambig_phi.append(phi)
#                         sp3_long_unambig_psi.append(psi)
#                         sp3_long_unambig_chi.append(chi)
#                 else:
#                     print('unknown hybridization state for contact in '+str(clust))
#         for resacc_contact in contacts_thisclust['resacc'].keys():
#             contactdata=contacts_thisclust['resacc'][resacc_contact]
#             contact_resname=contactdata[0]
#             dist_this_contact=contactdata[2]
#             hyb_this_contact=contactdata[-1]
#             psi=contactdata[3]
#             phi=contactdata[4]
#             chi=contactdata[5]
#             if contact_resname in ambiguous_res:
#                 if hyb_this_contact=='SP2':
#                     sp2d_ambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp2_short_ambig_phi.append(phi)
#                         sp2_short_ambig_psi.append(psi)
#                         sp2_short_ambig_chi.append(chi)
#                     else:
#                         sp2_long_ambig_phi.append(phi)
#                         sp2_long_ambig_psi.append(psi)
#                         sp2_long_ambig_chi.append(chi)
#                 elif hyb_this_contact=='SP3':
#                     sp3d_ambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp3_short_ambig_phi.append(phi)
#                         sp3_short_ambig_psi.append(psi)
#                         sp3_short_ambig_chi.append(chi)
#                     else:
#                         sp3_long_ambig_phi.append(phi)
#                         sp3_long_ambig_psi.append(psi)
#                         sp3_long_ambig_chi.append(chi)
#                 else:
#                     print('unknown hybridization state for contact in '+str(clust))
#             else:
#                 if hyb_this_contact=='SP2':
#                     sp2d_unambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp2_short_unambig_phi.append(phi)
#                         sp2_short_unambig_psi.append(psi)
#                         sp2_short_unambig_chi.append(chi)
#                     else:
#                         sp2_long_unambig_phi.append(phi)
#                         sp2_long_unambig_psi.append(psi)
#                         sp2_long_unambig_chi.append(chi)
#                 elif hyb_this_contact=='SP3':
#                     sp3d_unambig.append(dist_this_contact)
#                     if dist_this_contact<=2.1:
#                         sp3_short_unambig_phi.append(phi)
#                         sp3_short_unambig_psi.append(psi)
#                         sp3_short_unambig_chi.append(chi)
#                     else:
#                         sp3_long_unambig_phi.append(phi)
#                         sp3_long_unambig_psi.append(psi)
#                         sp3_long_unambig_chi.append(chi)
#                 else:
#                     print('unknown hybridization state for contact in '+str(clust))
#     try:
#         distance_plots(sp2d_unambig,'sp2d_unambig_'+str(frag))
#     except:
#         pass
#     try:
#         distance_plots(sp3d_unambig,'sp3d_unambig_'+str(frag))
#     except:
#         pass
#     try:
#         distance_plots(sp2d_ambig,'sp2d_ambig_'+str(frag))
#     except:
#         pass
#     try:
#         distance_plots(sp3d_ambig,'sp3d_ambig_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_short_unambig_phi,'phi','sp2_short_unambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_short_unambig_phi,'phi','sp3_short_unambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_long_unambig_phi,'phi','sp2_long_unambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_long_unambig_phi,'phi','sp3_long_unambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_short_ambig_phi,'phi','sp2_short_ambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_short_ambig_phi,'phi','sp3_short_ambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_long_ambig_phi,'phi','sp2_long_ambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_long_ambig_phi,'phi','sp3_long_ambig_phi_'+str(frag))
#     except:
#         pass
#     try:
#     #
#         angle_plots(sp2_short_unambig_psi,'psi','sp2_short_unambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_short_unambig_psi,'psi','sp3_short_unambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_long_unambig_psi,'psi','sp2_long_unambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_long_unambig_psi,'psi','sp3_long_unambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_short_ambig_psi,'psi','sp2_short_ambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_short_ambig_psi,'psi','sp3_short_ambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp2_long_ambig_psi,'psi','sp2_long_ambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         angle_plots(sp3_long_ambig_psi,'psi','sp3_long_ambig_psi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp2_short_unambig_chi,'chi','sp2_short_unambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp3_short_unambig_chi,'chi','sp3_short_unambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp2_long_unambig_chi,'chi','sp2_long_unambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp3_long_unambig_chi,'chi','sp3_long_unambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp2_short_ambig_chi,'chi','sp2_short_ambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp3_short_ambig_chi,'chi','sp3_short_ambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp2_long_ambig_chi,'chi','sp2_long_ambig_chi_'+str(frag))
#     except:
#         pass
#     try:
#         torsion_plots(sp3_long_ambig_chi,'chi','sp3_long_ambig_chi_'+str(frag))
#     except:
#         pass
#



'''
CLUSTERING
'''
# dd=[i for i in range(1,21)]
# aa=[i for i in range(1,18)]
# aa2=[i for i in range(1,18)]
# tt=[i for i in range(1,36)]
# lol=[dd,aa,aa2,tt]
# loll=list(itertools.product(*lol))
#
# all_clusters={}
# #
# for frag in fraglist:
#     fragclusters=defaultdict(list)
#     contacts_thisfrag=contacts_allfrags[frag]
#     for clust in contacts_thisfrag.keys():
#         contacts_thisclust=contacts_thisfrag[clust]
#         for resdonor_contact in contacts_thisclust['resdonors'].keys():
#             contactdata=contacts_thisclust['resdonors'][resdonor_contact]
#             contact_resname=contactdata[0]
#             dist_this_contact=contactdata[2]
#             hyb_this_contact=contactdata[-1]
#             psi=contactdata[3]
#             phi=contactdata[4]
#             chi=contactdata[5]
#             sres1=contactdata[8]
#             sres2=contactdata[9]
#             sres3=contactdata[10]
#             slig1=contactdata[11]
#             slig2=contactdata[12]
#             slig3=contactdata[13]
#             dbin=math.ceil((dist_this_contact-1.4)/0.2)
#             # print(dist_this_contact)
#             # print(dbin)
#             # print(phi)
#             phibin=math.ceil(phi/20)
#             # print(phibin)
#             psibin=math.ceil(psi/20)
#             if hyb_this_contact=='SP2':
#                 chibin=math.ceil((chi+180)/20)
#                 cluster_id='_'.join([sres1,sres2,sres3,slig1,slig2,slig3,contact_resname,str(dbin),str(phibin),str(psibin),str(chibin)])
#                 fragclusters[cluster_id].append((clust,resdonor_contact))
#             else:
#                 cluster_id='_'.join([sres1,sres2,sres3,slig1,slig2,slig3,contact_resname,str(dbin),str(phibin),str(psibin)])
#                 fragclusters[cluster_id].append((clust,resdonor_contact))
#         for resacc_contact in contacts_thisclust['resacc'].keys():
#             contactdata=contacts_thisclust['resacc'][resacc_contact]
#             contact_resname=contactdata[0]
#             dist_this_contact=contactdata[2]
#             hyb_this_contact=contactdata[-1]
#             psi=contactdata[3]
#             phi=contactdata[4]
#             chi=contactdata[5]
#             sres1=contactdata[8]
#             sres2=contactdata[9]
#             sres3=contactdata[10]
#             slig1=contactdata[11]
#             slig2=contactdata[12]
#             slig3=contactdata[13]
#             dbin=math.ceil((dist_this_contact-1.4)/0.2)
#             # print(dist_this_contact)
#             # print(dbin)
#             # print(phi)
#             phibin=math.ceil(phi/20)
#             # print(phibin)
#             psibin=math.ceil(psi/20)
#             if hyb_this_contact=='SP2':
#                 chibin=math.ceil((chi+180)/20)
#                 cluster_id='_'.join([sres1,sres2,sres3,slig1,slig2,slig3,contact_resname,str(dbin),str(phibin),str(psibin),str(chibin)])
#                 fragclusters[cluster_id].append((clust,resacc_contact))
#             else:
#                 cluster_id='_'.join([sres1,sres2,sres3,slig1,slig2,slig3,contact_resname,str(dbin),str(phibin),str(psibin)])
#                 fragclusters[cluster_id].append((clust,resacc_contact))
#     all_clusters[frag]=fragclusters
#
#
# #plotting clust pop distributions for each frag
# clust_output_dict={}
# for frag in list(all_clusters.keys()):
#     pops=[]
#     for clustbin in all_clusters[frag].keys():
#         pops.append((len(all_clusters[frag][clustbin]),clustbin))
#     pops=sorted(pops,reverse=True,key=lambda x: x[0])
#     # xs=[i for i in range(1,len(pops)+1)]
#     # ys=[i[0] for i in pops]
#     # fig = plt.figure()
#     # ax = fig.add_subplot()
#     # ax.bar(xs,ys)
#     # ax.set_title('Cluster Populations: '+str(frag))
#     # ax.set_xlabel('Cluster Index')
#     # ax.set_ylabel('Population')
#     # # ax.axvline(int(0.05*len(xs)),c='g',linestyle='dashed')
#     # ax.axvline(30,c='g',linestyle='dashed')
#     # pdf.savefig()
#     # plt.clf()
#     # plt.close()
#     clust_output_dict[frag]=pops
#
#
#
# # pdf.close()
#
# #output clusters to json, this will be input for cst generation
# out_filejson2 = open(parent_directory+'_polar_contacts_clustered_liberal.json','w')
# json.dump(clust_output_dict, out_filejson2, indent = 6)
# out_filejson2.close()
