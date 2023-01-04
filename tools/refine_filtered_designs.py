'''
/Users/student/Desktop/esl4rmfilt_testfirstgo
In [2]: len(l)
Out[2]: 1233

print(tf.create_task_and_apply_taskoperations(test_pose))



12.24 STATE OF THE CODE
ive got it such that ligand polar atoms that dont have an hbond are identified
    all fs res that dont already hb with lig then mutated to styk, minimized
        if now hbonds with lig, mutation saved
            the most favorable energetically of these mutations is applied to the pose
then, all polar res in fs that dont hb with lig identified
    check that they dont hb with anything that DOES hb with lig (preserve networks)
        if not, mutate to nonpolar res ilvf, minimize
            most favorable of these is saved and applied to pose
these mutated poses are then output to a new directory
all structures without any of these issues copied to this new directory
I AM THINKING, SHOULD I CHANGE THE MUTATION PREFERENCES TO ONLY FAVORED
RES FOR THE SECONDARY STRUCTURE THAT THE RESIDUE IS A PART OF?
-----------------------------------------------------------------------------
sstype=pyrosetta.rosetta.core.scoring.dssp.Dssp(p).get_dssp_secstruct(resnum)
G = 3-turn helix (310 helix). Min length 3 residues.
H = 4-turn helix (alpha helix). Min length 4 residues.
I = 5-turn helix (pi helix). Min length 5 residues.
T = hydrogen bonded turn (3, 4 or 5 turn)
E = beta sheet in parallel and/or anti-parallel sheet conformation (extended strand). Min length 2 residues.
B = residue in isolated beta-bridge (single pair beta-sheet hydrogen bond formation)
S = bend (the only non-hydrogen-bond based assignment)

Ala, Glu, Leu, and Met are most often found in helices whereas, Gly, Tyr, Ser, and Pro are less likely to be seen

Large aromatic residues (tyrosine, phenylalanine, tryptophan) and β-branched amino
acids (threonine, valine, isoleucine)
are favored to be found in β-strands in the middle of β-sheets.
Different types of residues (such as proline) are likely to be found in the edge strands in β-sheets, presumably to avoid the "edge-to-edge" association between proteins that might lead to aggregation and amyloid formation
-----------------------------------------------------------------------------
IN ADDITION, THERE ARE THE MATTERS OF ADDRESSING
    exposed hydrophobics
    making small res bigger (perhaps also taking into account ss pref here)
-----------------------------------------------------------------------------
FINALLY, THIS IS ONLY CAPABLE OF FIXING 1 MISSING LIG HB AT THE MOMENT
    because it doesnt actually check for hbonds between ptn and lig missing atom,.
    it just finds a mutation that creates a new ptn-lig hbond
1.3.2023
    CHANGED IT SO THAT IT ITERATES THRU THE HB MUT PROCESS
    FOR AS MANY TIMES AS THERE ARE UNSAT LIG ATOMS
    IN THIS WAY I THINK IT IS NOW TECHNICALLY CAPABLE OF FINDING MUTATIONS TO
    SATISFY ALL SUCH CASES
-----------------------------------------------------------------------------



'''


import sys
prm1=sys.argv[1]
pdb=sys.argv[2]
curr_soln_dirname=sys.argv[3]
#
import os
# l=[i for i in os.listdir() if i[-3:]=='pdb']
#
from pyrosetta import *
from pyrosetta.toolbox import mutate_residue
init('-load_PDB_components False -ex1 -ex2')

#
sf=get_fa_scorefxn()
#
polar_res=['S','T','Y','K'] #res to mutate to looking for new hbonds
polar_residues=['SER','THR','ARG','LYS','ASN','GLN','ASP','GLU'] #res to find in bs and try mutate to np if not hbonding
npres=['L','I','V','F'] #npres to mutate to
#
#
# test_count=1
#
noligunsats=[]
nouseless_buried_polar=[]
if not os.path.exists(curr_soln_dirname):
    os.makedirs(curr_soln_dirname,exist_ok=True)
#load pose
lig=[prm1]
p=Pose()
generate_nonstandard_residue_set(p,lig)
pose_from_file(p, pdb)
p.update_residue_neighbors()
#identify first shell residues
ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
neighborhood_selector_bool = neighborhood_selector.apply(p)
neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
first_shell_res=list(neighborhood_residues_resnums)
#get the sasa of first shell res to make sure they're facing binding site
print(first_shell_res)
for residd in first_shell_res:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=12.:
        first_shell_res.remove(residd)
print(first_shell_res)
#this is one of those stupid instances wbere i have to filter a list twice
#for no apparent reason, pisses me off -.-
for residd in first_shell_res:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=12.:
        first_shell_res.remove(residd)
print(first_shell_res)
#find whether all ligand polar atoms hbond
hbond_set = rosetta.core.scoring.hbonds.HBondSet()
p.update_residue_neighbors()
rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
#get atom indices of ligand polar atoms in network pose
lig_pol_atoms=[]
for i in range(1,len(p.residue(p.total_residue()).atoms())+1):
    s=str(p.residue(p.total_residue()).atom_type(i))
    ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()
    if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
        lig_pol_atoms.append(i)
        print(ligatom)
#identify any ligand polar atoms with no hbonds
hbond_data=[]
residues_hb_with_lig=[]
s=p.sequence()
if hbond_set.nhbonds()>0:
    for hbond_index in range(1,hbond_set.nhbonds()+1):
        drip=hbond_set.hbond(hbond_index).don_res_is_protein()
        arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
        if drip==True and arip==True:
            continue
        donres_ind=int(hbond_set.hbond(hbond_index).don_res())
        accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
        acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
        donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
        don_atom_index=int(p.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
        if drip==True and accres_ind==p.total_residue():
            hbond_data.append(acc_atom_index)
            residues_hb_with_lig.append(donres_ind)
        elif drip==False and donres_ind==p.total_residue():
            hbond_data.append(don_atom_index)
            residues_hb_with_lig.append(accres_ind)
unsat_lig_atoms=[]
for atind in lig_pol_atoms:
    if atind not in hbond_data:
        unsat_lig_atoms.append(atind)
############################################################
############################################################
############################################################
############################################################
############################################################
ulp=0
ghb=0
gm=0
#mutate first shell residues not already hbonding with lig to polar, minimize,
#see if theres hbond now
solncount=0
if len(unsat_lig_atoms)==0:
    noligunsats.append(pdb)
else:
    for xx in range(len(unsat_lig_atoms)):
        good_hbond_mutations=[]
        for resid in first_shell_res:
            if resid not in residues_hb_with_lig:
                for polrestype in polar_res:
                    test_pose=p.clone()
                    start_energy=sf(test_pose)
                    mutate_residue(test_pose,resid,polrestype)
                    mm = MoveMap()
                    mm.set_bb(False)
                    mm.set_chi(False)
                    mm.set_bb(resid,True)
                    mm.set_chi(resid,True) ## For side chain minimization
                    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
                    minmover.movemap(mm)
                    minmover.score_function(sf)
                    minmover.apply(test_pose)
                    w=test_pose.energies().energy_graph().find_energy_edge(resid,test_pose.total_residue())
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5:
                        print('\n\n\n')
                        print(resid)
                        print(polrestype)
                        print(hbsc)
                        print('\n\n\n')
                        # not_fs_selector_res=list(not_fs_selector_resnums)
                        new_energy=sf(test_pose)
                        print('\n\n\n')
                        print(start_energy)
                        print(new_energy)
                        print('\n\n\n')
                        if new_energy<=start_energy+2.:
                            # temp_pdb_name=pdb.split('.')[0]+'_hbonds_fixed_'+str(solncount)+'.pdb'
                            solncount+=1
                            # test_pose.dump_pdb('temp_solutions/'+temp_pdb_name)
                            good_hbond_mutations.append((new_energy,p.residue(resid).name(),resid,polrestype))
        if len(good_hbond_mutations)>0:
            ghb+=1
            good_hbond_mutations=sorted(good_hbond_mutations, key=lambda first: first[0])
            best_hbond_mut=good_hbond_mutations[0]
            best_hbond_mut_resid=best_hbond_mut[2]
            best_hbond_mut_restype=best_hbond_mut[3]
            test_pose=p.clone()
            mutate_residue(test_pose,best_hbond_mut_resid,best_hbond_mut_restype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            p=test_pose.clone()
            residues_hb_with_lig.append(best_hbond_mut_resid)
        else:
            if len(residues_hb_with_lig)>=len(lig_pol_atoms):
                ulp+=1
                pass
            else:
                print('\n\n\nNO GOOD HB MUTS')
                sys.exit()
############################################################
############################################################
############################################################
############################################################
############################################################
#identify polar residues in binding site that are not hydrogen
#bonding with the ligand
useless_buried_polar=[]
for resnum in first_shell_res:
    restype=p.residue(resnum).name()[:3]
    if restype in polar_residues:
        if resnum not in residues_hb_with_lig:
            important_hb=0
            for resnum2 in residues_hb_with_lig:
                try:
                    w=p.energies().energy_graph().find_energy_edge(resnum,resnum2)
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5:
                        important_hb+=1
                except:
                    pass
            if important_hb>0:
                pass
            else:
                useless_buried_polar.append(resnum)
#mutate them to aliphatics
solncount2=0
if len(useless_buried_polar)==0:
    nouseless_buried_polar.append(pdb)
else:
    for resid in useless_buried_polar:
        good_mutations=[]
        for nprestype in npres:
            test_pose=p.clone()
            start_energy=sf(test_pose)
            mutate_residue(test_pose,resid,nprestype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            new_energy=sf(test_pose)
            print('\n\n\n')
            print(start_energy)
            print(new_energy)
            print('\n\n\n')
            if new_energy<=start_energy+2.:
                # temp_pdb_name=pdb.split('.')[0]+'_hbonds_fixed_'+str(solncount)+'_'+str(solncount2)+'.pdb'
                solncount2+=1
                # test_pose.dump_pdb('temp_solutions/'+temp_pdb_name)
                good_mutations.append((new_energy,p.residue(resid).name(),resid,nprestype))
        if len(good_mutations)>0:
            gm+=1
            good_mutations=sorted(good_mutations, key=lambda first: first[0])
            best_mut=good_mutations[0]
            best_mut_resid=best_mut[2]
            best_mut_restype=best_mut[3]
            test_pose=p.clone()
            mutate_residue(test_pose,best_mut_resid,best_mut_restype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            p=test_pose.clone()
        else:
            pass
if ulp>0:
    print('\n\n\nUNSAT LIG POL BUT HAS MULTI HB WITH ANOTHER\n\n\n')
if ghb>0:
    print('\n\n\nHB MUT\n\n\n')
if gm>0:
    print('\n\n\nEP MUT\n\n\n')
if len(unsat_lig_atoms)==0 and len(useless_buried_polar)==0:
    os.system('cp '+pdb+' '+curr_soln_dirname+'/'+pdb.split('.')[0]+'_unrefined.pdb')
    print('\n\n\nNO PROBLEMS\n\n\n')
elif len(unsat_lig_atoms)>0 and len(useless_buried_polar)==0:
    # os.system('cp '+pdb+' '+curr_soln_dirname+'/'+pdb.split('.')[0]+'_og.pdb')
    soln_pdb_name=pdb.split('.')[0]+'_refined.pdb'
    p.dump_pdb(curr_soln_dirname+'/'+soln_pdb_name)
    print('UNSAT ONLY')
elif len(unsat_lig_atoms)==0 and len(useless_buried_polar)>0:
    # os.system('cp '+pdb+' '+curr_soln_dirname+'/'+pdb.split('.')[0]+'_og.pdb')
    soln_pdb_name=pdb.split('.')[0]+'_refined.pdb'
    p.dump_pdb(curr_soln_dirname+'/'+soln_pdb_name)
    print('EXTRANEOUS POLAR ONLY')
elif len(unsat_lig_atoms)>0 and len(useless_buried_polar)>0:
    # os.system('cp '+pdb+' '+curr_soln_dirname+'/'+pdb.split('.')[0]+'_og.pdb')
    soln_pdb_name=pdb.split('.')[0]+'_refined.pdb'
    p.dump_pdb(curr_soln_dirname+'/'+soln_pdb_name)
    print('BOTH')



'''

unsat_lig_atoms
useless_buried_polar


nouseless_buried_polar=[]


useless_buried_polar=[]
for resnum in first_shell_res:
    restype=p.residue(resnum).name()[:3]
    if restype in polar_residues:
        if resnum not in residues_hb_with_lig:
            useless_buried_polar.append(resnum)
solncount2=1
if len(useless_buried_polar)==0:
    nouseless_buried_polar.append(pdb)
else:
    good_mutations=[]
    for resid in useless_buried_polar:
        for nprestype in npres:
            test_pose=p.clone()
            start_energy=sf(test_pose)
            mutate_residue(test_pose,resid,nprestype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            new_energy=sf(test_pose)
            print('\n\n\n')
            print(start_energy)
            print(new_energy)
            print('\n\n\n')
            if new_energy<=start_energy+2.:
                temp_pdb_name=pdb.split('.')[0]+'_hbonds_fixed_'+str(solncount)+'_'+str(solncount2)+'.pdb'
                solncount2+=1
                test_pose.dump_pdb('temp_solutions/'+temp_pdb_name)
                good_mutations.append((new_energy,p.residue(resid).name(),resid,nprestype,temp_pdb_name))
good_mutations=sorted(good_mutations, key=lambda first: first[0])
best_mut=good_mutations[0]
best_mut_resid=best_mut[2]
best_mut_restype=best_mut[3]
test_pose=p.clone()
mutate_residue(test_pose,best_mut_resid,best_mut_restype)
mm = MoveMap()
mm.set_bb(False)
mm.set_chi(False)
mm.set_bb(resid,True)
mm.set_chi(resid,True) ## For side chain minimization
minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
minmover.movemap(mm)
minmover.score_function(sf)
minmover.apply(test_pose)
p=test_pose.clone()





individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
individual_res_selector.set_index(residd)
individual_res_selector.apply(p)

for residd in first_shell_res:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=18:
        first_shell_res.remove(residd)
    print('\n')
    print(residd)
    print(p.residue(residd).name())
    print(this_res_sasa)
    print('\n')

test_pose.dump_pdb('temp_solutions/test.pdb')
scorefxn = create_score_function_ws_patch("standard","score12")
mm = MoveMap()
mm.set_bb(True)
pose_move_map.set_chi(True) ## For side chain minimization
mm.set_bb_true_range(1,4) ## Here range means 1 to 4 or the minimzation is applied only to residues 1 and4 ??

minmover = MinMover(mm, scorefxn, 'dfpmin', 10, True) ## I don't know the meaning of dfpmin,10,True. I saw it somewhere and used it

minmover.apply(pose)



#repacking tf
tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking()) # should be operation.PreventRepacking()?
prevent = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT() # should be operation.RestrictToRepackRLT() ?
restrict_to_focus = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent, neighborhood_selector, True)
tf.push_back(restrict_to_focus)
#repack and minimize first shell
packer=pyrosetta.rosetta.protocols.minimization_packing.MinPackMover(sf)
# packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sf)
packer.task_factory(tf)
print(packer)
packer.apply(test_pose)
new_energy=sf(test_pose)
print('\n\n\n')
print(start_energy)
print(new_energy)
print('\n\n\n')
if new_energy<=start_energy+10.0:
    #now check again and make sure no unsat lig pol atoms
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    test_pose.update_residue_neighbors()
    rosetta.core.scoring.hbonds.fill_hbond_set(test_pose, False, hbond_set)
    #get atom indices of ligand polar atoms in network pose
    lig_pol_atoms=[]
    for i in range(1,len(test_pose.residue(test_pose.total_residue()).atoms())+1):
        s=str(test_pose.residue(test_pose.total_residue()).atom_type(i))
        ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()
        if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
            lig_pol_atoms.append(i)
            print(ligatom)
    #identify any ligand polar atoms with no hbonds
    hbond_data=[]
    residues_hb_with_lig=[]
    s=test_pose.sequence()
    if hbond_set.nhbonds()>0:
        for hbond_index in range(1,hbond_set.nhbonds()+1):
            drip=hbond_set.hbond(hbond_index).don_res_is_protein()
            arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
            if drip==True and arip==True:
                continue
            donres_ind=int(hbond_set.hbond(hbond_index).don_res())
            accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
            acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
            donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
            don_atom_index=int(test_pose.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
            if drip==True and accres_ind==test_pose.total_residue():
                hbond_data.append(acc_atom_index)
                residues_hb_with_lig.append(donres_ind)
            elif drip==False and donres_ind==test_pose.total_residue():
                hbond_data.append(don_atom_index)
                residues_hb_with_lig.append(accres_ind)
    unsat_lig_atoms=[]
    for atind in lig_pol_atoms:
        if atind not in hbond_data:
            unsat_lig_atoms.append(atind)
    if len(unsat_lig_atoms)==0:
        temp_pdb_name=pdb.split('.')[0]+'hbonds_fixed_'+str(solncount)+'.pdb'
        solncount+=1
        test_pose.dump_pdb('temp_solutions/'+temp_pdb_name)
        good_hbond_mutations.append((resid,polrestype,temp_pdb_name))



if new_energy<=start_energy+10.0:
    #now check again and make sure no unsat lig pol atoms
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    test_pose.update_residue_neighbors()
    rosetta.core.scoring.hbonds.fill_hbond_set(test_pose, False, hbond_set)
    #get atom indices of ligand polar atoms in network pose
    lig_pol_atoms=[]
    for i in range(1,len(test_pose.residue(test_pose.total_residue()).atoms())+1):
        s=str(test_pose.residue(test_pose.total_residue()).atom_type(i))
        ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()
        if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
            lig_pol_atoms.append(i)
            print(ligatom)
    #identify any ligand polar atoms with no hbonds
    hbond_data=[]
    residues_hb_with_lig=[]
    s=test_pose.sequence()
    if hbond_set.nhbonds()>0:
        for hbond_index in range(1,hbond_set.nhbonds()+1):
            drip=hbond_set.hbond(hbond_index).don_res_is_protein()
            arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
            if drip==True and arip==True:
                continue
            donres_ind=int(hbond_set.hbond(hbond_index).don_res())
            accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
            acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
            donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
            don_atom_index=int(test_pose.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
            if drip==True and accres_ind==test_pose.total_residue():
                hbond_data.append(acc_atom_index)
                residues_hb_with_lig.append(donres_ind)
            elif drip==False and donres_ind==test_pose.total_residue():
                hbond_data.append(don_atom_index)
                residues_hb_with_lig.append(accres_ind)
    unsat_lig_atoms=[]
    for atind in lig_pol_atoms:
        if atind not in hbond_data:
            unsat_lig_atoms.append(atind)
    if len(unsat_lig_atoms)==0:
        temp_pdb_name=pdb.split('.')[0]+'hbonds_fixed_'+str(solncount)
        solncount+=1
        test_pose.dump_pdb('temp_solutions/'+temp_pdb_name)
        good_hbond_mutations.append((resid,polrestype,temp_pdb_name))


print(tf.create_task_and_apply_taskoperations(test_pose))

'''

'''
# Create TaskFactory
tf = TaskFactory()
tf.push_back(operation.InitializeFromCommandline())
tf.push_back(operation.RestrictToRepacking()) # should be operation.PreventRepacking()?
tf.push_back(restrict_to_focus)


   # Select specific residues
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus(str(mutant_position))
    nbr_selector.set_distance(pack_radius)
    nbr_selector.set_include_focus_in_subset(True)

    prevent = operation.PreventRepackingRLT() # should be operation.RestrictToRepackRLT() ?
    restrict_to_focus = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent, nb
r_selector, True)

    # Create TaskFactory
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking()) # should be operation.PreventRepacking()?
    tf.push_back(restrict_to_focus)

    # Repack
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(pack_scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)


#identify first shell residues open to mutation (does not hbond with ligand
#or have fantastic score)
#
#
# #find whether all ligand polar atoms hbond
# network_pose=Pose()
# ligand_pose=p.residue(p.total_residue()).clone()
# network_pose.append_residue_by_jump(ligand_pose, 1)
# res_pose1=p.residue(first_shell_res[0]).clone()
# network_pose.append_residue_by_jump(res_pose1, 1)
# c=1
# for i,resnum in enumerate(first_shell_res[1:]):
#     res_pose=p.residue(resnum).clone()
#     last_num=first_shell_res[i-1]
#     if last_num==resnum-1:
#         network_pose.append_residue_by_jump(res_pose, c)
#     else:
#         c+=1
#         network_pose.append_residue_by_jump(res_pose, c)
# hbond_set = rosetta.core.scoring.hbonds.HBondSet()
# network_pose.update_residue_neighbors()
# rosetta.core.scoring.hbonds.fill_hbond_set(network_pose, False, hbond_set)
# #get atom indices of ligand polar atoms in network pose
# lig_pol_atoms=[]
# for i in range(1,len(network_pose.residue(1).atoms())+1):
#     s=str(network_pose.residue(1).atom_type(i))
#     ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()
#     if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
#         lig_pol_atoms.append(i)
#         print(ligatom)
# #identify any ligand polar atoms with no hbonds
# hbond_data=[]
# s=network_pose.sequence()
# if hbond_set.nhbonds()>0:
#     for hbond_index in range(1,hbond_set.nhbonds()+1):
#         drip=hbond_set.hbond(hbond_index).don_res_is_protein()
#         arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
#         if drip==True and arip==True:
#             continue
#         donres_ind=int(hbond_set.hbond(hbond_index).don_res())
#         accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
#         acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
#         donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
#         don_atom_index=int(network_pose.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
#         if drip==True and arip==False:
#             hbond_data.append(acc_atom_index)
#         elif drip==False and arip==True:
#             hbond_data.append(don_atom_index)
# unsat_lig_atoms=[]
# for atind in lig_pol_atoms:
#     if atind not in hbond_data:
#         unsat_lig_atoms.append(atind)
#





#find whether all ligand polar atoms hbond
hbond_set = rosetta.core.scoring.hbonds.HBondSet()
p.update_residue_neighbors()
rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
#get atom indices of ligand polar atoms in network pose
lig_pol_atoms=[]
for i in range(1,len(p.residue(p.total_residue()).atoms())+1):
    s=str(p.residue(p.total_residue()).atom_type(i))
    ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()
    if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
        lig_pol_atoms.append(i)
        print(ligatom)
#identify any ligand polar atoms with no hbonds
hbond_data=[]
residues_hb_with_lig=[]
s=p.sequence()
if hbond_set.nhbonds()>0:
    for hbond_index in range(1,hbond_set.nhbonds()+1):
        drip=hbond_set.hbond(hbond_index).don_res_is_protein()
        arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
        if drip==True and arip==True:
            continue
        donres_ind=int(hbond_set.hbond(hbond_index).don_res())
        accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
        acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
        donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
        don_atom_index=int(p.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
        if drip==True and accres_ind==p.total_residue():
            hbond_data.append(acc_atom_index)
            residues_hb_with_lig.append(donres_ind)
        elif drip==False and donres_ind==p.total_residue():
            hbond_data.append(don_atom_index)
            residues_hb_with_lig.append(accres_ind)
unsat_lig_atoms=[]
for atind in lig_pol_atoms:
    if atind not in hbond_data:
        unsat_lig_atoms.append(atind)
#mutate first shell residues not already hbonding with lig to polar, minimize,
#see if theres hbond now
for resid in first_shell_res:
    if resid not in residues_hb_with_lig:
        for polrestype in polar_res:
            test_pose=p.clone()
            mutate_residue(test_pose,resid,polrestype)
            w=p.energies().energy_graph().find_energy_edge(resid,test_pose.total_residue())
            w.fill_energy_map()
            hbsc=w[rosetta.core.scoring.hbond_sc]
            print(resid)
            print(polrestype)
            print(hbsc)

from pyrosetta.toolbox import mutate_residue


polar_res=['S','T','Y','K']



def mutate_residue( pose , mutant_position , mutant_aa ,
        pack_radius = 4.0 , pack_scorefxn = '' ):
    """
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <pack_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose,30,A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence
    """
    #### a MutateResidue Mover exists similar to this except it does not pack
    ####    the area around the mutant residue (no pack_radius feature)
    #mutator = MutateResidue( mutant_position , mutant_aa )
    #mutator.apply( test_pose )

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = create_score_function( 'standard' )

    task = standard_packer_task( test_pose )

    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()
    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range( 1 , 21 ):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task( mutant_position
        ).restrict_absent_canonical_aas( aa_bool )

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1 , pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue( i ).nbr_atom_xyz() ) > pack_radius**2:
            task.nonconst_residue_task( i ).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover( pack_scorefxn , task )
    packer.apply( test_pose )

    return test_pose





network_pose.residue(1).atom(1).

#
lig_atom_type_hbonds=[]
hb_w_lig=[]
ptn_lig_hb=[]
for keyb in hbond_data.keys():
    l=hbond_data[keyb]
    if l[0]==1:
        hb_w_lig.append(l[1])
        network_members.append((l[0],l[1],l[-2],l[-1]))
        ptn_lig_hb.append((l[-2],l[-3],'ligdon'))
        lig_atom_type_hbonds.append((l[6],l[7],l[8]))
    if l[1]==1:
        hb_w_lig.append(l[0])
        network_members.append((l[0],l[1],l[-2],l[-1]))
        ptn_lig_hb.append((l[-1],l[-4],'ligacc'))
        lig_atom_type_hbonds.append((l[7],l[6],l[8]))
for i in hb_w_lig:
    hbws=[]
    for keyc in hbond_data.keys():
        l=hbond_data[keyc]
        if l[0]==i:
            if l[1]!=1:
                hbws.append(l[1])
                network_members.append((l[0],l[1],l[-2],l[-1]))
        if l[1]==i:
            if l[0]!=1:
                hbws.append(l[0])
                network_members.append((l[0],l[1],l[-2],l[-1]))
data_across_strc[sta]=[pdbid,ligname,naroatoms,netq,narofirstshell,
                       nposfirstshell,nnegfirstshell,bscharge,npotentialhb,
                       alltypes,
                       data_each_res,lig_atom_type_hbonds,network_members]



#now analysis of first shell residues
ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
neighborhood_selector_bool = neighborhood_selector.apply(p)
neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
first_shell_res=list(neighborhood_residues_resnums)
#
sf(p)
rosetta.core.pack.optimizeH(p, sf)
all_sec_res=[]
data_each_res=[]
narofirstshell=0
nposfirstshell=0
nnegfirstshell=0
for resnum in first_shell_res:
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
    if float(res_lig_score)<0.0:
        #residue name
        restype=p.residue(resnum).name()[:3]
        if restype in arores:
            narofirstshell+=1
        elif restype in posq:
            nposfirstshell+=1
        elif restype in negq:
            nnegfirstshell+=1
        sstype=pyrosetta.rosetta.core.scoring.dssp.Dssp(p).get_dssp_secstruct(resnum)
        renergy=float(p.energies().residue_total_energy(resnum))
        data_each_res.append((resnum,restype,sstype,renergy,res_lig_score))
        #identifying secondary residues
        curr_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(resnum)
        curr_res_neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(curr_res_selector, 8., False)
        curr_res_neighborhood_selector_bool = curr_res_neighborhood_selector.apply(p)
        curr_res_neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(curr_res_neighborhood_selector_bool)
        curr_res_contact_res=list(curr_res_neighborhood_residues_resnums)
        n_sec_res=len(curr_res_contact_res)
        if p.total_residue() in curr_res_contact_res:
            curr_res_contact_res.remove(p.total_residue())
        if resnum in curr_res_contact_res:
            curr_res_contact_res.remove(resnum)
    else:
        first_shell_res.remove(resnum)
        continue
n_first_shell_res=len(first_shell_res)
bscharge=nposfirstshell-nnegfirstshell
#now identify the network of hbonds about the ligand and analyze
network_res=[]
for i in all_sec_res:
    network_res.append(i)
for i in first_shell_res:
    network_res.append(i)
network_res=set(network_res)
network_res=sorted(network_res)
network_pose=Pose()
ligand_pose=p.residue(p.total_residue()).clone()
network_pose.append_residue_by_jump(ligand_pose, 1)
res_pose1=p.residue(network_res[0]).clone()
network_pose.append_residue_by_jump(res_pose1, 1)
c=1
for i,resnum in enumerate(network_res[1:]):
    res_pose=p.residue(resnum).clone()
    last_num=network_res[i-1]
    if last_num==resnum-1:
        network_pose.append_residue_by_jump(res_pose, c)
    else:
        c+=1
        network_pose.append_residue_by_jump(res_pose, c)
hbond_set = rosetta.core.scoring.hbonds.HBondSet()
network_pose.update_residue_neighbors()
rosetta.core.scoring.hbonds.fill_hbond_set(network_pose, False, hbond_set)
# rosetta.core.scoring.hbonds.fill_hbond_set_by_AHdist_threshold(network_pose,3.0,hbond_set)
s=network_pose.sequence()
hbond_data={}
if hbond_set.nhbonds()>0:
    for hbond_index in range(1,hbond_set.nhbonds()+1):
        drip=hbond_set.hbond(hbond_index).don_res_is_protein()
        arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
        donres_ind=int(hbond_set.hbond(hbond_index).don_res())
        accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
        donres=s[donres_ind-1]
        accres=s[accres_ind-1]
        acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
        donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
        hbdist=float(hbond_set.hbond(hbond_index).get_HAdist(network_pose))
        don_atom_index=int(network_pose.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
        acc_atom_data=str(network_pose.residue(accres_ind).atom_type(acc_atom_index))
        acc_atom=(acc_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
        don_atom_data=str(network_pose.residue(donres_ind).atom_type(don_atom_index))
        don_atom=(don_atom_data.split('\n'))[0].split('Atom Type:')[1].strip()
        acc_atom_name=str(network_pose.residue(accres_ind).atom_name(acc_atom_index)).strip(' ')
        don_atom_name=str(network_pose.residue(donres_ind).atom_name(don_atom_index)).strip(' ')
        hbond_data[hbond_index]=[donres_ind,accres_ind,donres,accres,drip,arip,don_atom,acc_atom,hbdist,don_atom_name,acc_atom_name]
#identify the network of hbonds about the ligand
#that is, residues hbonded with ligand, and other residues to which those residues are hbonded
network_members=[]
lig_atom_type_hbonds=[]
hb_w_lig=[]
ptn_lig_hb=[]
for keyb in hbond_data.keys():
    l=hbond_data[keyb]
    if l[0]==1:
        hb_w_lig.append(l[1])
        network_members.append((l[0],l[1],l[-2],l[-1]))
        ptn_lig_hb.append((l[-2],l[-3],'ligdon'))
        lig_atom_type_hbonds.append((l[6],l[7],l[8]))
    if l[1]==1:
        hb_w_lig.append(l[0])
        network_members.append((l[0],l[1],l[-2],l[-1]))
        ptn_lig_hb.append((l[-1],l[-4],'ligacc'))
        lig_atom_type_hbonds.append((l[7],l[6],l[8]))
for i in hb_w_lig:
    hbws=[]
    for keyc in hbond_data.keys():
        l=hbond_data[keyc]
        if l[0]==i:
            if l[1]!=1:
                hbws.append(l[1])
                network_members.append((l[0],l[1],l[-2],l[-1]))
        if l[1]==i:
            if l[0]!=1:
                hbws.append(l[0])
                network_members.append((l[0],l[1],l[-2],l[-1]))
data_across_strc[sta]=[pdbid,ligname,naroatoms,netq,narofirstshell,
                       nposfirstshell,nnegfirstshell,bscharge,npotentialhb,
                       alltypes,
                       data_each_res,lig_atom_type_hbonds,network_members]
# data_across_strc[sta]=
# [n_first_shell_res,bse,nethbe,data_each_res, network_members, ptn_lig_hb]
except:
print('\n\n\nproblem analyzing '+str(sta))
continue
'''
