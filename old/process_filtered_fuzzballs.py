import os
import sys
#
fogpath=sys.argv[1]
ligand_path=sys.argv[2]
target_name=sys.argv[3]
frag_name=sys.argv[4]
#
fog=open(fogpath,'r')
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
offbname=('_').join([target_name,frag_name,'polar_contact_fuzzball.pdb'])
of=open(offbname,'w')
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
offbname2=('_').join([target_name,frag_name,'nonpolar_contact_fuzzball.pdb'])
of=open(offbname2,'w')
for line in clean_fuzzball_lines:
    of.write(line)
of.close()
