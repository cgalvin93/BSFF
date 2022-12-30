#just execute in working directory (target parent w/ np fuzzball)
#time python ~/BSFF/check_fuzzball_vectors.py


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
    if calc==0.0:
        return True
    else:
        return False
badres=[]
# c=1
for i,e in enumerate(terindices[:-1]):
    currentresatoms=lines[e+1:terindices[i+1]]
    atomcoords=[]
    # print(c)
    # c+=1
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
                    if p1==p2:
                        badres.append((e+1,terindices[i+1]))
                    else:
                        if p1==p3:
                            badres.append((e+1,terindices[i+1]))
                        else:
                            if p2==p3:
                                badres.append((e+1,terindices[i+1]))
                            else:
                                if a(p1,p2,p3):
                                    # print('caught colinear')
                                    # print(p1)
                                    # print(p2)
                                    # print(p3)
                                    badres.append((e+1,terindices[i+1]))
                                else:
                                    pass

sbr=list(set(badres))
print('this many bad res caught:')
print(str(len(sbr)))
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
