#analysis of scores from json file
sfname='a8sr5shorthbfd1.json'


import json
scores2 = [json.loads(line) for line in open(sfname,'r')]
terms=list(scores[0].keys())
#make a pdf showing all score distributions
######
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
def plot_dists(terms,scores,outfilename):
    pdf = PdfPages(outfilename)
    for term in terms:
        allscores=[]#score
        if term != 'decoy':
            for d in scores:
                allscores.append(float(d[term]))
        fig,ax=plt.subplots()
        if len(allscores)!=0:
            ax.hist(allscores)#bins=int(len(allscores)/20)
            ax.set_title(term)
            ax.set_ylabel('frequency')
            ax.set_xlabel('score')
            pdf.savefig()
            plt.clf()
    pdf.close()

plot_dists(terms,scores,'a8sr5shorthbfd1.json')

def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        termval=float(d[term])
        if condition=='<':
            if termval<=float(threshold):
                filtered_scores.append(d)
            else:
                pass
        elif condition=='>':
            if termval>=float(threshold):
                filtered_scores.append(d)
            else:
                pass
    return filtered_scores

f1=return_filtered(scores,'buns2interface','<',7.0)
plot_dists(terms,f1,'a8sr5shorthbfd1.pdf')

f=open('a8sr5shorthbfd1.json','w')
for d in scores:
    f.write(str(d))
f.close()

f1=return_filtered(scores,'buns2interface','<',0.0)
85
f2=return_filtered(f1,'contact_molsurf','>',150)
In [9]: len(f2)
Out[9]: 4
f2=return_filtered(f1,'lighphobesasa','<',50)

filtered_strc=[]
for d in f2:
    filtered_strc.append(d['decoy'])
['UM_1_I42G43W159_1_clean_1xm3B_209_1fd1_0001',
 'UM_1_T68S65S33_1_model_104915_361_1fd1_0001',
 'UM_1_T30S32S50_1_5tpj_128385_362_1fd1_0001',
 'UM_1_T94S96S19_1_5tpj_117768_361_1fd1_0001',
 'UM_1_T83S96S19_1_5tpj_117768_362_1fd1_0001']

os.makedirs('filtered',exist_ok=True)
for i in filtered_strc:
    os.system('cp '+i+'.pdb filtered/'+i+'.pdb')


scp cgalvin@log2.wynton.ucsf.edu:a8sr5matches/a8sshorthb/design/a8sr5fd1/a8sr5shorthbfd1_dists.pdf ~/desktop/a8sr5shorthbfd1_dists.pdf
scp -r cgalvin@log2.wynton.ucsf.edu:a8sr5matches/a8sshorthb/design/a8sr5fd1/filtered/ ~/desktop/a8sr5fd1filtered
