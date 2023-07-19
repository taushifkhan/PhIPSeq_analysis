import pandas as pd
import gzip
import glob
import os

def nhits(listFile):
    k = pd.DataFrame()
    for i in listFile:
        if k.shape[1]==0:
            k = pd.read_csv(gzip.open(i)).set_index('id')
        else:
            _ = pd.read_csv(gzip.open(i)).set_index('id')
            k = k.join(_)
    print (k.shape,)
    return k.sum(axis=1)

beads = []
samples_qbb = []
sample_occ = []
sample_mis = []
for k in glob.glob('../phipDir/*.zhit.csv.gz'):
    if k[10] == 'B':
        beads.append(k)
    elif k[10:13] == 'SIR':
        samples_qbb.append(k)
    elif k[10:13] =='OCC':
        sample_occ.append(k)
    else:
       sample_mis.append(k)

print ("beads : {}\t qbb: {}\t occ: {}\trest: {}".format(len(beads),len(samples_qbb),len(sample_occ),len(sample_mis)))
 

k_bead = nhits(beads)
print (k_bead[k_bead>2].shape)
k_bead.to_csv('./nhits_beads.csv')
os.system("gzip -c ./nhits_beads.csv > ./nhits_beads.csv.gz")

k_samples_qbb = nhits(samples_qbb)
print (k_samples_qbb[k_samples_qbb<2].shape)
k_samples_qbb.to_csv('./nhits_samples_qbb.csv')
os.system("gzip -c ./nhits_samples_qbb.csv > ./nhits_sample_qbb.csv.gz")


k_samples_occ = nhits(sample_occ)
print (k_samples_occ[k_samples_occ<2].shape)
k_samples_occ.to_csv('./nhits_samples_occ.csv')

k_samples_mis = nhits(sample_mis)
print (k_samples_mis[k_sample_mis<2].shape)
k_samples_mis.to_csv('./nhits_samples_mis.csv')
#k.sum(axis=1).to_csv('allbeads.csv')
