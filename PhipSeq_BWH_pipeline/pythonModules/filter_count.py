import pandas as pd
import sys

if not (sys.argv[1] and sys.argv[2]):
    print ("give arguments : sampleMetadata, sample statistics") 
    sys.exit()

p_corr_ub = 0.6
p_corr_lb = 0
sm = pd.read_csv(sys.argv[1])
fl = pd.read_csv(sys.argv[2])
print ("input samples",fl.shape[0])

print (fl.head())
print (sm.head())
#import ipdb; ipdb.set_trace()
p_corr_lb = fl.pearson_coff.min()
print ("UB: {} LB: {}".format(p_corr_ub, p_corr_lb))

ok_sample = fl[(fl.pearson_coff < p_corr_ub)&(fl.pearson_coff >= p_corr_lb)].SampleID.values
sm = sm.set_index('sid')
print (sm.loc[ok_sample].head())
print (len(ok_sample))

outsample = "{}_{}_{}.csv".format(sys.argv[1],p_corr_lb,p_corr_ub)
print ("output written in :",outsample)
sm.loc[ok_sample].to_csv(outsample)
