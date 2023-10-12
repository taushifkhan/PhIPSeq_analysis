# Loading the required packages

import re
import sys
import functools
import collections

import numpy as np
import scipy as sp
import scipy.stats
import scipy.optimize
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab

import gzip

#load annotation file
annotations = pd.read_csv('VIR3_clean.csv.gz', index_col=0)
scores_sidra = pd.read_csv('scores.csv', index_col=0)
public_epitopes = pd.read_csv('public_epitope_annotations.csv', index_col=0)
public_epitopes_clean = public_epitopes.ix[public_epitopes['fraction_positive'] > 0.3]
hits_sidra = pd.read_csv('hits.csv', index_col=0)
thresholds = pd.read_csv('virscan-thresholds-spp-sh.csv', index_col=0)

#set up dataframe

public_epitope_df = pd.DataFrame(index = set(public_epitopes['species']), columns = hits_sidra.columns)
for spe in public_epitope_df.index:
    temp = public_epitopes.ix[public_epitopes['species'] == spe]
    public_epitope_df.ix[spe] = hits_sidra.ix[temp.index].sum(axis = 0)

public_epitope_clean_df = pd.DataFrame(index = set(public_epitopes_clean['species']), columns = hits_sidra.columns)
for spe in public_epitope_clean_df.index:
    temp = public_epitopes_clean.ix[public_epitopes_clean['species'] == spe]
    public_epitope_clean_df.ix[spe] = hits_sidra.ix[temp.index].sum(axis = 0)

#incorporate criteria for public epitode and virus threshold    
vir2_calls_public = (scores_sidra.ix[thresholds.index].gt(thresholds['threshold'], axis = 0)) & (public_epitope_df > 0)
vir2_calls_public_clean = (scores_sidra.ix[thresholds.index].gt(thresholds['threshold'], axis = 0)) & (public_epitope_clean_df > 0)


#listing the results
vir2_calls_public.sum(axis = 1).sort_values(ascending = False)[0:200] / len(vir2_calls_public.columns) *100
vir2_calls_public_clean.sum(axis = 1).sort_values(ascending = False)[0:200] / len(vir2_calls_public.columns) *100