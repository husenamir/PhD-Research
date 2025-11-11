#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 18:18:16 2020

@author: jamunoz
"""

import yaml, sys
import pandas as pd
import seaborn as sns
import matplotlib as plt
from matplotlib.pyplot import plot
import sklearn.mixture
import numpy as np
from scipy.optimize import curve_fit

root = 'B2_2.78_1500K'
#in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/'+root+'/'
in_path = '/Users/csloaner/Documents/FeV_lammps/phonopy/'+root+'/'

# min_e = 0 #THz
# max_e = 15.0 #THz
# nbins = 100

# min_e = -3.0 #THz
# max_e = 15.0 #THz
# nbins = 120

min_e = -3.5 #THz
max_e = 16.5 #THz
nbins = 130

def do_dispersion_df(df_hsd, min_e, max_e, nbins):
    range_e = max_e - min_e
    heatmap_dict = {}
    for j in range(201):
        hit_dict = {}
        for b in range(1, nbins+1):
            bound = range_e/nbins * b + min_e
            hit_dict[bound] = 0
        for i in df_hsd.iloc[:,j].values:
            for key in hit_dict.keys():
                if i < key:
                    hit_dict[key] = hit_dict[key] + 1
                    break
        
        heatmap_dict[j] = pd.Series(hit_dict)
        
    df = pd.DataFrame(heatmap_dict)
    df = df.sort_index(axis=0 ,ascending=False)
    
    return df

print(in_path)

# df_hsd  = pd.read_csv(in_path+root+'_G_X_500_100.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# heatmap = sns.heatmap(df, vmin=0, vmax=100, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_G_X_500_100.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_X_M_500_100.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# heatmap = sns.heatmap(df, vmin=0, vmax=100, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_X_M_500_100.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_M_G_500_100.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# heatmap = sns.heatmap(df, vmin=0, vmax=100, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_M_G_500_100.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_G_R_500_100.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# heatmap = sns.heatmap(df, vmin=0, vmax=100, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_G_R_500_100.png', dpi=300)

df_hsd  = pd.read_csv(in_path+root+'_M_G_3000_999.csv', index_col=0) # hsd = high symmetry direction
df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
df.T.iloc[0,:].to_csv(in_path+root+'_M_3000_999_intensity.csv')
ax1 = df.T.iloc[0,:].plot()
ax1.set_ylim(0,1000)
# df.T.iloc[0,70:120].plot()

sys.exit()

hist = df.T.iloc[0,70:120].values
bin_edges = df.T.iloc[0,70:120].index
bin_centres = bin_edges

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


p0 = [1., 0., 1.]
coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

hist_fit = gauss(bin_centres, *coeff)

print('Fitted mean = ', coeff[1])
print('Fitted standard deviation = ', coeff[2])

plot(bin_centres, hist, label='Test data')
plot(bin_centres, hist_fit, label='Fitted data')

# print(df.T.iloc[0])

#df_hsd.iloc[0].plot()


sys.exit()

######
# 2.90, 300K; mu = 4.28 THz, std = 0.55 THz
# 
#
######

# _G_X
# _X_M
# _M_G
# _G_R