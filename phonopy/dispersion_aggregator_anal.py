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


root = 'B2_2.88_1500K'
# in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/'+root+'/'
in_path = '/Users/csloaner/Documents/FeV_lammps/phonopy/'+root+'/'
# in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/'+root+'_7x7x7/'

# min_e = 0 #THz
# max_e = 15.0 #THz
# nbins = 100

min_e = -4.5 #THz
max_e = 16.5 #THz
nbins = 140

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

# df_hsd  = pd.read_csv(in_path+root+'_G_X_3000_999.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# df.T.sort_index(axis=1, ascending=True).to_csv(in_path+root+'_G_X_3000_999_heatmap.csv')
# heatmap = sns.heatmap(df, vmin=0, vmax=800, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_G_X_3000_999.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_X_M_3000_999.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# df.T.sort_index(axis=1, ascending=True).to_csv(in_path+root+'_X_M_3000_999_heatmap.csv')
# heatmap = sns.heatmap(df, vmin=0, vmax=800, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_X_M_3000_999.png', dpi=300)

df_hsd  = pd.read_csv(in_path+root+'_M_G_3000_999.csv', index_col=0) # hsd = high symmetry direction
df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
df.T.sort_index(axis=1, ascending=True).to_csv(in_path+root+'_M_G_3000_999_heatmap.csv')
heatmap = sns.heatmap(df, vmin=0, vmax=800, cbar=False, xticklabels=False, yticklabels=10)
fig = heatmap.get_figure()
fig.savefig(in_path+root+'_M_G_3000_999.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_G_R_3000_999.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# df.T.sort_index(axis=1, ascending=True).to_csv(in_path+root+'_G_R_3000_999_heatmap.csv')
# heatmap = sns.heatmap(df, vmin=0, vmax=800, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_G_R_3000_999.png', dpi=300)

# df_hsd  = pd.read_csv(in_path+root+'_R_M_3000_999.csv', index_col=0) # hsd = high symmetry direction
# df = do_dispersion_df(df_hsd=df_hsd, min_e=min_e, max_e=max_e, nbins=nbins)
# df.T.sort_index(axis=1, ascending=True).to_csv(in_path+root+'_R_M_3000_999_heatmap.csv')
# heatmap = sns.heatmap(df, vmin=0, vmax=800, cbar=False, xticklabels=False, yticklabels=10)
# fig = heatmap.get_figure()
# fig.savefig(in_path+root+'_R_M_3000_999.png', dpi=300)


sys.exit()

# _G_X
# _X_M
# _M_G
# _G_R