#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 18:18:16 2020

@author: jamunoz
"""

import yaml, sys, os
import pandas as pd

init_file_number = 0
ndisps = 1
root = 'B2_2.86_10K'
# in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/' + root + '/'
#in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/' + root + '_7x7x7/'
in_path = 'D:\\Research\\phonopy\\'+root+'\\'
suffix = '_' + str(init_file_number) + '_' + str(ndisps)
directions = ['_G_X', '_X_M', '_M_G', '_G_R', '_R_M']

# dones = []
# files = os.listdir(in_path)
# for file in files:
#     if '.yaml' in file:
#         # print(file, file[5:-5])
#         dones.append(int(file[5:-5]))

# kill = False        
# for file_number in range(init_file_number, init_file_number+ndisps):
#     if file_number in dones:
#         continue
#     else:
#         kill = True
#         print(file_number)
    
# if kill:
#     print('Per your request, this was killed')
#     sys.exit()
    
# sys.exit()

out_filenames = []
print(in_path)
for direction in directions:
    out_filenames.append(direction+suffix+'.csv')


freq_dict = {}  

freq_dict_G_X = {}
freq_dict_X_M = {}
freq_dict_M_G = {}
freq_dict_G_R = {}  
freq_dict_R_M = {}  

for file_number in range(init_file_number, init_file_number+ndisps):    
    in_filename = 'band_' + str(file_number) + '.yaml'
    print(in_filename)
    with open(in_path+in_filename) as file:
        yaml_file = yaml.full_load(file)
    
    phonons = yaml_file['phonon']
    
    for i in range(201):
        distance = phonons[0]['distance']
        band_contents = phonons[i]['band']
        freqs = []
        for j in range(6):
            freqs.append(band_contents[0]['frequency'])
            
        if file_number == init_file_number:
            freq_dict_G_X[distance] = freqs
        else:
            freq_dict_G_X[distance].extend(freqs)
            
    for i in range(201,402):
        distance = phonons[i]['distance']
        band_contents = phonons[i]['band']
        freqs = []
        for j in range(6):
            freqs.append(band_contents[j]['frequency'])
            
        if file_number == init_file_number:
            freq_dict_X_M[distance] = freqs
        else:
            freq_dict_X_M[distance].extend(freqs)
            
    for i in range(402,603):
        distance = phonons[i]['distance']
        band_contents = phonons[i]['band']
        freqs = []
        for j in range(6):
            freqs.append(band_contents[j]['frequency'])
            
        if file_number == init_file_number:
            freq_dict_M_G[distance] = freqs
        else:
            freq_dict_M_G[distance].extend(freqs)
            
    for i in range(603,804):
        distance = phonons[i]['distance']
        band_contents = phonons[i]['band']
        freqs = []
        for j in range(6):
            freqs.append(band_contents[j]['frequency'])
            
        if file_number == init_file_number:
            freq_dict_G_R[distance] = freqs
        else:
            freq_dict_G_R[distance].extend(freqs)
            
    for i in range(804,1005):
        distance = phonons[i]['distance']
        band_contents = phonons[i]['band']
        freqs = []
        for j in range(6):
            freqs.append(band_contents[j]['frequency'])
            
        if file_number == init_file_number:
            freq_dict_R_M[distance] = freqs
        else:
            freq_dict_R_M[distance].extend(freqs)
            
    print(file_number)

df_G_X = pd.DataFrame(freq_dict_G_X)
df_X_M = pd.DataFrame(freq_dict_X_M)
df_M_G = pd.DataFrame(freq_dict_M_G)
df_G_R = pd.DataFrame(freq_dict_G_R)
df_R_M = pd.DataFrame(freq_dict_R_M)

freq_dfs = [df_G_X, df_X_M, df_M_G, df_G_R, df_R_M]

for c, out_filename in enumerate(out_filenames):
    freq_dfs[c].to_csv(in_path+root+out_filename)
    
# sys.exit()
    

# df_G_X.to_csv(in_path+root+'_G_X_3000_100.csv')
# df_X_M.to_csv(in_path+root+'_X_M_3000_100.csv')
# df_M_G.to_csv(in_path+root+'_M_G_3000_100.csv')
# df_G_R.to_csv(in_path+root+'_G_R_3000_100.csv')
# df_R_M.to_csv(in_path+root+'_R_M_3000_100.csv')