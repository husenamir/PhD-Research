# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 22:30:42 2023

@author: ahusen
"""
#import seaborn as sb
#import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, isclose
#import h5py
import os, sys, yaml
import pandas as pd
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from geneticalgorithm import geneticalgorithm as ga


# parameter setting
root = 'B2_2.86_10K'
n = 250 # Number of atoms
ndisps = 5 # Number of dispersions to generate (equal to number of yalm files)
lattice_value = 2.86
system_size = 5
alat = lattice_value * system_size # the size of the box
total_atoms = 2 * (system_size ** 3)
init_no = 0


# Defining unit cell

unit_cell_length = lattice_value #lattice_parameter#*2

basis = np.eye(3)*unit_cell_length

base_atoms = np.array([[0.0, 0.0, 0.0],
                       [0.5, 0.5, 0.5]])*unit_cell_length

# Generate atom positions
in_path = '/Users/amirhusen/Desktop/Amir/Research/pressure_runs/'+root+'/'
in_filename = 'ideal_FeV_'+root+'_unit.txt' # Modify pos file
ide_lat = np.genfromtxt(in_path+in_filename)
#print(ide_lat.shape)
ide_lat = ide_lat[:n, 1:]
a_val = alat / ((n / 2.)**(1./3.)) #Side length of each cube in bcc

ideal_distances = np.zeros((n, 3, n))
ideal_dist_sca = np.zeros((n, n)) # sca stands for scalar

for i in range(n):
    for j in range(n):
        if i != j:
            ideal_distances[i, :, j] = ide_lat[j, :] - ide_lat[i, :] #Ideal distance vector from i to j

            #Apply mic: 
            for d in range(3):
                if ideal_distances[i, d, j] > (alat / 2):
                    ideal_distances[i, d, j] -= alat

                elif ideal_distances[i, d, j] <= (-alat / 2):
                    ideal_distances[i, d, j] += alat

            ideal_dist_sca[i, j] = np.linalg.norm(ideal_distances[i, :, j]) #Distance scalar

neighbors = [] #List of dictionaries holding a list of the IDs for all types of nn

#Distances to each type of neighbor:
first = (sqrt(3) * a_val) / 2
second = a_val
third = sqrt(2) * a_val
fourth = (sqrt(11) * a_val) / 2
fifth = sqrt(3) * a_val
nn_dist = [first, second, third, fourth, fifth]

for i in range(n):
    i_neigh = [[i]] 
    for prox in range(5):
        i_neigh.append([]) #Empty placeholder lists for neighbors of i
    
    #Add each list of nn to the placeholder dictionary, append the dictionary to the list:
    for j in range(n):
        if i != j:
            for prox in range(1,6):
                if isclose(ideal_dist_sca[i, j], nn_dist[prox-1], rel_tol=.001):
                    i_neigh[prox].append(j)

    neighbors.append(i_neigh)

fc_mat = np.zeros((n, n, 3, 3))
base_mat = np.zeros((6, 3, 3))
nn_matrix = .5 * a_val * np.array([[0.,0.,0.],
                                   [1.,1.,1.],
                                   [2.,0.,0.],
                                   [0.,2.,2.],
                                   [3.,1.,1.],
                                   [2.,2.,2.]])

# print(nn_matrix)
# sys.exit()
fc_arr = np.array([[0, 0, 14, 14],# alfa, beta, gama, delta
                   [1, 1, 2, 2],
                   [3, 4, 14, 14],
                   [5, 6, 14, 7],
                   [8, 9, 10, 11],
                   [12, 12, 13, 13]])


def f(X):
    # α0 + 8α1 + 2α2 + 4β2 + 4α3 + 8β3 + 8α4 + 16β4 + 8α5 = 0 and α1 = X[0] = [-0.85, -0.55]
    α0 = X[0]
    α1 = X[1]
    β1 = -0.9375
    α2 = -0.650625
    β2 = -0.245625 
    α3 = 0.0535625 
    β3 = -0.0975
    γ3 = -0.06625
    α4 = -0.0473125
    β4 = -0.0011875
    γ4 = -0.0015
    δ4 = -0.03025
    α5 = 0.018125
    β5 = -0.02
    #α0 =  (8*α1 + 2*α2 + 4*β2 + 4*α3 + 8*β3 + 8*α4 + 16*β4 + 8*α5)
    #α0 = -(8*α1 + 2*α2 + 4*β2 + 4*α3 + 8*β3 + 8*α4 + 16*β4 + 8*α5)
    fc = [α0, α1, β1, α2, β2, α3, β3, γ3, α4, β4, γ4, δ4, α5, β5, 0.0]
    
    for prox in range(6):
        base_mat[prox, :, :] = np.array([[fc[fc_arr[prox, 0]], fc[fc_arr[prox, 2]], fc[fc_arr[prox, 2]]],
                                         [fc[fc_arr[prox, 2]], fc[fc_arr[prox, 1]], fc[fc_arr[prox, 3]]],
                                         [fc[fc_arr[prox, 2]], fc[fc_arr[prox, 3]], fc[fc_arr[prox, 1]]]])
    
    for i in range(250):
        for prox in range(6):
            for j in neighbors[i][prox]:
                fc_mat[i, j, :, :] = base_mat[prox, :, :]
                if prox in [2, 3, 4]:
                    for r in range(1,3):
                        if isclose(abs(ideal_distances[i, r, j]), nn_matrix[prox, 0], rel_tol=.001):
                            fc_mat[i, j, [r, 0], :] = fc_mat[i, j, [0, r], :]
                            fc_mat[i, j, :, [r, 0]] = fc_mat[i, j, :, [0, r]]
                            fc_mat[i, j, :, [r, 0]] = 0
                            
                for r in range(3):
                    if ideal_distances[i, r, j] < 0:
                        fc_mat[i, j, r, :] = -fc_mat[i, j, r, :]
                        fc_mat[i, j, :, r] = -fc_mat[i, j, :, r]
                        
    np.around(fc_mat, decimals=5)
    
    
    out_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
    out_filename = 'FORCE_CONSTANTS_GA'
    with open(out_path+out_filename, 'w') as file:
        print(str(n) + ' ' + str(n), file=file)
        for i in range(n):
            for j in range(n):
                print(str(i+1) + ' ' + str(j+1), file=file)
                np.savetxt(file, fc_mat[i, j, :, :])
    
    
    #out_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/B2_2.86_10K'
    os.chdir(out_path)
    command = 'phonopy -p -s band.conf'
    
    # print('here 1', m)
    os.system(command)
    # print('here 2', m)
    command = 'mv band.yaml band_' + 'ga' + '.yaml'
    os.system(command)
    command = 'mv band.pdf band_' + 'ga' + '.pdf'
    os.system(command)
    command = 'rm phonopy.yaml'
    os.system(command)

    
    
    
    """
    in_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
    path = [[[0, 0, 0], [0.5, 0, 0.0], [0.5, 0.5, 0.0],
            [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0]]]
    labels = ["$\\Gamma$", "X", "M",  "$\\Gamma$", "R", "M"]
    qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=201)
    phonons = phonopy.load(supercell_matrix=[5, 5, 5],
                  primitive_matrix='auto',
                  unitcell_filename= in_path+"POSCAR",
                  force_constants_filename= in_path+"FORCE_CONSTANTS_GA")
    phonons.run_band_structure(qpoints, path_connections=connections, labels=labels)
    #phonons.plot_band_structure().show()
    #sys.exit()
    res=phonons.get_band_structure_dict()['frequencies']
    #print(len(res))
    freqs=[]
    for i in range(5):
        for j in range(201):
            freqs.append(res[i][j][0])
    """
            
            
            
    in_path_y = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
    in_filename = 'band_' + 'ga' + '.yaml'
    #print(in_filename)
    with open(in_path_y+in_filename) as file:
        yaml_file = yaml.full_load(file)
    phonons = yaml_file['phonon']
    phonon_no = len(phonons)
    band_contents = []
    freqs = []
    for j in range(phonon_no):
        band_contents.append(phonons[j]['band'])
        freqs.append(band_contents[j][0]['frequency'])
            
    # Taking negative frequencies
    neg_freqs = []
    for i in range(len(freqs)):
        if freqs[i] < 0:
            neg_freqs.append(freqs[i])
    
    squared_freqs = [number**2 for number in freqs] # for all
    rmse = sqrt(sum(squared_freqs)/len(freqs))
    pen1 = 0
    total = (α0 + 8*α1 + 2*α2 + 4*β2 + 4*α3 + 8*β3 + 8*α4 + 16*β4 + 8*α5)
    if total != 0:
        pen1 = total
    pen2 = 0
    if len(neg_freqs) > 0:
        pen2 = len(neg_freqs) * 1000000
    
    obj = rmse + pen2 + pen1
    return obj 

#print(f([-0.1]))
#sys.exit()

varbound=np.array([[9, 10],[-1.00, 0.00]])

algorithm_param = {'max_num_iteration': 500,\
                   'population_size':10,\
                   'mutation_probability':0.01,\
                   'elit_ratio': 0.01,\
                   'crossover_probability': 0.5,\
                   'parents_portion': 0.3,\
                   'crossover_type':'uniform',\
                   'max_iteration_without_improv':None}

model=ga(function=f,\
            dimension=2,\
            variable_type='real',\
            function_timeout=30,\
            variable_boundaries=varbound,\
            algorithm_parameters=algorithm_param)


model.run()