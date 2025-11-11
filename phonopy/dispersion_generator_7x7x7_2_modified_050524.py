# Date: 07/21/2023
# Amir Husen
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, isclose
import h5py
import os, sys, yaml
import pandas as pd
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections


#######

# init_file_number = 3000
# ndisps = 500
# root = 'B2_2.88_300K'
# in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/phonopy/' + root + '/'

# dones = []
# files = os.listdir(in_path)
# for file in files:
#     if '.yaml' in file:
#         # print(file, file[5:-5])
#         dones.append(int(file[5:-5]))

# left = []     
# for file_number in range(init_file_number, init_file_number+ndisps):
#     if file_number in dones:
#         continue
#     else:
#         left.append(file_number)
    
    
# print(left)
# 
####### C:\MyFiles\pressure_runs\B2_2.86_10K


# Remember to modify POSCAR file to correct lattice parameter
# Setup and Created Cristal 
root = 'B2_2.86_10K'
n = 250 # Number of atoms
ndisps = 5 # Number of dispersions to generate (equal to number of yalm files)
lattice_value = 2.86
system_size = 5
alat = lattice_value * system_size # the size of the box
total_atoms = 2 * (system_size ** 3)
init_no = 0

"""
lattice_parameters = []
for i in range(1):
    lattice_parameters.append(lattice_value_for)
    lattice_value_for += 0.01
"""
#lattice_parameters = [lattice_value]

#%% Calculates the directory for each lattice parameter, creates the atoms' 
#   positions, and runs the simulations 

#for lattice_parameter in lattice_parameters:
#%lattice_parameters has single element, so skip loop%
# Unit cell length for FeV

""""
for lattice_parameter in lattice_parameters:

    # Unit cell length for FeV
    unit_cell_length = lattice_parameter#*2
    
    # Defining unit cell
    basis = np.array([[1.0, 0.0, 0.0], 
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])*unit_cell_length
"""
unit_cell_length = lattice_value #lattice_parameter#*2

# Defining unit cell
basis = np.eye(3)*unit_cell_length

base_atoms = np.array([[0.0, 0.0, 0.0],
                       [0.5, 0.5, 0.5]])*unit_cell_length

in_path = '/Users/amirhusen/Desktop/Amir/Research/pressure_runs/'+root+'/'
in_filename = 'ideal_FeV_'+root+'_unit.txt' # Modify pos file
ide_lat = np.genfromtxt(in_path+in_filename)
#print(ide_lat.shape)
ide_lat = ide_lat[:n, 1:]
"""# Generate atom positions
positions = []
for i in range(system_size): 
    for j in range(system_size):
        for k in range(system_size):
            base_position = np.array([i, j, k]) #Local variable
            cart_position = np.inner(basis.T, base_position) #Local
            for atom in base_atoms:
                positions.append(cart_position + atom)
ide_lat = np.array(positions)"""
"""
#in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/pressure_runs/'+root
in_path = 'C:\\Research\\pressure_runs\\'+root+'\\'
in_filename = 'ideal_FeV_'+root+'_unit.txt' # Modify pos file
ide_lat = np.genfromtxt(in_path+in_filename)
#print(ide_lat.shape)
ide_lat = ide_lat[:n, 1:]

for i in range(250):
    for j in range(250):
        if np.all(ide_lat[i] == ide_lat_n[j]):
            print(i,j)

#print(ide_lat, ide_lat_n)

fig = plt.figure()
ax = fig.add_subplot(projection = '3d')
x1 = ide_lat_n[:,0]
y1 = ide_lat_n[:,1]
z1 = ide_lat_n[:,2]
ax.scatter(x1, y1, z1)
ax.view_init(elev = 0, azim = 0, roll = 0)
plt.show()
"""
#sys.exit()


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


#print(fc_arr)
#sys.exit()
# For B2
#fc_in_path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/FeV_lammps/force_constants/'+root+'_7x7x7/'
#fc_in_path = 'D:\\Research\\force_constants\\'+root+'\\'

#fc_filename = 'fc_'+root+'_4000.csv'

#fc_per_ts_df = pd.read_csv(fc_in_path+fc_filename, skiprows=2, header = None)
#Taken input manually
data = np.array([[9.262,		-0.95,	-0.9375,	-0.650625,	-0.245625,		0.0535625,	-0.0975,	-0.06625,		-0.0473125,	-0.0011875,	-0.0015,	-0.03025,		0.018125,	-0.02],
                 [9.262,		-0.85,	-0.9375,	-0.650625,	-0.245625,		0.0535625,	-0.0975,	-0.06625,		-0.0473125,	-0.0011875,	-0.0015,	-0.03025,		0.018125,	-0.02],
                 [9.262,		-0.75,	-0.9375,	-0.650625,	-0.245625,		0.0535625,	-0.0975,	-0.06625,		-0.0473125,	-0.0011875,	-0.0015,	-0.03025,		0.018125,	-0.02],
                 [9.262,		-0.65,	-0.9375,	-0.650625,	-0.245625,		0.0535625,	-0.0975,	-0.06625,		-0.0473125,	-0.0011875,	-0.0015,	-0.03025,		0.018125,	-0.02],
                 [9.262,		-0.55,	-0.9375,	-0.650625,	-0.245625,		0.0535625,	-0.0975,	-0.06625,		-0.0473125,	-0.0011875,	-0.0015,	-0.03025,		0.018125,	-0.02]])
fc_per_ts_df =pd.DataFrame(data, columns=['alfa0', 'alfa1',	'beta1',	'alfa2',	'beta2',	'alfa3', 'beta3',	'gama3',	'alfa4',	'beta4',	'gama4',	'delta4',	'alfa5',	'beta5'])

#print(data[2,2]*2+data[1,1]*2)
#sys.exit()     
      

rmses=[]


for m in range(init_no, init_no+ndisps):
# for m in left:
    for a in range(1,6):

        fc = [fc_per_ts_df.iloc[m, 0]+ a/10, fc_per_ts_df.iloc[m, 1], fc_per_ts_df.iloc[m, 2], fc_per_ts_df.iloc[m, 3], fc_per_ts_df.iloc[m, 4], fc_per_ts_df.iloc[m, 5], 
              fc_per_ts_df.iloc[m, 6], fc_per_ts_df.iloc[m, 7], fc_per_ts_df.iloc[m, 8], fc_per_ts_df.iloc[m, 9], fc_per_ts_df.iloc[m, 10], fc_per_ts_df.iloc[m, 11], 
              fc_per_ts_df.iloc[m, 12], fc_per_ts_df.iloc[m, 13], 0.0]
        
    
        
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
                                
                    for r in range(3):
                        if ideal_distances[i, r, j] < 0:
                            fc_mat[i, j, r, :] = -fc_mat[i, j, r, :]
                            fc_mat[i, j, :, r] = -fc_mat[i, j, :, r]
                          
        
     
        np.around(fc_mat, decimals=5)
        out_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
        out_filename = 'FORCE_CONSTANTS'
        with open(out_path+out_filename, 'w') as file:
            print(str(n) + ' ' + str(n), file=file)
            for i in range(n):
                for j in range(n):
                    print(str(i+1) + ' ' + str(j+1), file=file)
                    np.savetxt(file, fc_mat[i, j, :, :])
                    
            
        out_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/B2_2.86_10K'
        os.chdir(out_path)
        command = 'phonopy -p -s band.conf'
        
        # print('here 1', m)
        os.system(command)
        # print('here 2', m)
        command = 'mv band.yaml band_' + str(m) +'_'+ str(a)+ '.yaml'
        os.system(command)
        command = 'mv band.pdf band_' + str(m) + '_'+ str(a)+ '.pdf'
        os.system(command)
        command = 'rm phonopy.yaml'
        os.system(command)

    
      
#sys.exit()  

              
        in_path_y = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
        in_filename = 'band_' + str(m) +'_'+ str(a) + '.yaml'
        print(in_filename)
        with open(in_path_y+in_filename) as file:
            yaml_file = yaml.full_load(file)
        phonons = yaml_file['phonon']
        phonon_no = len(phonons)
        band_contents = []
        freqs = []
        for j in range(phonon_no):
            band_contents.append(phonons[j]['band'])
            freqs.append(band_contents[j][0]['frequency'])
            
        #print(freqs[:5])
        #sys.exit() 
        
        """
        in_path = '/Users/amirhusen/Desktop/Amir/Research/Phonopy/'+root+'/'
        path = [[[0, 0, 0], [0.5, 0, 0.0], [0.5, 0.5, 0.0],
                [0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0]]]
        labels = ["$\\Gamma$", "X", "M",  "$\\Gamma$", "R", "M"]
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=201)
        phonons = phonopy.load(supercell_matrix=[5, 5, 5],
                      primitive_matrix='auto',
                      unitcell_filename= in_path+"POSCAR",
                      force_constants_filename= in_path+"FORCE_CONSTANTS")
        phonons.run_band_structure(qpoints, path_connections=connections, labels=labels)
        phonons.plot_band_structure().show()
        #sys.exit()
        res=phonons.get_band_structure_dict()['frequencies']
        #print(len(res))
        freqs=[]
        for i in range(5):
            for j in range(201):
                freqs.append(res[i][j][0])"""
                
        """# Taking negative frequencies
        neg_freqs = []
        for i in range(len(freqs)):
            if freqs[i] < 0:
                neg_freqs.append(freqs[i])"""
        #print(len(neg_freqs))
        #sys.exit()
        """ #*freqs are printed on a text file*
        a_freqs=np.array(freqs).reshape(1005,1)
        with open(out_path+'new_res', 'w') as file:
            for i in range(1005):
                np.savetxt(file,a_freqs[i])
        """
                
        squared_freqs = [number**2 for number in freqs] #  for all
        rmse = sqrt(sum(squared_freqs)/len(freqs))
        """if len(neg_freqs)==0:
            rmse = 0 
        else:
            rmse = sqrt(sum(squared_freqs)/len(neg_freqs))"""
        rmses.append(rmse)
        
rmses_array=np.array(rmses).reshape(a,m+1)        
#print(np.array(rmses).reshape(a,m+1))    
sb.heatmap(rmses_array, annot= False, cmap="BuPu")  
plt.xlabel("m values")
plt.ylabel("a values")
plt.show()
       
sys.exit()

""" Study Links
http://phonopy.github.io/phonopy/phonopy-module.html
https://github.com/phonopy/phonopy/blob/develop/phonopy/phonon/band_structure.py
https://github.com/phonopy/phonopy/blob/develop/phonopy/harmonic/dynamical_matrix.py

"""