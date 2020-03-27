#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 18:09:13 2020

@author: maksbh
"""

import numpy as np
import h5py



#TODO: Edit the voxel properties
num_materials = 2
Nx = 10
Ny = 10
Nz = 10
filename = "mytestfile.h5"

#Allocating space. Do not change this block
size = Nx*Ny*Nz
Allignmentdtype = [('Sx', (np.float32)), ('Sy', np.float32), ('Sz', np.float32)]
volume_fraction = np.zeros(((num_materials,size)),dtype = np.float32)
allignment_vector = np.zeros(((num_materials,size)),dtype = Allignmentdtype)
unaligned_fraction = np.zeros(((num_materials,size)),dtype = np.float32)



## TODO: Edit enteries of volume fractions.
for i in range(0,num_materials):
    volume_fraction[i,:] = 0


## TODO: Edit enteries of alignment vectors.
for i in range(0,num_materials):
    allignment_vector[i]['Sx'] = 1  
    allignment_vector[i]['Sy'] = 1
    allignment_vector[i]['Sz'] = 1

## TODO: Edit enteries of unaligned fractions.
for i in range(0,num_materials):
    unaligned_fraction[i] = 1





#Do not change below this.
f = h5py.File(filename, "w")
g0 = f.create_group('igor_parameters')
dset = g0.create_dataset("igormaterialnum", data = num_materials,dtype = np.int32)
g1 = f.create_group('vector_morphology')
for i in range(0,num_materials):
    name = "Mat_"+ str(i+1) + "_unaligned"
    unalignedFraction = np.reshape(volume_fraction[i]*unaligned_fraction[i],(Nx,Ny,Nz))
    dset = g1.create_dataset(name, data = (unalignedFraction),dtype = np.float32)
    
AlignedFraction = np.zeros((size,3),dtype = np.float32)
for i in range(0,num_materials):
    name = "Mat_"+ str(i+1) + "_alignment"
    
    for k in range(size):
        AlignedFraction[k][0] = volume_fraction[i][k]*allignment_vector[i][k]["Sx"]
        AlignedFraction[k][1] = volume_fraction[i][k]*allignment_vector[i][k]["Sy"]
        AlignedFraction[k][2] = volume_fraction[i][k]*allignment_vector[i][k]["Sz"]
    
    
    dset = g1.create_dataset(name, data = (AlignedFraction))
    
    
f.close()