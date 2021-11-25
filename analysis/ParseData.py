# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:02:26 2021

@author: Gabriel Carvalho
"""

import os
import numpy as np
from astropy.io import fits
from pathlib import Path

def generate_fits(data, filename, dim):
    # Convert array to correct format
    positions = np.empty([dim, data.shape[0], int(data.shape[1] / dim)])
    for i in range(int(data.shape[1] / dim)):
        for j in range(dim):
            positions[j,:,i] = data[:,dim*i+j]
    
    # Remove old file
    try:
        os.remove(filename)
    except OSError:
        pass
    
    if os.path.exists(filename):
        print('Failed to remove old .fits file')
        return
    
    # Make a new FITS file
    hdu = fits.PrimaryHDU(positions)
    hdu.writeto(filename)
    
    print("Exported " + filename)
    
def generate_pointcloud_xyz(data, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    for i in range(data.shape[0]):
        filename = directory + '{:04d}'.format(i) + '.xyz'
        with open(filename, 'w') as f:
            for j in range(int(data.shape[1] / 3)):
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,3*j]) + "," + str(data[i,3*j+2]) + "," + str(data[i,3*j+1]) + '\n')
                
    print("Exported pointcloud to " + directory)
    
def generate_pointcloud_xyzw(data, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    for i in range(data.shape[0]):
        filename = directory + '{:04d}'.format(i) + '.xyz'
        with open(filename, 'w') as f:
            for j in range(int(data.shape[1] / 4)):
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,4*j]) + "," + str(data[i,4*j+2]) + "," + str(data[i,4*j+1]) + "," + str(data[i,4*j+3]) + '\n')
                
    print("Exported pointcloud to " + directory)

def main():
    
    dataDM = np.loadtxt('output/dataDM.np', float, delimiter=',')
    dataGas = np.loadtxt('output/dataGas.np', float, delimiter=',')
    
    if dataDM.shape[0] > 0: generate_fits(dataDM, 'output/dataDM.fits', 3)
    if dataGas.shape[0] > 0: generate_fits(dataGas, 'output/dataGas.fits', 4)
    
    if dataDM.shape[0] > 0: generate_pointcloud_xyz(dataDM, 'output/pointcloudDM/')
    if dataGas.shape[0] > 0: generate_pointcloud_xyzw(dataGas, 'output/pointcloudGas/')            

    if True:
        dataDMVel = np.loadtxt('output/dataDM_velocity.np', float, delimiter=',')
        dataGasVel = np.loadtxt('output/dataGas_velocity.np', float, delimiter=',')
        
        if dataDMVel.shape[0] > 0: generate_fits(dataDMVel, 'output/dataDM_velocity.fits', 3)
        if dataGasVel.shape[0] > 0: generate_fits(dataGasVel, 'output/dataGas_velocity.fits', 3)

if __name__ == '__main__':
    main()