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
  
def generate_fits_singlestep(data, filename, dim):
    # Convert array to correct format
    positions = np.empty([int(data.shape[0] / dim), dim])
    for i in range(int(data.shape[0] / dim)):
        positions[i,:] = data[i*dim:i*dim+dim]
    
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

def data_exists(filename):
    return os.path.isfile(filename) and os.stat(filename).st_size > 0

def main():
    
    output_dir = '../output/'
    
    # position
    filename = 'DarkMatter_position'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits(data, output_dir + filename + '.fits', 3)
        generate_pointcloud_xyz(data, output_dir + 'pointcloudDM/')
    
    filename = 'Gas_position'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits(data, output_dir + filename + '.fits', 4)
        generate_pointcloud_xyzw(data, output_dir + 'pointcloudGas/')
    
    # final position
    filename = 'DarkMatter_position_final'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 3)
    
    filename = 'Gas_position_final'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 4)
    
    # velocity
    filename = 'DarkMatter_velocity'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 3)
    
    filename = 'Gas_velocity'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 3)
    
    # density
    filename = 'DarkMatter_density'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 1)
    
    filename = 'Gas_density'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 1)
        
    print('Done!')

if __name__ == '__main__':
    main()