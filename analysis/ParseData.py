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

def generate_pointcloud_xyz(data, directory, clear_folder, start_index=0):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    if clear_folder == True:
        for f in os.listdir(directory):
            os.remove(os.path.join(directory, f))
    
    for i in range(data.shape[0]):
        filename = directory + '{:06d}'.format(start_index + i) + '.xyz'
        with open(filename, 'w') as f:
            for j in range(int(data.shape[1] / 3)):
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,3*j]) + "," + str(data[i,3*j+2]) + "," + str(data[i,3*j+1]) + '\n')
                
    print("Exported pointcloud to " + directory)
    
def generate_pointcloud_xyzw(data, directory, clear_folder, start_index=0):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    if clear_folder == True:
        for f in os.listdir(directory):
            os.remove(os.path.join(directory, f))
    
    for i in range(data.shape[0]):
        filename = directory + '{:06d}'.format(start_index + i) + '.xyz'
        with open(filename, 'w') as f:
            for j in range(int(data.shape[1] / 4)):
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,4*j]) + "," + str(data[i,4*j+2]) + "," + str(data[i,4*j+1]) + "," + str(data[i,4*j+3]) + '\n')
                
    print("Exported pointcloud to " + directory)

def data_exists(filename):
    return os.path.isfile(filename) and os.stat(filename).st_size > 0

def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def generate_position_file_parts(filename, num_parts, fits_directory, xyz_directory, clear_folder=True):
    Path(fits_directory).mkdir(parents=True, exist_ok=True)
    Path(xyz_directory).mkdir(parents=True, exist_ok=True)
    
    if clear_folder == True:
        for f in os.listdir(fits_directory):
            os.remove(os.path.join(fits_directory, f))
            
    if clear_folder == True:
        for f in os.listdir(xyz_directory):
            os.remove(os.path.join(xyz_directory, f))
    
    num_steps = file_len(filename)
    steps_per_part = int(num_steps / num_parts)
    generated_steps = 0
    if num_parts > 1:
        for i in range(num_parts - 1):
            data = np.genfromtxt(filename, float, delimiter=',', skip_header=generated_steps, max_rows=steps_per_part)
            generate_fits(data, fits_directory + '{:06d}'.format(generated_steps) + '.fits', 4)
            generate_pointcloud_xyzw(data, xyz_directory, False, generated_steps)
            generated_steps += steps_per_part
            print('[' + str(i + 1) + "/" + str(num_parts) + ']')
        
    data = np.genfromtxt(filename, float, delimiter=',', skip_header=generated_steps) # read any remaining steps
    generate_fits(data, fits_directory + '{:06d}'.format(generated_steps) + '.fits', 4)
    generate_pointcloud_xyzw(data, xyz_directory, False, generated_steps)
    print('[' + str(num_parts) + "/" + str(num_parts) + ']')
    print('Finished parsing positions')

def main():
    
    output_dir = '../output/'
    
    # position
    filename = 'DarkMatter_position'
    if data_exists(output_dir + filename + '.np'):
        generate_position_file_parts(output_dir + filename + '.np', 35, output_dir + 'positions_fits/', output_dir + 'pointcloudDM/', True)
    
    filename = 'Gas_position'
    if data_exists(output_dir + filename + '.np'):
        generate_position_file_parts(output_dir + filename + '.np', 40, output_dir + 'positions_fits/', output_dir + 'pointcloudGas/', True)
    
    # final position
    filename = 'DarkMatter_position_final'
    if data_exists(output_dir + filename + '.np'):
        data = np.loadtxt(output_dir + filename + '.np', float, delimiter=',')
        generate_fits_singlestep(data, output_dir + filename + '.fits', 4)
    
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