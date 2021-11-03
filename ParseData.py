# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:02:26 2021

@author: Gabriel Carvalho
"""

import os
import numpy as np
from astropy.io import fits
from pathlib import Path

fits_filename = 'output/data.fits'
xyz_directory = 'output/pointcloud/'

def main():
    
    generate_fits = True
    generate_xyz = True
    
    data = np.loadtxt('output/data.np', float, delimiter=',')
    
    if generate_fits:
        # Convert array to correct format
        positions = np.empty([3, data.shape[0], int(data.shape[1] / 3)])
        for i in range(int(data.shape[1] / 3)):
            positions[0,:,i] = data[:,3*i]
            positions[1,:,i] = data[:,3*i+1]
            positions[2,:,i] = data[:,3*i+2]
        
        # Remove old file
        try:
            os.remove(fits_filename)
        except OSError:
            pass
        
        if os.path.exists(fits_filename):
            print('Failed to remove old .fits file')
            return
        
        # Make a new FITS file
        hdu = fits.PrimaryHDU(positions)
        hdu.writeto(fits_filename)
    
    if generate_xyz:
        Path(xyz_directory).mkdir(parents=True, exist_ok=True)
        
        for i in range(data.shape[0]):
            filename = xyz_directory + '{:04d}'.format(i) + '.xyz'
            with open(filename, 'w') as f:
                for j in range(int(data.shape[1] / 3)):
                    # write x,y,z as x,z,y (since y is up in TouchDesigner)
                    f.write(str(data[i,3*j]) + "," + str(data[i,3*j+2]) + "," + str(data[i,3*j+1]) + '\n')
                    
            

if __name__ == '__main__':
    main()