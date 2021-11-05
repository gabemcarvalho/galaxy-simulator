# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:02:26 2021

@author: Gabriel Carvalho
"""

import os
import numpy as np
from astropy.io import fits
from pathlib import Path

def generate_fits(data, filename):
    # Convert array to correct format
    positions = np.empty([3, data.shape[0], int(data.shape[1] / 3)])
    for i in range(int(data.shape[1] / 3)):
        positions[0,:,i] = data[:,3*i]
        positions[1,:,i] = data[:,3*i+1]
        positions[2,:,i] = data[:,3*i+2]
    
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
    
def generate_pointcloud(data, directory):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    for i in range(data.shape[0]):
        filename = directory + '{:04d}'.format(i) + '.xyz'
        with open(filename, 'w') as f:
            for j in range(int(data.shape[1] / 3)):
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,3*j]) + "," + str(data[i,3*j+2]) + "," + str(data[i,3*j+1]) + '\n')
                
    print("Exported pointcloud to" + directory)

def main():
    
    dataDM = np.loadtxt('output/dataDM.np', float, delimiter=',')
    dataGas = np.loadtxt('output/dataGas.np', float, delimiter=',')
    
    generate_fits(dataDM, 'output/dataDM.fits')
    generate_fits(dataGas, 'output/dataGas.fits')
    
    generate_pointcloud(dataDM, 'output/pointcloudDM/')
    generate_pointcloud(dataGas, 'output/pointcloudGas/')            

if __name__ == '__main__':
    main()