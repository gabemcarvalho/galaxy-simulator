import os, math
import numpy as np
from pathlib import Path

def gausskernel(r,h):
    return (h * math.sqrt(math.pi))**-3 * math.exp(-r * r / (h * h))

def generate_pointcloud_with_density(data, directory, clear_folder, start_index=0):
    Path(directory).mkdir(parents=True, exist_ok=True)
    
    if clear_folder == True:
        for f in os.listdir(directory):
            os.remove(os.path.join(directory, f))
    
    for i in range(data.shape[0]):
        filename = directory + '{:06d}'.format(start_index + i) + '.xyz'
        with open(filename, 'w') as f:
            N = int(data.shape[1] / 4)
            
            for j in range(N):
                d = density(j, h=0.125, M=2, N=N, step=i, data=data)
                
                # write x,y,z as x,z,y (since y is up in TouchDesigner)
                f.write(str(data[i,4*j]) + "," + str(data[i,4*j+2]) + "," + str(data[i,4*j+1]) + "," + str(d) + '\n')

def density(particle, h, M, N, step, data):
    rho = 0
    m = M/N
    
    x0 = data[step,4*particle+0]
    y0 = data[step,4*particle+1]
    z0 = data[step,4*particle+2]
    
    for i in range(N):
        x = data[step,4*i+0]
        y = data[step,4*i+1]
        z = data[step,4*i+2]
        r = math.sqrt((x0-x)**2+(y0-y)**2+(z0-z)**2)
        rho += m * gausskernel(r, h)
    return rho

data = np.loadtxt('../output/old_output/13_toy_star_final_400_parts/Gas_position.np', float, delimiter=',')
generate_pointcloud_with_density(data, '../output/old_output/13_toy_star_final_400_parts/pointcloudDensity/', True, 0)
