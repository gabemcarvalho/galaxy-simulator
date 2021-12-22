import os, math
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../output/old_output/13_toy_star_final_400_parts/pointcloudDensity/000999.xyz', float, delimiter=',')
x = data[:,0]
y = data[:,1]
z = data[:,2]
d = data[:,3]

fig = plt.figure(dpi=300, figsize=(20,6))
ax = fig.add_subplot(111, projection='3d')

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

p = ax.scatter(x, y, z, c=d, cmap=plt.cm.copper, s=30)

fig.colorbar(p, ax=ax, shrink=0.7, label='Density')
plt.show()