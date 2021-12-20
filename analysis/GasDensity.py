from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sp
import os


def initializePosition(file,extension):
    if extension == '.fits':
        hdulist = fits.open(file)
        print("Reading fits file:", file)
        Positionarray = np.array(hdulist[0].data)
        hdulist.close()
    else:
        print("Reading file:", file)
        Positionarray = np.genfromtxt(file, delimiter=',')

    return Positionarray


def gausskernel(r,h):
    return (h * math.sqrt(math.pi))**-3 * math.exp(-r * r / (h * h))

def kernel(r, h):
    q = r / h
    if (q < 0.5):
        return 8.0 /math.pi/ (h*h*h) * (1.0 - 6.0 * q * q + 6.0 * q**3)   
    elif (q < 1.0):
        return 16.0 / math.pi / (h*h*h) * (1.0 - q)**3
    return 0.0

def toykernel(r, h):
    q = r / h
    Cs = 1/(4*math.pi*(h*h*h))
    if (q < 1):
        return Cs * ((2 - q)**3 - (4*(1 - q)**3))   
    elif (q < 2):
        return Cs * ((2 - q)**3)
    return 0.0


def density(x0, y0, z0, h, M, positions):

  rho = 0
  N = len(positions[:,0])
  m = M/N

  for i in range(N):
    x = positions[i,0]
    y = positions[i,1]
    z = positions[i,2]
    r = math.sqrt((x0-x)**2+(y0-y)**2+(z0-z)**2)
    rho += m * gausskernel(r, h)
  return rho

m = 2/1000
h = 0.125
y = 0
z = 0
M = 2
R = 0.75
n = 1

file = '474final\Data\Gas_position_final.fits'
name, extension = os.path.splitext(file)
Positionarray = initializePosition(file,extension)

rlist = []
dlist = []
validlist = []
jlist = []

for i in np.arange(0,1,0.01):
    dx = density(i,0,0, h, M, Positionarray)
    dy = density(0,i,0, h, M, Positionarray)
    dz = density(0,0,i, h, M, Positionarray)

    valid = (M*math.pi**(-3/2)/(R**3))*(sp.special.gamma(5/2+n)/sp.special.gamma(1+n))*(1-(i**2/(R**2)))**n
    if valid >= 0:
        jlist.append(i)
        validlist.append(valid)
    rlist.append(i)
    dlist.append((dx+dy+dz)/3)

plt.scatter(rlist, dlist, s=2, c= 'k')
plt.plot(jlist,validlist, c = 'r')
plt.title("Density")
plt.ylim(0,3)
plt.xlabel("r")
plt.ylabel("Density")
plt.savefig("GasDensity.png")