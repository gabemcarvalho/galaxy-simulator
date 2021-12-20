from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import math
import matplotlib.pyplot as plt
from statistics import mean

#fits_image_filename = fits.util.get_testdata_filepath()
def initializePosition():
    hdulist = fits.open('474final\Data\DarkMatter_position_final.fits')
    hdulist.info()

    Positionarray = np.array(hdulist[0].data)
    hdulist.close()
    return Positionarray

def initializeVelocity():
    vdulist = fits.open('474final\Data\DarkMatter_velocity.fits')
    vdulist.info()

    Velocityarray = np.array(vdulist[0].data)
    vdulist.close()
    return Velocityarray  

Velocityarray = initializeVelocity()
Positionarray = initializePosition()

xmean = mean(Positionarray[:,0].tolist())
ymean = mean(Positionarray[:,1].tolist())
zmean = mean(Positionarray[:,2].tolist())
smoothing = Positionarray[:,3].tolist()

rlist = []
circvlist = []
for i in range(len(Positionarray[:,0])):
    x = Positionarray[i,0] - xmean
    y = Positionarray[i,1] - ymean
    z = Positionarray[i,2] - zmean

    xv = Velocityarray[i,0]
    yv = Velocityarray[i,1]
    zv = Velocityarray[i,2]

    r = math.sqrt(x**2+y**2+z**2)
    circv = (xv*y - x*yv)/(x**2+y**2)
    rlist.append(r)
    circvlist.append(abs(circv))


plt.scatter(rlist, circvlist, s=2, c= 'k')
plt.title("Rotation Curve of Dark Matter")
plt.xlabel("R in kpc")
plt.xlim(0, 3)
plt.ylim(0,20)
plt.ylabel("Circular Velocity (kpc/Gyr)")
plt.savefig("DarkMatterRotationCurve.png")