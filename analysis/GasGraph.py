from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import math
import matplotlib.pyplot as plt
from statistics import mean

#fits_image_filename = fits.util.get_testdata_filepath()
def initializePosition():
    hdulist = fits.open('dataGas.fits')
    hdulist.info()

    Positionarray = np.array(hdulist[0].data)
    FinalPositionarray = Positionarray[:,-1,:]
    hdulist.close()
    return FinalPositionarray

def initializeVelocity():
    vdulist = fits.open('dataGas_velocity.fits')
    vdulist.info()

    Velocityarray = np.array(vdulist[0].data)
    FinalVelocityarray = Velocityarray[:,-1,:]
    vdulist.close()
    return FinalVelocityarray  

Positionarray = initializePosition()
Velocityarray = initializeVelocity()

xmean = mean(Positionarray[0].tolist())
ymean = mean(Positionarray[1].tolist())
zmean = mean(Positionarray[2].tolist())

rlist = []
circvlist = []
for i in range(len(Positionarray[0])):
    x = Positionarray[0,i] - xmean
    y = Positionarray[1,i] - ymean
    z = Positionarray[2,i] - zmean

    xv = Velocityarray[0,i]
    yv = Velocityarray[1,i]
    zv = Velocityarray[2,i]

    r = math.sqrt(x**2+y**2+z**2)
    circv = (xv*y - x*yv)/(x**2+y**2)
    rlist.append(r)
    circvlist.append(abs(circv))

print(rlist)
print(circvlist)

plt.scatter(rlist, circvlist, s=2, c= 'k')
plt.title("Rotation Curve of Galactic Gas")
plt.xlabel("R in kpc")
plt.xlim(0, 20)
plt.ylim(0,250)
plt.ylabel("Circular Velocity (kpc/Gyr)")
plt.savefig("GasRotationCurve.png")