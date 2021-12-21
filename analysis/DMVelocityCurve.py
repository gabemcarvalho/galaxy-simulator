from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import math
import matplotlib.pyplot as plt
from statistics import mean

#fits_image_filename = fits.util.get_testdata_filepath()
def initializePosition():
    hdulist = fits.open('474final\Data\DarkMatter_position_final.fits')
    Positionarray = np.array(hdulist[0].data)
    hdulist.close()
    return Positionarray

def initializeVelocity():
    vdulist = fits.open('474final\Data\DarkMatter_velocity.fits')
    Velocityarray = np.array(vdulist[0].data)
    vdulist.close()
    return Velocityarray  

def averaging(circv):
    xlist = []
    average = []
    for i in np.arange(0,40,0.2):
        ylist = []
        for j in range(len(circv[:,0])):
            if circv[j,0] <= i and circv[j,0] >= i-0.2:
                ylist.append(circv[j,1])
        
        if len(ylist) > 0:
            average.append(sum(ylist)/len(ylist))
            xlist.append(i)
        ylist = []
    return xlist, average

def rollingaverage(xlist,average):
    rllavg = []
    for i in range(len(average)):
        if i < 5:
            rllavg.append(average[i])
        elif i < len(average) - 5:
            rollavg = (average[i-5]+ average[i-4] + average[i-3] + average[i-2] + average[i-1] + 
            average[i] + average[i+1] + average[i+2] + average[i+3]+ average[i+4]+ average[i+5])/11
            rllavg.append(rollavg)
        else:
            rllavg.append(average[i])

    return xlist, rllavg


Velocityarray = initializeVelocity()
Positionarray = initializePosition()

xmean = mean(Positionarray[:,0].tolist())
ymean = mean(Positionarray[:,1].tolist())
zmean = mean(Positionarray[:,2].tolist())
smoothing = Positionarray[:,3].tolist()

rlist = []
circvlist = []
dispvlist = []
for i in range(len(Positionarray[:,0])):

    x = Positionarray[i,0]
    y = Positionarray[i,1]
    z = Positionarray[i,2]

    xv = Velocityarray[i,0]
    yv = Velocityarray[i,1]
    zv = Velocityarray[i,2]

    dispersionvels = xv**2+yv**2+zv**2
    dispvlist.append(dispersionvels)
    r = math.sqrt(x**2+y**2+z**2)
    circv = (xv*y - x*yv)/(x**2+y**2)*r
    rlist.append(r)
    circvlist.append(abs(circv))

sigmadisp = sum(dispvlist)
totdispersion = math.sqrt(sigmadisp/(3*len(Positionarray[:,0])))
print("The total velocity dispersion is ", totdispersion)
rotationarray = np.array(list(zip(rlist, circvlist)))
sortedarray = rotationarray[rotationarray[:, 0].argsort()]
avgr, avgy = averaging(sortedarray)
avgr, avgy = rollingaverage(avgr, avgy)

plt.scatter(sortedarray[:,0], sortedarray[:,1], s=2, c= 'k', alpha = 0.1)
plt.plot(avgr, avgy, c= 'r')
plt.title("Rotation Curve of Dark Matter")
plt.xlabel("R in kpc")
plt.xlim(0,30)
plt.ylim(0,1.5)
plt.ylabel("Circular Velocity (kpc/Gyr)")
plt.savefig("DarkMatterRotationCurve.png")