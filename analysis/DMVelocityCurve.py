from re import A
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import math
import matplotlib.pyplot as plt
from statistics import mean
from scipy.signal import savgol_filter
import os

def main(file1,file2,file3):
    Positionarray = initializePosition(file2)
    XVel, YVel, ZVel = initializeVelocity(file1,file2,file3)
    avgr, avgy = equation(Positionarray, XVel, YVel, ZVel)
    return avgr,avgy

def initializePosition(file):
    Positionarray = np.genfromtxt(file, delimiter=',')
    return Positionarray

def initializeVelocity(file1,file2,file3):
    Positionarray = np.genfromtxt(file1, delimiter=',')
    Positionarray3 = np.genfromtxt(file3, delimiter=',')
    XVel = (Positionarray3[:,0] - Positionarray[:,0])/(2*0.04)
    YVel = (Positionarray3[:,2] - Positionarray[:,2])/(2*0.04)
    ZVel = (Positionarray3[:,1] - Positionarray[:,1])/(2*0.04)
    return XVel, YVel, ZVel 

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

def rollingaverage(average):
    y = savgol_filter(average,19,3)
    return  y

def equation(Positionarray, XVel, YVel, ZVel):
    rlist = []
    circvlist = []

    for i in range(len(Positionarray[:,0])):

        x = Positionarray[i,0]
        y = Positionarray[i,2]
        z = Positionarray[i,1]

        xv = XVel[i]
        yv = YVel[i]
        zv = ZVel[i]

        r = math.sqrt(x**2+y**2+z**2)
        circv = (xv*y - x*yv)/(x**2+y**2)*r
        rlist.append(r)
        circvlist.append(abs(circv))

    rotationarray = np.array(list(zip(rlist, circvlist)))
    sortedarray = rotationarray[rotationarray[:, 0].argsort()]
    avgr, avgy = averaging(sortedarray)
    avgy = rollingaverage(avgy)

    return avgr, avgy

file1 = '474final\Data\position1.xyz'
file2 = '474final\Data\position2.xyz'
file3 = '474final\Data\position3.xyz'
avgr1, avgy1 = main(file1,file2,file3)
plt.plot(avgr1, avgy1, c= 'k',alpha= 0.6)

file1 = '474final\Data\position4.xyz'
file2 = '474final\Data\position5.xyz'
file3 = '474final\Data\position6.xyz'
Positionarray = initializePosition(file2)
avgr2, avgy2 = main(file1,file2,file3)

plt.plot(avgr2, avgy2,'k')


V0 = 0.7
c = 30
r0 = 0.05

jlist = []
klist = []
halolist = []
cmlist = []

for i in np.arange(0,30,0.1):
    halo = i*V0/(math.sqrt(i**2+c**2))
    CentMass = math.sqrt(0.014*i**2/(i**2+r0**2)**(3/2))
    if halo >= 0:
        jlist.append(i)
        halolist.append(halo)
    if CentMass >= 0:
        klist.append(i)
        cmlist.append(CentMass)

plt.plot(jlist,halolist,"k--")
plt.plot(klist,cmlist, 'k:')
plt.xlabel("Radius")
plt.xlim(0,30)
plt.ylim(0,0.8)
plt.ylabel("Circular Velocity")
plt.savefig("DarkMatterRotationCurve.png")
