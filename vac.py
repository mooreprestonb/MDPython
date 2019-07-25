#!/usr/bin/python3

import numpy
import sys
import time
import math
from scipy.fftpack import dct

file = open(sys.argv[1],"r")

lines = file.readlines()

nwindow = 2**10

print("Called with:",sys.argv)
print("First line of",sys.argv[0],lines[0])
print("Window length:",nwindow)

natoms = lines[0].split()[1]
natoms = int(natoms)
dt = float(lines[natoms+1].split()[3]) - float(lines[0].split()[3])
print(natoms,nwindow,dt)

nconf = len(lines)//(natoms+1)
print("Number of configurations: ",nconf)

data = numpy.zeros((nwindow,natoms,3))
vac = numpy.zeros(nwindow)

# print(lines[0:10])
for i in range(nwindow): # get initial data
    for j in range(natoms):
        line  = lines[i*(natoms+1)+j+1].split()
        #print(i,j,line)
        data[i][j] = line

itime = 10
tnow = time.time()
ttime = tnow

print("Correlating")
for i in range(nconf-nwindow):
    # Correlate
    #vac += numpy.einsum('ijk,jk',data,data[0])
    for k in range(nwindow):
        for j in range(natoms):
            vac[k] += numpy.dot(data[k][j],data[0][j])
            
    # get next data window
    for j in range(nwindow-1): # move array down
        data[j] = data[j+1]
    for j in range(natoms): # read data into end of array
        data[-1][j] = lines[i*(natoms+1)+j+1].split()

    if(itime < time.time()-tnow): # report where we are
        print('conf = {}/{} = {:.4f}%, time = {:g}'.format(i,nconf,i/nconf*100,time.time()-ttime))
        tnow = time.time()

vac /= (nconf-nwindow)*natoms # normalize

v2 = vac[0]
vac /= vac[0] # set initial value to 1

vac_dct = dct(vac)
#print (vac_dct)

print("Writing",sys.argv[2])
file = open(sys.argv[2],"w")
file.write("# vac of velocity data, v2 = " + str(v2) +"\n") 
for i in range(nwindow):
    line = str(dt*i) + " " + str(vac[i]) + "\n"
    file.write(line)

file.write("\n# cosine transform\n")
for i in range(nwindow):
    line = str(i*math.pi/(dt*nwindow))+ " " + str(vac_dct[i]) + "\n"
    file.write(line)

