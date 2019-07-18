#!/usr/bin/python3

import numpy
import sys
import time
from scipy.fftpack import dct

file = open("vel.dat","r")

lines = file.readlines()

nwindow = 2**10
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

# vac /= vac[0] # set initial value to 1

vac_dct = dct(vac)
print (vac_dct)

print("Writing vac.dat")
file = open("vac.dat","w")
for i in range(nwindow):
    line = str(dt*i) + " " + str(vac[i]) + " " + str(vac_dct[i]) + "\n"
    file.write(line)

#vac_fft = numpy.fft.fft(vac)
#Issues:
    #1.) correlation not going negative
    #2.) FT
