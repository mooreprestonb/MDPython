from __future__ import print_function

import numpy as np
from scipy.integrate import simps
from numpy import trapz

file = open("well.dat","r")

lines = file.readlines()
data = np.zeros((len(lines)))
print("")

for i in range(len(lines)):
    data[i] = float(lines[i].split()[1])
     
print("Data: ",data)
print("Length of array:",len(data))
# composite trapezoidal rule.
area = trapz(data, dx=0.0993957)

dt = data[1]-data[0]
print("dt: ",dt)

xmin = float(lines[0].split()[0])
print("xmin: ",xmin)
xmax = float(lines[-1].split()[0])
print("xmax: ",xmax)


mom1 = (1/(xmax-xmin)*area)
print("First Moment: ", mom1)

mom2 = np.std(data)
print("Second Moment: ", mom2)
