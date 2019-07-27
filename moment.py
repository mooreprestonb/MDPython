#takes a file of two columns, x,y sorted by x without comments at the top

#need a way to specify which x-values to search
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

threshold = 0.20
samples = []
for dF in data:
    if dF and dF > threshold:
        samples.append(dF)

print("Data: ",samples)
print("Length of array:",len(data))

dt = data[1]-data[0]
print("dt: ",dt)
xmin = float(lines[0].split()[0])
print("xmin: ",xmin)
xmax = float(lines[-1].split()[0])
print("xmax: ",xmax)

# composite trapezoidal rule.
area = trapz(samples, dx=dt)
mom1 = (1/(xmax-xmin)*area)
print("First Moment: ", mom1)

mom2 = np.std(data)
print("Second Moment: ", mom2)

#sample data
#1.38149 0.0869713
#1.40832 0.186367
#1.43515 0.173943
#1.46198 0.198791
