#takes a file of two columns, x,y sorted by x without comments at the top
from __future__ import print_function
import math
import numpy as np
from scipy.integrate import simps
from numpy import trapz

file = open("harmonic.eig.dat","r") # Figure out how to use wildcard

lines = file.readlines()
print("")

data = []
domain = []
rangemin = 1.
rangemax = 2.5
for i in lines: #Grab the data whose y-values are relevant, make float
    if rangemax > float(i.split()[0]) > rangemin:
        domain.append(float(i.split()[0]))
        data.append(float(i.split()[1]))

x = np.array(domain)
sig0 = np.array(data)   # f_x
sig1 = sig0*x   # x*f_x
sig2 = sig1*x   # x^2*f_x

print("x-values: ",x)
print("y-values: ",sig0)

dt = domain[2]-domain[1] # Calculate integration step; changes with bin count
print("dx: ",dt)
xmin = domain[0] #First one
print("xmin: ",xmin)
xmax = domain[-1] # Last one
print("xmax: ",xmax)

# composite trapezoidal rule.
area = trapz(sig1, dx=dt)
integral = trapz(sig0, dx=dt)
extra_int = trapz(sig2, dx=dt)

print("Distribution Integral: ", integral)

mom1 = ((1/(integral))*area)
print("First Moment: ", mom1)

mom2 = 2*math.sqrt(((1/(integral))*extra_int) - mom1**2) # Variance
print("Second Moment: ", mom2)

file = open("moment.dat","w")
#line = str("XXX")+" "+str(mom1)+" "+str(mom2)+ "\n"
line = str(mom1)
file.write(line)

#sample data
#1.38149 0.0869713
#1.40832 0.186367
#1.43515 0.173943
#1.46198 0.198791
