#!/usr/bin/python

import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math

#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Autocorrelate data (time data1 data2 ...)')
parser.add_argument('input', help="input file")
parser.add_argument('output', nargs='?', help="output file")
parser.add_argument('-n',type=int,default=-1,help="#'s to correllate")
parser.add_argument('-ts',type=float,help="set time step i.e dt")
parser.add_argument('-rp',type=int,default=10,help="Report progress every \"rp\" steps")
parser.add_argument("-av",help="do NOT subtract average",action="store_true")
parser.add_argument("-s",help="Scale C[0] to 1",action="store_true")
parser.add_argument("-gr",help="Graph results",action="store_true")
parser.add_argument("-td",help="Average Correlation functions",action="store_true")
args = parser.parse_args()

print args

data = np.loadtxt(args.input) # read in special numpy

#with open(args.input) as f:
#    lines = f.readlines()

#print len(lines),"lines read in, parsing"

time = 0
tol = 1e-8

x = data[:,0] # First column
dy = np.array(data[:,1:]) # all other columns
time = x.shape[0]

print "time data",dy.shape

dt = x[1]-x[0]
for i in range(1,time):
    if ( abs(((x[i]-x[i-1])-dt)/dt) >tol):
        print "Warning steps not equal at step i",i,x[i],x[i-1],dt
        
if (args.n == -1): # set numbers to correlate
    mnum = time
else:
    mnum = min(args.n,time)

if(args.ts): # set new time step if needed
    dt = args.ts
    
corr = np.zeros((mnum,dy.shape[1])) # initalize array
xt = np.arange(mnum)*dt # initalize array

avg = np.zeros((dy.shape[1]))
avg = np.average(dy,axis=0) # average down columns?

print "Averages: ",avg
print "Variance: ",np.var(dy)

if(not args.av):
    print "Subtracting off average"
    dy = dy - avg # subtract off average

outp = max(time/args.rp,1)
hdr = " ".join(sys.argv) # create header from sys.argv

print "Correlating: Reporting progress every",outp,"steps, with a window of",mnum,"steps and a total time of ",time

for i in range(time):
    if (not (i%outp)):
        print i,"/",time

    # slice out what we want i to i+mnum (time window)
    for j in range(dy.shape[1]):
        y = dy[i:min(time,i+mnum),j]
        y0 = y[0] # slice out step i
        corrd = np.dot(y,y0)  # correlate step i using numpy dot product with i .. i+mnum
        corr[:corrd.shape[0],j] += corrd # add to corr column j with correct shape

#print "Corr:",corr

print "Done Correlating: writing out"
    
for i in range(mnum): # average 
    corr[i] = corr[i]/(time-i) #(c[0] has ntime data, c[1] has ntime-1 data etc)

if(args.td):
    print "Summing"
    td = np.zeros(mnum) # initalize array
    td = np.sum(corr,axis=1) #
    
# print corr,corrtd
print ", Now C[0] = ",corr[0]

if(args.s): # scale to 1
    print "Scaling corr[0] to 1, C[0] = ",corr[0]
    corr = corr/corr[0]
    if(args.td):
        td = td/td[0] # scale to 1

if(args.td):
    outa = np.column_stack((xt,corr,td))  # need to pass one argument to get x y
else:
    outa = np.column_stack((xt,corr))  # need to pass one argument to get x y

if(args.output==None):
    output = args.input + ".corr"
else:
    output = args.output

np.savetxt(output,outa,fmt="%f",header=hdr) # output with x y z format
    
if(args.gr): # use graphs
    print "Creating graphs"
    plt.plot(xt,corr,label="Correlated data")
    if(args.td):
        plot(xt,corr,label="Total")
    plt.legend()
    plt.show()

print "Done"
