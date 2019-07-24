#!/usr/bin/python3
# run MD using lammps style input

import sys
import numpy
import time
import re
import math

#import mdglobal  # file with global variables
import mdinput   # file with input routines
import mdoutput  # file with output routines
#import mdlj      # file with non-bonded routines
import mdbond    # file with bonding routines
#import mdbend    # file with bending routines
import mdstep    # file with integration routines

# assigned global variables
global natoms       # number of atoms
global atypes       # number of atom types
global nbonds       # number of bonds
global tbonds       # number of bond types
global box          # box coordinates
global pot          # potential components
global mass         # masses of each type
global masses       # broadcast of mass type to atoms
global pos          # positions of each atom
global vel          # velocities of each atom
global acc          # acceleration of each atom
global aatype       # array of atom types
global bonds        # bonds (array with type, ibond, jbond)
global hessian      # hessian matrix
global abtype       # array of bond types
global logfile      # file to output thermodata
global ensamble

box = numpy.zeros(3)
pot = numpy.zeros(6)

#-------------------------------------------
def readin(): # read lammps like infile

    global nsteps, dt, initfile, ithermo, idump, dumpfile, bond_style, bondcoeff
    global logfile, inmfile, inmo
    global bondcoeff, reps, ensamble, fix_type
    global natoms, atypes, nbonds, tbonds, box

    # print lines
    data, bond_style, bondcoeff, reps, fix_type, var_lst = mdinput.readsysvals(sys.argv[1]) # read lammps in file
    dt, initfile, bond_styles, idump, dumpfile, ithermo, logfile, inmfile, inmo, nsteps = data
    print("dt, initfile, bond_styles, idump, dumpfile, ithermo, logfile, inmfile, inmo, nsteps",data)
    
    natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = mdinput.readinvals(initfile) 
    print("Natoms",natoms," Atypes",atypes," Bonds",nbonds," Btypes",tbonds)
    print("Box",box)

    # allocate arrays from data
    global mass, aatype, pos, vel, acc, masses, bonds, hessian, zeta, Q

    acc = numpy.zeros((natoms,3))

    mass, aatype, pos, vel, masses, bonds = mdinput.make_arrays(initfile,reps)

    if re.search('nve',fix_type,flags=re.IGNORECASE):
        ensamble = 0
    elif re.search('nvt',fix_type,flags=re.IGNORECASE):
        global T, Tdamp, Q
        var = [float(num) for num in var_lst[:3]]
        T, T, Tdamp = var
        Q = numpy.array([3*natoms*T*Tdamp*Tdamp]*2)
        ensamble = 1
    else:
        print("Error: no ensamble? ",fix_type)
        exit(1)

#-----------------------------------------------------------

# read command line for input file
if (len(sys.argv) != 2):  # error check that we have an input file
    print("No input file? or wrong number of arguments")
    exit(1)

if not re.search(re.compile(r'.+\.in'),sys.argv[1]):
    print('Incorrect input file type.')
    exit(1)

print (sys.argv)

readin() # read infile

# inital force and adjustments
mdstep.zero_momentum(masses,vel)  # zero the momentum
mdstep.force(natoms,pot,nbonds,bond_style,bondcoeff,bonds,masses,pos,vel,acc)
teng = mdoutput.write_thermo(logfile,0,natoms,masses,pos,vel,pot)

itime = 1
tnow = time.time()
ttime = tnow
tol = 1e-8
dump_vel = 5 # how often to dump velocities
if(fix_type == 'nvt'):
    ensamble = 1 # nvt
else:
    ensamble = 0 # nve

print("Running dynamics")

eig_array = [] # empty array for the eigenvalues

for istep in range(1,nsteps+1):

    mdstep.step(natoms,ensamble,dt,pot,nbonds,bond_style,bondcoeff,bonds,masses,pos,vel,acc)

    if(istep%inmo==0): # get instantaneous normal modes
        hessian = numpy.zeros((pos.size,pos.size))
        mdbond.inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)

        # print(hessian)
        w,v = numpy.linalg.eig(hessian)
        # remove 3 lowest eigegvalues (should be translations of entire system)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        eig_array.extend(w.real) # only get real part of array - imag do to round off error is small so we throw away.
        # mdoutput.write_inm(istep,hessian)

    if(istep%ithermo==0): # write out thermodynamic data
        teng = mdoutput.write_thermo(logfile,istep,natoms,masses,pos,vel,pot)

    if(istep%idump==0): # dump to xyz file so we can see this in lammps
        mdoutput.write_dump(dumpfile,istep,natoms,pos,aatype)

    if(istep%dump_vel==0):
        mdoutput.write_dump_vel("vel.dat",istep*dt,natoms,vel)

    if(itime < time.time()-tnow): # report where we are
        print('step = {}/{} = {:.4f}%, teng = {:g}, time = {:g}'.format(istep,nsteps,istep/nsteps*100,teng,time.time()-ttime))
        tnow = time.time()

print('Done dynamics! total time = {:g} seconds'.format(time.time()-ttime))
mdoutput.write_init("test.init",istep-1,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype)

#Create histogram!
nconf = len(eig_array)
if(nconf==0):
    print("No configurations calculated eigenvalues! thus NOT calculating historgram")
else:
    print("Creating Histogram with",len(eig_array),"configurations")
    earray = numpy.array(eig_array)
    for i in range(earray.size): # sqrt eigenvaules (from w^2 to w)
        if(earray[i] < 0) :
            earray[i] = -math.sqrt(-earray[i])
        else:
            earray[i] = math.sqrt(earray[i])

    bin_ct = int(numpy.log2(nconf) + 1)
    bin_ct = 50
    histo,histedge = numpy.histogram(earray,bins=bin_ct,density=True)
    histdat = numpy.zeros((histo.size,2))
    for i in range(histo.size):
        histdat[i][0] = (histedge[i]+histedge[i+1])/2
        histdat[i][1] = histo[i]
        #print(histo,histedge,histdat)
    head = "Histogram of eigenvalues " + sys.argv[0] + " " + str(len(eig_array))
    numpy.savetxt(inmfile,(histdat),header=head,fmt="%g")

print("Done!")
exit(0)
