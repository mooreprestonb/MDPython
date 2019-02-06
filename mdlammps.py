#!/usr/bin/python3
# run MD using lammps style input

import sys
import numpy

import mdinput   # file with input routines
import mdoutput  # file with output routines
import mdbond    # file with bonding routines
import mdlj      # file with non-bonded routines

# assigned global variables
global mass         # masses of each type
global pos          # positions of each atom
global vel          # velocities of each atom
global acc          # acceleration of each atom
global box          # box coordinates
global natoms       # number of atoms
global atypes       # number of atom types
global aatype       # array of atom types    
global nbonds       # number of bonds
global tbonds       # number of bond types
global bonds        # bonds
global abtype       # array of bond types    
global masses       # broadcast of mass type to atoms

box = numpy.zeros(3)  # initialize box array
    
#-------------------------------------------------
def readinit(datafile): # read lammps init data file
    
    global natoms, atypes, nbonds, tbonds, box
    
    print ("Reading",datafile)
    fi = open(datafile,"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()

    data =[0]*7
    mdinput.readinvals(lines,data)

    # unpack or destructuring
    natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = data
    print("Natoms",natoms," Atypes",atypes," Bonds",nbonds," Btypes",tbonds)
    print("Box",box)

    # allocate arrays from data
    global mass, aatype, pos, vel, acc, masses, bonds
    
    mass = numpy.zeros(atypes)
    aatype = numpy.zeros(natoms,dtype=int)
    pos = numpy.zeros((natoms,3))
    vel = numpy.zeros((natoms,3))
    acc = numpy.zeros((natoms,3))
    masses = numpy.zeros((natoms,3))
    bonds =numpy.zeros((nbonds,3),dtype=int)

    mdinput.getmasses(lines,atypes,mass)
    mdinput.getatoms(lines,natoms,aatype,pos,mass,masses)
    mdinput.getvel(lines,natoms,vel)
    mdinput.getbonds(lines,nbonds,bonds)
    
#-------------------------------------------
def readin():

    global nsteps, dt, initfile, ithermo, idump, dumpfile, bond_style, bondcoeff

    print ("Reading",sys.argv[1])
    fi = open(sys.argv[1],"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()
    
    # print lines
    data =[0]*7
    mdinput.readsysvals(lines,data) # read lammps in file
    nsteps, dt, initfile, ithermo, idump, dumpfile, bond_styles = data
    print("nsteps, dt, initfile, ithermo, idump, dumpfile, bond_styles",data)
    
    readinit(initfile)  # read initfile to get atoms bonds and types

    global bondcoeff 
    bondcoeff = numpy.zeros((tbonds,3))  # read in kbond r0
    bond_style = mdinput.getbondcoeff(lines,bond_styles,tbonds,bondcoeff)
        
#-----------------------------------------------------------
def force():
    global masses
    global pbond, nbonds, bonds, bondcoeff
    global pos, vel, acc

    acc.fill(0) # zero out forces/acceration
    # lj
    # bonds
    pbond = 0
    if bond_style == 0: # harmonic bonds
        pbond = mdbond.bond_harm(nbonds,bonds,bondcoeff,pos,acc)
    elif bond_style == 1: # morse bonds
        pbond = mdbond.bond_morse(nbonds,bonds,bondcoeff,pos,acc)
    else:
        print ("Error bond style not found!")
        exit(1)
    # bend
    # torsion
    acc /= masses  # change forces into accelerations

#-----------------------------------------------------------    
def step():
    global pos, vel, acc, dt
    # print istep,pos,vel,acc
    vel += acc*dt/2.0
    pos += vel*dt
    force()
    vel += acc*dt/2.0

#-----------------------------------------------------------    

# read command line for input file
if (len(sys.argv) == 2):  # error check that we have an input file
    infile = sys.argv[1]
else:
    print("No input file?")
    exit(1)

print (sys.argv)
readin()
print (nsteps, dt)

# inital force and adjustments
mdlj.zero_momentum(masses,vel)  # zero the momentum
force()

print("Running dynamics")
for istep in range(nsteps):

    if(istep%ithermo==0):
        ke = .5*numpy.sum(masses*vel*vel)
        print (istep,ke+pbond,ke,pbond)

    if(istep%idump==0):
        mdoutput.write_dump(dumpfile,istep,natoms,pos,aatype)

    step()
        
print("Done dynamics!")

mdoutput.write_init("test.init",istep,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype)
exit(0)
