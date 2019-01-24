#!/usr/bin/python
# run MD using lammps style input

import sys
import numpy
import re
import math

#-------------------------------------------------
def findline1(lines,val): # find line with first word val
    for ln in range(len(lines)):  # loop over lines 
        words = lines[ln].split() # get words from line
        #print words
        if(len(words)>0): # make sure it is not a blank line
            if (words[0] == val):
                return ln
        
    print "Word",val," not found!"
    return -1

#-------------------------------------------------
def findline2(lines,val): # find line with word anywhere
    for ln in range(len(lines)):  # loop over lines
        # do we have val in line
        if re.search(val,lines[ln],flags=re.IGNORECASE): 
            return ln
            
    print "Word",val,"not found!"
    return -1

#-------------------------------------------
def readin():
    print "Reading",sys.argv[1]
    fi = open(sys.argv[1],"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()
    
    # print lines

    global nsteps
    ln = findline1(lines,"run")
    nsteps = int(lines[ln].split()[1])

    global dt
    ln = findline1(lines,"timestep")
    dt = float(lines[ln].split()[1])

    global initfile
    ln = findline1(lines,"read_data")
    initfile = (lines[ln].split()[1])

    global ithermo
    ln = findline1(lines,"thermo")
    ithermo = int(lines[ln].split()[1])

    global idump
    global dumpfile
    ln = findline1(lines,"dump")
    idump = int(lines[ln].split()[4])
    dumpfile = lines[ln].split()[5]

#-------------------------------------------------
def readinit(datafile): # read lammps init data file
    global mass         # masses of each type
    global pos          # positions of each atom
    global vel          # velocities of each atom
    global box          # box coordinates
    global natoms       # number of atoms
    global atypes       # number of atom types
    global aatype       # array of atom types    
    global nbonds       # number of bonds
    global tbonds       # number of bond types
    global bonds        # bonds
    global abtype       # array of bond types    
    global masses       # broadcast of mass type to atoms
    
    print "Reading",datafile
    fi = open(datafile,"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()

    ln = findline2(lines,"atoms")
    natoms = int(lines[ln].split()[0])

    ln = findline2(lines,"atom types")
    atypes = int(lines[ln].split()[0])

    ln = findline2(lines,"bonds")
    nbonds = int(lines[ln].split()[0])

    ln = findline2(lines,"bond types")
    tbonds = int(lines[ln].split()[0])

    box = numpy.zeros(3)
    ln = findline2(lines,"xlo xhi")
    xlo = float(lines[ln].split()[0])
    xhi = float(lines[ln].split()[1])
    box[0] = xhi-xlo

    ln = findline2(lines,"ylo yhi")
    ylo = float(lines[ln].split()[0])
    yhi = float(lines[ln].split()[1])
    box[1] = yhi-ylo

    ln = findline2(lines,"zlo zhi")
    zlo = float(lines[ln].split()[0])
    zhi = float(lines[ln].split()[1])
    box[2] = zhi-zlo

    print natoms, atypes, nbonds, tbonds

    mass = numpy.zeros(atypes)
    aatype = numpy.zeros(natoms,dtype=int)
    pos = numpy.zeros((natoms,3))
    masses = numpy.zeros((natoms,3))
    vel = numpy.zeros((natoms,3))
    bonds =numpy.zeros((nbonds,3),dtype=int)
    
    ln = findline1(lines,"Masses")
    if (len(lines[ln+1].split())==0):
        ioff = ln+2
    else:
        ioff = ln+1
        
    for i in range(atypes):
        words = lines[ioff+i].split()
        if (int(words[0]) != i+1):
            print "Error: while assigning masses at line",ln
            print "Expecting",i+1,"got",words[0],"from",lines[ln+2]
            exit(1)
        mass[i] = float(words[1])
    print "Assigned masses"

    ln = findline1(lines,"Atoms")
    if (len(lines[ln+1].split())==0):
        ioff = ln+2
    else:
        ioff = ln+1

    for i in range(natoms):
        words = lines[ioff+i].split()
        if (int(words[0]) != i+1):
            print "Error: while assigning atoms at line",ln
            print "Expecting",i+1,"got",words[0],"from",lines[ln+2]
            exit(1)
        aatype[i] = int(words[1])
        pos[i][0] = float(words[4])
        pos[i][1] = float(words[5])
        pos[i][2] = float(words[6])
        masses[i][0] = mass[aatype[i]-1] 
        masses[i][1] = mass[aatype[i]-1] 
        masses[i][2] = mass[aatype[i]-1] 
    print "Assigned atom types and positions"

    ln = findline1(lines,"Velocities")
    if (ln!=-1):
        if (len(lines[ln+1].split())==0):
            ioff = ln+2
        else:
            ioff = ln+1
        for i in range(natoms):
            words = lines[ioff+i].split()
            if (int(words[0]) != i+1):
                print "Error: while assigning velocities at line",ln
                print "Expecting",i+1,"got",words[0],"from",lines[ln+2]
                exit(1)
            vel[i][0] = float(words[1])
            vel[i][1] = float(words[2])
            vel[i][2] = float(words[3])
        print "Assigned Velocities"
    else:
        print "Velocities will start at zero"

    ln = findline1(lines,"Bonds")
    if (ln!=-1):
        if (len(lines[ln+1].split())==0):
            ioff = ln+2
        else:
            ioff = ln+1

        for i in range(nbonds):
            words = lines[ioff+i].split()
            if (int(words[0]) != i+1):
                print "Error: while assigning velocities at line",ln
                print "Expecting",i+1,"got",words[0],"from",lines[ln+2]
                exit(1)
            bonds[i][0] = int(words[1])
            bonds[i][1] = int(words[2])-1
            bonds[i][2] = int(words[3])-1
        print "Assigned Bonds"
    else:
        print "No bonds found"

#----------------------------------------------------------
def write_dump(dumpfile,istep):

    fo = open(dumpfile,"a+")
    fo.write("%d\n"% natoms)
    fo.write("# step %d\n" % istep)
    for i in range(natoms):
        line = "Ar " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        fo.write(line)
    fo.close()
    return 0
#-----------------------------------------------------------
def bond(): # harmonic
    bondk = 2
    r0 = 1.5
    global pbond
    pbond = 0
    for i in range(nbonds):
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0
        pot = bondk*dr**2
        pbond += pot
        dudr = 2.*bondk*dr
        dpos = (dudr/r)*dpos
        acc[bonds[i][1]] += dpos
        acc[bonds[i][2]] -= dpos
        
#-----------------------------------------------------------
def force():
    global acc
    global masses
    acc.fill(0)
    # lj
    # bond
    bond()
    # bend
    # torsion
    acc /= masses

#-----------------------------------------------------------    
# read command line for input file

if (len(sys.argv) == 2):  # error check that we have an input file
    infile = sys.argv[1]
else:
    print "No input file?"
    exit(1)

print sys.argv
readin()
print nsteps, dt, initfile

readinit(initfile)

global acc
acc = numpy.zeros((natoms,3))
force()

print "Running dynamics"
for istep in range(nsteps):

    if(istep%ithermo==0):
        ke = .5*numpy.sum(masses*vel*vel)
        print istep,ke+pbond,ke,pbond

    if(istep%idump==0):
        write_dump(dumpfile,istep)

    # print istep,pos,vel,acc
    vel += acc*dt/2.0
    pos += vel*dt
    force()
    vel += acc*dt/2.0
        
print "Done!"
exit(0)
