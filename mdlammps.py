#!/usr/bin/python3
# run MD using lammps style input

import sys
import numpy
import re
import math

import mdoutput  # file with output routines

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

#-------------------------------------------------
def findline1(lines,val): # find line with first word val
    for ln in range(len(lines)):  # loop over lines 
        words = lines[ln].split() # get words from line
        #print words
        if(len(words)>0): # make sure it is not a blank line
            if (words[0] == val):
                return ln
        
    print ("Word",val," not found!")
    return -1

#-------------------------------------------------
def findline2(lines,val): # find line with word anywhere
    for ln in range(len(lines)):  # loop over lines
        # do we have val in line
        if re.search(val,lines[ln],flags=re.IGNORECASE): 
            return ln
            
    print ("Word",val,"not found!")
    return -1

#-------------------------------------------
def readin():
    print ("Reading",sys.argv[1])
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
    readinit(initfile)        # read initfile to get atoms bonds and types

    global ithermo
    ln = findline1(lines,"thermo")
    ithermo = int(lines[ln].split()[1])

    global idump
    global dumpfile
    ln = findline1(lines,"dump")
    idump = int(lines[ln].split()[4])
    dumpfile = lines[ln].split()[5]

    global bond_style
    ln = findline1(lines,"bond_style")
    bond_styles = lines[ln].split()[1]
    global bondcoeff 

    if re.search("harmonic",bond_styles,flags=re.IGNORECASE): 
        print("Reading in harmonic bonds coefficents for",tbonds,"types")
        bond_style = 0 #harmonic bond style
        bondcoeff = numpy.zeros((tbonds,2))  # read in kbond and r0
        ln = findline1(lines,"bond_coeff")
        for i in range(tbonds):
            btype = lines[ln+i].split()[1]
            if (int(btype) != i+1):
                print("Error wrong type in bond_coeff at line",ln+i)
                print(lines[ln+i])
                exit(1)
            else:
                bondcoeff[i][0] = float(lines[ln+i].split()[2])
                bondcoeff[i][1] = float(lines[ln+i].split()[3])
    elif re.search("morse",bond_styles,flags=re.IGNORECASE): 
        print("Reading in morse bonds coefficents for",tbonds,"types")
        bond_style = 1 # morse bond style
        bondcoeff = numpy.zeros((tbonds,3))  # read in D alpha and r0
        ln = findline1(lines,"bond_coeff")
        for i in range(tbonds):
            btype = lines[ln+i].split()[1]
            if (int(btype) != i+1):
                print("Error wrong type in bond_coeff at line",ln+i)
                print(lines[ln+i])
                exit(1)
            else:
                bondcoeff[i][0] = float(lines[ln+i].split()[2])
                bondcoeff[i][1] = float(lines[ln+i].split()[3])
                bondcoeff[i][2] = float(lines[ln+i].split()[4])
    else:
        print ("Error bond_style \"",bond_styles,"\" not implemented")
        exit(1)
        
#-------------------------------------------------
def readinit(datafile): # read lammps init data file
    
    print ("Reading",datafile)
    fi = open(datafile,"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()

    global natoms
    ln = findline2(lines,"atoms")
    natoms = int(lines[ln].split()[0])

    global atypes
    ln = findline2(lines,"atom types")
    atypes = int(lines[ln].split()[0])

    global nbonds
    ln = findline2(lines,"bonds")
    nbonds = int(lines[ln].split()[0])

    global tbonds
    ln = findline2(lines,"bond types")
    tbonds = int(lines[ln].split()[0])

    global box
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

    print ("Natoms",natoms," Atypes",atypes," Bonds",nbonds," Btypes",tbonds)

    global mass
    global aatype
    global pos
    global vel
    global masses
    global bonds
    
    mass = numpy.zeros(atypes)
    aatype = numpy.zeros(natoms,dtype=int)
    pos = numpy.zeros((natoms,3))
    vel = numpy.zeros((natoms,3))
    masses = numpy.zeros((natoms,3))
    bonds =numpy.zeros((nbonds,3),dtype=int)
    
    ln = findline1(lines,"Masses")
    if (len(lines[ln+1].split())==0):
        ioff = ln+2
    else:
        ioff = ln+1
        
    for i in range(atypes):
        words = lines[ioff+i].split()
        if (int(words[0]) != i+1):
            print ("Error: while assigning masses at line",ln)
            print ("Expecting",i+1,"got",words[0],"from",lines[ln+2])
            exit(1)
        mass[i] = float(words[1])
    print ("Assigned masses")

    ln = findline1(lines,"Atoms")
    if (len(lines[ln+1].split())==0):
        ioff = ln+2
    else:
        ioff = ln+1

    for i in range(natoms):
        words = lines[ioff+i].split()
        if (int(words[0]) != i+1):
            print ("Error: while assigning atoms at line",ln)
            print ("Expecting",i+1,"got",words[0],"from",lines[ln+2])
            exit(1)
        aatype[i] = int(words[1])
        pos[i][0] = float(words[4])
        pos[i][1] = float(words[5])
        pos[i][2] = float(words[6])
        masses[i][0] = mass[aatype[i]-1] 
        masses[i][1] = mass[aatype[i]-1] 
        masses[i][2] = mass[aatype[i]-1] 
    print ("Assigned atom types and positions")

    ln = findline1(lines,"Velocities")
    if (ln!=-1):
        if (len(lines[ln+1].split())==0):
            ioff = ln+2
        else:
            ioff = ln+1
        for i in range(natoms):
            words = lines[ioff+i].split()
            if (int(words[0]) != i+1):
                print ("Error: while assigning velocities at line",ln)
                print ("Expecting",i+1,"got",words[0],"from",lines[ln+2])
                exit(1)
            vel[i][0] = float(words[1])
            vel[i][1] = float(words[2])
            vel[i][2] = float(words[3])
        print ("Assigned Velocities")
    else:
        print ("Velocities will start at zero")

    ln = findline1(lines,"Bonds")
    if (ln!=-1):
        if (len(lines[ln+1].split())==0):
            ioff = ln+2
        else:
            ioff = ln+1

        for i in range(nbonds):
            words = lines[ioff+i].split()
            if (int(words[0]) != i+1):
                print("Error: while assigning velocities at line",ln)
                print("Expecting",i+1,"got",words[0],"from",lines[ln+2])
                exit(1)
            bonds[i][0] = int(words[1])-1   # type
            bonds[i][1] = int(words[2])-1  # atom 1
            bonds[i][2] = int(words[3])-1  # atom 2
        print("Assigned Bonds")
    else:
        print("No bonds found")

def zero_momentum():  # zero the liniar momentum
    global masses, vel  
    mom = masses*vel
    tmom = numpy.sum(mom,axis=0)/numpy.sum(masses,axis=0) # total mom/mass
    vel -= tmom # zero out
    
#-----------------------------------------------------------
def bond_harm(): # harmonic bondk = 2 r0 = 1.5
    #Harmonic	k*(r-r0)^2
    #d/dr	2*k*(r-r0)
    #d^2/dr^2   2*k
	
    global pbond
        
    for i in range(nbonds):  # loop over bonds, 
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        bondk = bondcoeff[bonds[i][0]][0] # use type to bond params
        r0 = bondcoeff[bonds[i][0]][1]
        
        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0

        pot = bondk*dr**2
        dudr = 2.*bondk*dr
        # du2dr2 = 2.*bondk
        pbond += pot             # total bond

        dpos = (dudr/r)*dpos
        acc[bonds[i][1]] += dpos  # add forces back 
        acc[bonds[i][2]] -= dpos

#-----------------------------------------------------------
def bond_morse(): # Morse harmonic bondk = 2 r0 = 1.5
    # Morse EQ	D*(1-exp(-a(r-r0)))^2
    # d/dr	2 D a exp(-a(r-r0))*(1-exp(-a(r-r0)))
    # d^2/dr^2	4 a^2 D exp(-2 a (r-r0)) - 2 a^2 D exp(-a (r-r0))

    global pbond
        
    for i in range(nbonds):  # loop over bonds, 
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        # use type to bond params
        D = bondcoeff[bonds[i][0]][0]
        alpha = bondcoeff[bonds[i][0]][1] 
        r0 = bondcoeff[bonds[i][0]][2]
        
        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0
        
        expar = math.exp(-alpha*dr)
        pot = D*(1-expar)**2
        pbond += pot             # total bond

        dudr = 2.*D * alpha * expar *(1-expar)
        dpos = (dudr/r)*dpos
        
        acc[bonds[i][1]] += dpos  # add forces back 
        acc[bonds[i][2]] -= dpos

#-----------------------------------------------------------
def force():
    global acc
    global masses
    global pbond
    acc.fill(0) # zero out forces/acceration
    # lj
    # bonds
    pbond = 0
    if bond_style == 0:
        bond_harm()
    elif bond_style == 1:
        bond_morse()
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

zero_momentum()  # zero the momentum
acc = numpy.zeros((natoms,3)) # allocate and zero the acceleration
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
