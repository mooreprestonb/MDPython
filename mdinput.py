# routines for input MDPython

import re

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

#-------------------------------------------------
def readinvals(lines,data): # read lammps init data file
    
    # global natoms
    ln = findline2(lines,"atoms")
    data[0] = int(lines[ln].split()[0])

    # global atypes
    ln = findline2(lines,"atom types")
    data[1] = int(lines[ln].split()[0])

    # global nbonds
    ln = findline2(lines,"bonds")
    data[2] = int(lines[ln].split()[0])

    #global tbonds
    ln = findline2(lines,"bond types")
    data[3] = int(lines[ln].split()[0])

    # global box
    # box = numpy.zeros(3)
    ln = findline2(lines,"xlo xhi")
    xlo = float(lines[ln].split()[0])
    xhi = float(lines[ln].split()[1])
    data[4] = xhi-xlo

    ln = findline2(lines,"ylo yhi")
    ylo = float(lines[ln].split()[0])
    yhi = float(lines[ln].split()[1])
    data[5] = yhi-ylo

    ln = findline2(lines,"zlo zhi")
    zlo = float(lines[ln].split()[0])
    zhi = float(lines[ln].split()[1])
    data[6] = zhi-zlo

    #print("Natoms",data[0]," Atypes",data[1]," Bonds",data[2]," Btypes",data[3])
    #print("Box",data[4],data[5],data[6])

#----------------------------------------------------------------
def getmasses(lines,atypes,mass):
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

#----------------------------------------------------------------
def getatoms(lines,natoms,aatype,pos,mass,masses):
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

#----------------------------------------------------------------
def getvel(lines,natoms,vel):
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
#----------------------------------------------------------------
def getbonds(lines,nbonds,bonds):
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
            bonds[i][0] = int(words[1])-1  # type
            bonds[i][1] = int(words[2])-1  # atom 1
            bonds[i][2] = int(words[3])-1  # atom 2
        print("Assigned Bonds")
    else:
        print("No bonds found!")
        
#----------------------------------------------------------------
def getbondcoeff(lines,bond_styles,tbonds,bondcoeff):
    if re.search("harmonic",bond_styles,flags=re.IGNORECASE): 
        print("Reading in harmonic bonds coefficents for",tbonds,"types")
        bond_style = 0 #harmonic bond style
        ln = findline1(lines,"bond_coeff")
        for i in range(tbonds):
            if(lines[ln+i].split()[0] != "bond_coeff"):
                print("Error not enough bond_coeff's")
                exit(1)
            btype = lines[ln+i].split()[1]
            if (int(btype) != i+1):
                print("Error wrong type in bond_coeff at line",ln+i)
                print(lines[ln+i])
                exit(1)
            else:
                bondcoeff[i][0] = float(lines[ln+i].split()[2])
                bondcoeff[i][1] = float(lines[ln+i].split()[3])
        if(len(lines[ln+i+1].split())>0):
            if (lines[ln+i+1].split()[0] == "bond_coeff"):
                print("Error: too MANY bond_coeffs for types")
                print(i,tbonds)
                exit(1)

    elif re.search("morse",bond_styles,flags=re.IGNORECASE): 
        print("Reading in morse bonds coefficents for",tbonds,"types")
        bond_style = 1 # morse bond style
        # bondcoeff = numpy.zeros((tbonds,3))  # read in D alpha r0
        ln = findline1(lines,"bond_coeff")
        for i in range(tbonds):
            if(lines[ln+i].split()[0] != "bond_coeff"):
                print("Error not enough bond_coeff's")
                exit(1)
            btype = lines[ln+i].split()[1]
            if (int(btype) != i+1):
                print("Error wrong type in bond_coeff at line",ln+i)
                print(lines[ln+i])
                exit(1)
            else:
                bondcoeff[i][0] = float(lines[ln+i].split()[2])
                bondcoeff[i][1] = float(lines[ln+i].split()[3])
                bondcoeff[i][2] = float(lines[ln+i].split()[4])
        if(len(lines[ln+i+1].split())>0):
            if (lines[ln+i+1].split()[0] == "bond_coeff"):
                print("Error: too MANY bond_coeffs for types")
                print(i,tbonds)
                exit(1)
    else:
        print ("Error bond_style \"",bond_styles,"\" not implemented")
        exit(1)

    return bond_style
#----------------------------------------------------------------
def readsysvals(lines,data): # read lammps init data file

    # global nsteps
    ln = findline1(lines,"run")
    data[0] = int(lines[ln].split()[1])

    # global dt
    ln = findline1(lines,"timestep")
    data[1] = float(lines[ln].split()[1])

    # global initfile
    ln = findline1(lines,"read_data")
    data[2] = (lines[ln].split()[1])

    # global ithermo
    ln = findline1(lines,"thermo")
    data[3] = int(lines[ln].split()[1])

    #global idump
    # global dumpfile
    ln = findline1(lines,"dump")
    data[4] = int(lines[ln].split()[4])
    data[5] = lines[ln].split()[5]

    # global bond_style
    ln = findline1(lines,"bond_style")
    data[6] = lines[ln].split()[1]

    # global logfile
    ln = findline1(lines,"log")
    if(ln!=-1):
        data[7] = lines[ln].split()[1]

