# mdbond.py  # bonding routines for mdpython

import numpy
import math

#-------------- # harmonic bonds
def bond_harm(nbonds,bonds,bondcoeff,pos,acc): 
    #Harmonic	k*(r-r0)^2
    #d/dr	2*k*(r-r0)
    #d^2/dr^2   2*k

    pbond = 0
    for i in range(nbonds):  # loop over bonds, 
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        bondk = bondcoeff[bonds[i][0]][0] # use type to bond params
        r0 = bondcoeff[bonds[i][0]][1]
        
        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0

        pot = bondk*dr**2
        pbond += pot             # total bond

        dudr = 2.*bondk*dr
        dpos = (dudr/r)*dpos
        # du2dr2 = 2.*bondk

        acc[bonds[i][1]] += dpos  # add forces back 
        acc[bonds[i][2]] -= dpos

    return(pbond)

#--------------------INM for harmonic-----------------------------
def inm(istep,natoms,masses,pos,pot,hessian):
    # create hessian
    #Harmonic    k*(r-r0)^2
    #d/dr    2*k*(r-r0)
    #d^2/dr^2   2*k

    pbond = 0

hessian[0][0] = 
hessian[0][1] = 
hessian[0][2] =
hessian[0][3] =
hessian[0][4] =
hessian[0][5] = 

hessian[1][0] = 
hessian[1][1] = 


#----------------Morse potential---------------------------
def bond_morse(nbonds,bonds,bondcoeff,pos,acc): 
    # Morse EQ	D*(1-exp(-a(r-r0)))^2
    # d/dr	2 D a exp(-a(r-r0))*(1-exp(-a(r-r0)))
    # d^2/dr^2	4 a^2 D exp(-2 a (r-r0)) - 2 a^2 D exp(-a (r-r0))

    pbond = 0
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

    return(pbond)

#-------------------------------------------------
def bond(bond_style,nbonds,bonds,bondcoeff,pos,acc):
    pbond = 0
    if bond_style == 0: # harmonic bonds
        pbond = bond_harm(nbonds,bonds,bondcoeff,pos,acc)
    elif bond_style == 1: # morse bonds
        pbond = bond_morse(nbonds,bonds,bondcoeff,pos,acc)
    else:
        print ("Error bond style not found!")
        exit(1)
    return pbond



    
    

