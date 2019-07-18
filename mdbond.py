# mdbond.py  # bonding routines for mdpython

import numpy
import math

kb = 1.38064852e-23
w = [1.0/(2.0 - 2.0**(1/3)),0,0]
w[2] = w[0]
w[1] = 1.0 - 2.0*w[0]

def bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc):

    pbond = 0
    for i in range(nbonds): 
        itype = bonds[i][0]
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]

        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))

        if(bond_style==0):  # Harmonic k*(r-r0)^2, d/dr=2*k*(r-r0), d^2/dr^2=2*k
            bondk = bondcoeff[itype][0] # use type to bond params
            r0 = bondcoeff[itype][1]
            dr = r-r0
            pot = bondk*dr*dr
            dudr = 2.*bondk*dr
            dpos = (dudr/r)*dpos
        elif (bond_style==1): # Morse
            D = bondcoeff[itype][0]
            alpha = bondcoeff[itype][1]
            r0 = bondcoeff[itype][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)
            exparm1 = 1-expar
            pot = D*exparm1*exparm1
            dudr = 2.0*D*exparm1*alpha*expar
            dpos = (dudr/r)*dpos
        else:
            print ("Error in bond_style? in routine bond_force\n")
            exit(1)

        pbond += pot             # sum bond potential

        acc[bonds[i][1]] += dpos  # add forces to particles
        acc[bonds[i][2]] -= dpos

    return(pbond)

#--------------------INM -----------------------------
def inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian):


    for i in range(nbonds):  # loop over bonds,
        itype = bonds[i][0]  # bond type (in this case harmonic
        idx = bonds[i][1] # which atoms are involved in this bonds
        jdx = bonds[i][2]
        posi = pos[idx]
        posj = pos[jdx]

        rv = posj-posi
        r2 = numpy.dot(rv,rv)
        r = math.sqrt(r2)

        # print (itype,idx,jdx,posi,posj,k,r0,r)

        # d^2/ dx0 dxi
        idx3 = idx*3
        jdx3 = jdx*3

        #d^2 /dx^2 U(r(x,y)) = r" U' + r'^2 U"
        #d^2/ dx dy U(r(x,y)) = d^2/dxdy r dU/dr + dr/dx dr/dy d^2 U/dr^2
        if(bond_style==0):  # Harmonic
            k0 = bondcoeff[itype][0] # use type to bond params
            r0 = bondcoeff[itype][1]
            dudr = 2*k0*(r-r0)
            du2dr2 = 2*k0
            #print("Harmonic stuff")
            #print(k0,r0,dudr,du2dr2)

        if(bond_style==1): #Morse
            D = bondcoeff[itype][0] #idx D alpha r0
            alpha = bondcoeff[itype][1]
            r0 = bondcoeff[itype][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)
            dudr = 2.0*D*alpha*expar*(1.0-expar)
            du2dr2 = (2.0*D*alpha*alpha)*(2*expar*expar - expar)
            #print("Morse Stuff")
            #print(D,alpha,r0,dudr,du2dr2)

        rr = 1./r
        r3 = 1./(r*r*r)
        for k in range(3): # x y z of particle with index idx
            ii = idx*3
            jj = jdx*3
            diagelm = dudr*rr
            hessian[ii+k][ii+k] += diagelm
            hessian[ii+k][jj+k] -= diagelm
            hessian[jj+k][jj+k] += diagelm
            for l in range(3): # x y z of particle with index jdx
                elmij = -(rv[k]*rv[l])*r3*dudr + du2dr2*rv[k]*rv[l]/r2
                hessian[ii+k][ii+l] += elmij
                hessian[ii+k][jj+l] -= elmij
                hessian[jj+k][jj+l] += elmij
#------
    # mass weight
    ma = masses.reshape(pos.size) # make it easy to mass weight
    for i in range(pos.size):
        hessian[i][i] /= ma[i] # diagonal elements
        for j in range(i+1,pos.size):
            mw = math.sqrt(ma[i]*ma[j])
            hessian[i][j] /= mw # off diagonal
            hessian[j][i] = hessian[i][j]

