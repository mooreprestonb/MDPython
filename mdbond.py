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
        pbond += pot             # total bond potential

        dudr = 2.*bondk*dr
        dpos = (dudr/r)*dpos
        # du2dr2 = 2.*bondk

        acc[bonds[i][1]] += dpos  # add forces back
        acc[bonds[i][2]] -= dpos

    return(pbond)

#--------------------INM for harmonic-----------------------------
def inm(istep,natoms,masses,pos,pot,hessian):

#0 1 2
#3 4 5
#6 7 8

    hessian[0][0] = k-((r0*k*(y0**2)+(z0**2))/(r**3))
    hessian[0][1] = r0*k*x0*y0/(r**3)
    hessian[0][2] = r0*k*x0*z0/(r**3)
    hessian[0][3] = r0*k*x0*y0/(r**3)
    hessian[0][4] = k-((r0*k*(x0**2)+(z0**2))/(r**3))
    hessian[0][5] = r0*k*z0*y0/(r**3)
    hessian[0][6] = r0*k*x0*z0/(r**3)
    hessian[0][7] = r0*k*z0*y0/(r**3)
    hessian[0][8] = k-((r0*k*(y0**2)+(x0**2))/(r**3))

    hessian[1][0] = k-((r0*k*(y1**2)+(z1**2))/(r**3)) 
    hessian[1][1] = r0*k*x1*y1/(r**3)
    hessian[1][2] = r0*k*x1*z1/(r**3)
    hessian[1][3] = r0*k*x1*y1/(r**3)
    hessian[1][4] = k-((r0*k*(x1**2)+(z1**2))/(r**3))
    hessian[1][5] = r0*k*z1*y1/(r**3)
    hessian[1][6] = r0*k*x1*z1/(r**3)
    hessian[1][7] = r0*k*z1*y1/(r**3)
    hessian[1][8] = k-((r0*k*(y1**2)+(x1**2))/(r**3))


#---------

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

#--------------Numerical forces  dU/dx -------------------------------
def bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc):

    pbond = 0
    print(pos.size)
    acc_num = numpy.zeros(pos.size)
    h = .00000001

    print("Numerical derivatives")
    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    # uses df(x)/dx = (f(x+h)-f(x-h))/(2h)
    for i in range(post.size):
        tpos = post[i]
        post[i] = tpos - h
        pbondm1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        #pbondm1 = bond_morse(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        post[i] = tpos + h
        pbond1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        #pbond1 = bond_morse(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        post[i] = tpos
        acc_num[i] = -(pbond1-pbondm1)/(h*2.0)

    print(acc_num.reshape(-1,3))
    return pbond

#--------------Numerical forces  dU/dx -------------------------------
def bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    pbond = 0
    print(pos.size)
    hess_num = numpy.zeros((pos.size,pos.size))
    h = .000001
    ma = masses.reshape(pos.size)

    print("Numerical 2nd derivatives")
    # uses d^2f(x)/ dx dy = (f(x+h,y+h)+f(x-h,y-h)-f(x+h,y-h)-f(x-h,y+h))/(4hh)
    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    for i in range(post.size):
        tpos = post[i]
        pbond = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        #pbondm1 = bond_morse(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        post[i] = tpos - h
        pbondm1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        #pbondm1 = bond_morse(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        post[i] = tpos + h
        pbond1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        #pbond1 = bond_morse(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
        post[i] = tpos
        hess_num[i][i] = (pbond1+pbondm1-2.0*pbond)/(h*h*ma[i]) #diagonal
        for j in range(i+1,post.size):  # off diagonals
            tposi = post[i]
            tposj = post[j]
            post[i] = tposi - h
            post[j] = tposj - h
            pbondm1m1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
            post[j] = tposj + h
            pbondm11 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
            post[i] = tposi + h
            pbond11 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)
            post[j] = tposj - h
            pbond1m1 = bond_harm(nbonds,bonds,bondcoeff,post.reshape(-1,3),acc)

            hess_num[i][j] = (pbond11+pbondm1m1-pbondm11-pbond1m1)/(h*h*4.0)
            hess_num[i][j] /= math.sqrt(ma[i]*ma[j])
            hess_num[j][i] = hess_num[i][j]
            post[i] = tposi
            post[j] = tposj

    print(hess_num)
    w,v = numpy.linalg.eig(hess_num)
    print(w)
    print(v)
    return pbond

#-------------------------------------------------
def bond(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):
    pbond = 0
    if bond_style == 0: # harmonic bonds
        pbond = bond_harm(nbonds,bonds,bondcoeff,pos,acc)
    elif bond_style == 1: # morse bonds
        pbond = bond_morse(nbonds,bonds,bondcoeff,pos,acc)
    else:
        print ("Error bond style not found!")
        exit(1)
#check forces
    print(acc)
    bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc)
    bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    exit(1)
#endcheck
    return pbond
