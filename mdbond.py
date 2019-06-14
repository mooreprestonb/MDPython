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
def inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian):

    # create hessian
    # Harmonic    k*(r-r0)^2
    # d/dr    2*k*(r-r0)
    # d^2/dr^2   2*k

    pbond = 0

    for i in range(nbonds):  # loop over bonds,
        itype = bonds[i][0]  # bond type (in this case harmonic
        idx = bonds[i][1] # which atoms are involved in this bonds
        jdx = bonds[i][2]
        posi = pos[idx]
        posj = pos[jdx]


        rv = posj-posi
        r2 = numpy.dot(rv,rv)
        r = math.sqrt(r2)

#        print (itype,idx,jdx,posi,posj,k,r0,r)

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
            print("Harmonic stuff")
            print(k0,r0)

        if(bond_style==1): #Morse
            D = bondcoeff[bonds[i][0]][0] #idx D alpha r0
            alpha = bondcoeff[bonds[i][0]][1]
            r0 = bondcoeff[bonds[i][0]][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)
            dudr = 2.0*D * alpha * expar *(1.0-expar)
            du2dr2 = (2*D*alpha*alpha)*((-expar) + (2*expar*expar))
            print("Morse Stuff")
            print(D,alpha,r0)

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

#        x0 = ipos[0]
#        y0 = ipos[1]
#        z0 = ipos[2]
#        x1 = jpos[0]
#        y1 = jpos[1]
#        z1 = jpos[2]
         #r = math.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)

#        drdx = rv[0]/r
#        drdy = rv[1]/r
#        drdz = rv[2]/r
#        dr2dxz = -(rv[0]*rv[2])/r**3
#        dr2dyz = -(rv[1]*rv[2])/r**3
#        dr2dxy = -(rv[0]*rv[1])/r**3
#        dr2dx2 = (1/r)-(rv[0]*rv[0])/r**3
#        dr2dy2 = (1/r)-(rv[1]*rv[1])/r**3
#        dr2dz2 = (1/r)-(rv[2]*rv[2])/r**3

#        hessian[0][0] = dr2dx2*dudr + drdx*drdx*du2dr2
#        hessian[0][1] = dr2dxy*dudr + drdx*drdy*du2dr2
#        hessian[0][2] = dr2dxz*dudr + drdx*drdz*du2dr2
#        hessian[0][3] = -dr2dx2*dudr - drdx*drdx*du2dr2
#        hessian[0][4] = -dr2dxy*dudr - drdy*drdx*du2dr2
#        hessian[0][5] = -dr2dxz*dudr - drdz*drdx*du2dr2

#        hessian[1][1] = dr2dy2*dudr + drdy*drdy*du2dr2
#        hessian[1][2] = dr2dyz*dudr + drdy*drdz*du2dr2
#        hessian[1][3] = -dr2dxy*dudr - drdx*drdy*du2dr2
#        hessian[1][4] = -dr2dy2*dudr - drdy*drdy*du2dr2
#        hessian[1][5] = -dr2dyz*dudr - drdz*drdy*du2dr2

#        hessian[2][2] = dr2dz2*dudr + drdz*drdz*du2dr2
#        hessian[2][3] = -dr2dxz*dudr - drdx*drdz*du2dr2
#        hessian[2][4] = -dr2dyz*dudr - drdy*drdz*du2dr2
#        hessian[2][5] = -dr2dz2*dudr - drdz*drdz*du2dr2

#        hessian[3][3] = dr2dx2*dudr + drdx*drdx*du2dr2
#        hessian[3][4] = dr2dxy*dudr + drdy*drdx*du2dr2
#        hessian[3][5] = dr2dxz*dudr + drdz*drdx*du2dr2

#        hessian[4][4] = dr2dy2*dudr + drdy*drdy*du2dr2
#        hessian[4][5] = dr2dyz*dudr + drdy*drdz*du2dr2

#        hessian[5][5] = dr2dz2*dudr + drdz*drdz*du2dr2

#    print(hessian)
#------
    # mass weight
    ma = masses.reshape(pos.size) # make it easy to mass weight
    for i in range(pos.size):
        hessian[i][i] /= ma[i] # diagonal elements
        for j in range(i+1,pos.size):
            mw = math.sqrt(ma[i]*ma[j])
            hessian[i][j] /= mw # off diagonal
            hessian[j][i] = hessian[i][j] # off diagonal
    mu = ma[0]*ma[3]/(ma[0]+ma[3])

    print("omega-squared")
    if bond_style==0: #Harmonic
        print(2*k0/mu)
        #freq_lst = [k/m for m in mu]
        #unique_freqs = []
        #for freq in freq_lst:
        #    if freq not in unique_freqs:
        #        unique_freqs.append(freq)
        #print(unique_freqs)
    if bond_style==1: #Morse:
        expar = math.exp(-alpha*dr)        
        print((2*alpha*alpha*D*(-1*expar+2*expar*expar))/mu)
        #freq_lst = [alpha/m for m in mw]
        #unique_freqs = []
        #for freq in freq_lst:
        #    if freq not in unique_freqs:
        #        unique_freqs.append(freq)
        #print(unique_freqs)

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

        dudr = 2.0*D * alpha * expar *(1.0-expar)
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
    # uses d^2f(x)/ dx dx = (f(x+h,y)+f(x-h,y)-2f(x,y))/(4hh)
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
#    print(v)
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

    print(pos.size)
    hessian = numpy.zeros((pos.size,pos.size))
    inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)
    print("")
    print("Analytical Hessian")
    print(hessian)
    w,v = numpy.linalg.eig(hessian)
    #    print(w)
    #    print(v)
    for i in range(pos.size):
        print("")
        print("Eigenvalue", i)
        print(i,w[i],v[:,i])

    f = open("check.vmd","w")
    f.write("mol new\n")
    f.write("draw color blue\n")
    f.write("draw sphere {%f %f %f} radius .2\n" % (pos[0][0], pos[0][1], pos[0][2]))
    f.write("draw color red\n")
    f.write("draw sphere {%f %f %f} radius .2\n" % (pos[1][0], pos[1][1], pos[1][2]))
    for i in range(6):
        f.write("draw color %d\n" %(i+1))
        f.write("draw line {%f %f %f} {%f %f %f} width 3\n" % (pos[0][0], pos[0][1], pos[0][2], pos[0][0]+v[0][i],pos[0][1]+v[1][i], pos[0][2]+v[2][i]))
        f.write("draw line {%f %f %f} {%f %f %f} width 3\n" % (pos[1][0], pos[1][1], pos[1][2], pos[1][0]+v[3][i],pos[1][1]+v[4][i], pos[1][2]+v[5][i]))
    f.close()


    exit(1)
#endcheck
    return pbond
