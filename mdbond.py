# mdbond.py  # bonding routines for mdpython

import numpy
import math

#-------------- # harmonic bonds
def bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc):

    pbond = 0
    for i in range(nbonds):  # loop over bonds,
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

#--------------Numerical forces  dU/dx -------------------------------
def bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc):

    pbond = 0
    
#    print("Analytical forces",pos.size)
#    print(acc)
    print("Calculating Numerical Forces")
    acc_num = numpy.zeros(pos.size)
    acc_d = numpy.zeros(acc.shape)
    h = .000001

    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    # uses df(x)/dx = (f(x+h)-f(x-h))/(2h)
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)

    for i in range(post.size):
        tpos = post[i]
        post[i] = tpos + h
        pbond1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos - h
        pbondm1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos
        
        acc_num[i] = -(pbond1-pbondm1)/(2.0*h)
        #print(i,pbond,pbond1,pbondm1)
        
    numpy.copyto(acc,acc_num.reshape(-1,3))
    return pbond

#--------------------INM for harmonic-----------------------------
def inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian):

    # create hessian
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
            D = bondcoeff[bonds[i][0]][0] #idx D alpha r0
            alpha = bondcoeff[bonds[i][0]][1]
            r0 = bondcoeff[bonds[i][0]][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)

            dudr = 2.0*D * alpha * expar * (1.0-expar)
            du2dr2 = (2.0*D*alpha*alpha) * ( 2*expar*expar - expar)
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
            hessian[j][i] = hessian[i][j] # off diagonal

#-----------------------Numerical second derivative --------------------
def bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses,hess_num):

    pbond = 0
    h = .0001
    ma = masses.reshape(pos.size)
    acc_d = numpy.zeros(acc.shape)

    print("Calculating Numerical 2nd derivatives")
    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    for i in range(post.size):
        tpos = post[i]
        # uses d^2f(x)/ dx dx = (f(x+h,y)+f(x-h,y)-2f(x,y))/(hh)
        pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos - h
        pbondm1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos + h
        pbond1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos

        hess_num[i][i] = (pbond1+pbondm1-2.0*pbond)/(h*h*ma[i]) #diagonal with mass weight
        
        #print(i,pbond1,pbondm1,pbond,h,ma[i])
        #print(i,pbond,hess_num[i][i],post)
        for j in range(i+1,post.size):  # off diagonals
            # uses d^2f(x)/ dx dy = (f(x+h,y+h)+f(x-h,y-h)-f(x+h,y-h)-f(x-h,y+h))/(4hh)
            tposi = post[i]
            tposj = post[j]
            post[i] = tposi - h
            post[j] = tposj - h
            pbondm1m1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[j] = tposj + h
            pbondm11 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[i] = tposi + h
            pbond11 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[j] = tposj - h
            pbond1m1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[i] = tposi
            post[j] = tposj

            hess_num[i][j] = (pbond11+pbondm1m1-pbondm11-pbond1m1)/(h*h*4.0)
            hess_num[i][j] /= math.sqrt(ma[i]*ma[j])
            hess_num[j][i] = hess_num[i][j]

#-------------------------------------------------
def bond(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    #print("Calculating bond forces")
    pbond = 0
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)

    return pbond
    # check routines...
    #check_forces(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    #print(acc)
    #check_inm(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    #exit(1)

#----------------------------------------------------------
def bond_hess(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian):
    inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)

#----------------------------------------------------------
def check_forces(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):
    #check forces

    acc.fill(0) # rezero forces
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)

    acc_num = numpy.copy(acc)
    bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc_num)
    
    tol = 1e-6
    diff = acc_num-acc
    mv = max(diff.max(),abs(diff.min()))
    
    if(mv < tol):
        print("Forces Match, pbond =",pbond,mv)
    else:
        print("Forces DO NOT Match!")
        print("Analytical")
        print(acc)
        print("Numerical")
        print(acc_num)
        print("Diff = ",diff)

#----------------------------------------------------------
def check_inm(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    print("Calculating Hessian")
    hessian = numpy.zeros((pos.size,pos.size))
    bond_hess(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)
    #print(hessian)
            
    #print(pos.size)
    hess_num = numpy.zeros((pos.size,pos.size))
    bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses,hess_num)
    #print(hess_num)
    
    hdiff = hess_num-hessian
    mv = max(hdiff.max(),abs(hdiff.min()))/hessian.max()
    
    rdiff = numpy.sum(numpy.sum(hdiff*hdiff,axis=1))
    tol = 1e-5
    print(rdiff,mv)
    if(rdiff<tol and mv < tol):
        pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)
        print("Hessians Match, pbond =",pbond,mv,rdiff)
    else:
        print("Hessians DO NOT Match!")
        print(hessian)
        print(hess_num)
        print("Differences!")
        print("Diff = ",mv,rdiff,hdiff)
        print("Bummer!")
        exit(1)
    
    print("omega-squared")
    for i in range(nbonds):  # loop over bonds,
        itype = bonds[i][0]  # bond type
        mi = masses[bonds[i][1]][0]
        mj = masses[bonds[i][2]][0]
        mu =mi*mj/(mi+mj) # reduced mass
        if(bond_style==0): #Harmonic
            k0 = bondcoeff[itype][0] # use type to bond params
            r0 = bondcoeff[itype][1]
            print("Harmonic",2*k0/mu)
        elif(bond_style==1): #Morse:
            D = bondcoeff[bonds[i][0]][0] #idx D alpha r0
            alpha = bondcoeff[itype][1]
            r0 = bondcoeff[itype][2]
            ipos = pos[bonds[i][1]]
            jpos = pos[bonds[i][2]]
            dpos = jpos-ipos
            r =  math.sqrt(numpy.dot(dpos,dpos))
            dr = r-r0
            expar = math.exp(-alpha*dr)
            print("Morse",(2*alpha*alpha*D*(-1*expar+2*expar*expar))/mu)
        else: #??
            print("ERROR in bond_style")
            exit(1)
        
    w,v = numpy.linalg.eig(hessian)
    print("eigenvalus:",w)
    #print(v)

    w,v = numpy.linalg.eig(hess_num)
    print("eigenvalues num:", w)
    #print(v)
    
    print("Eigenvectors")
    for i in range(pos.size):
        print("")
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
