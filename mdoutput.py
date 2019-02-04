# output for mdpython
def write_dump(dumpfile,istep,natoms,pos,aatype):
    global fowrite
    
    try: fowrite  # has fowrite be assigned yet?
    except NameError:
        fowrite = open(dumpfile,"w")

    fowrite.write("%d\n"% natoms)
    fowrite.write("# step %d\n" % istep)
    for i in range(natoms):
        if(aatype[i]==1):
            line = "O " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        elif (aatype[i]==2):
            line = "C " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        else :
            print("Error type not found?!")
            exit(1)
        fowrite.write(line)
    # fo.close()
    return 0

#---------------------------------------------------------------------------
def write_init(initfile,istep,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype):

    foinit = open(initfile,"w")

    foinit.write("Lammps init file written by MDPython after %d\n"% istep)
    foinit.write("\n")
    foinit.write("%d atoms\n"% natoms)
    foinit.write("%d atom types\n"% atypes)
    foinit.write("%d bonds\n"% nbonds)
    foinit.write("%d bond types\n"% tbonds)
    foinit.write("\n")
    line = str(-box[0]/2) + " " + str(box[0]/2) + " xlo xhi\n"
    foinit.write(line)
    line = str(-box[1]/2) + " " + str(box[1]/2) + " ylo yhi\n"
    foinit.write(line)
    line = str(-box[2]/2) + " " + str(box[2]/2) + " zlo zhi\n"
    foinit.write(line)
    foinit.write("\n")
    foinit.write("Masses\n")
    foinit.write("\n")
    for i in range(atypes):
        line = str(i+1) + " " + str(mass[i]) + "\n"
        foinit.write(line)
    
    foinit.write("\n")
    foinit.write("Atoms # full style: index type group charge x y z")
    foinit.write("\n")
    for i in range(natoms):
        line = str(i+1) + " " + str(aatype[i]) + " 1 0 " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.write("Velocities # vx vy vz\n")
    foinit.write("\n")
    for i in range(natoms):
        line = str(i+1) + " " + str(vel[i][0]) + " " + str(vel[i][1]) + " " + str(vel[i][2]) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.write("Bonds # type i j\n")
    foinit.write("\n")
    for i in range(nbonds):
        line = str(i+1) + " " + str(bonds[i][0]+1) + " " + str(bonds[i][1]+1) + " " + str(bonds[i][2]+1) + "\n"
        foinit.write(line)
        
    foinit.write("\n")
    foinit.close()
    return 0
