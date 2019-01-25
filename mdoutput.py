def write_dump(dumpfile,istep,natoms,pos):
    global fowrite
    
    try: fowrite  # has fowrite be assigned yet?
    except NameError:
        fowrite = open(dumpfile,"w")

    fowrite.write("%d\n"% natoms)
    fowrite.write("# step %d\n" % istep)
    for i in range(natoms):
        line = "Ar " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        fowrite.write(line)
    # fo.close()
    return 0
