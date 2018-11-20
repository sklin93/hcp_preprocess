for i in range(5):
    file = open("runall_female"+str(i)+".sh", "w+")
    with open("female"+str(i)+".txt","r") as f:
        line=f.readline()
        while (line!=""):
            file.write("python randWalk.py "+line)
            line=f.readline()
    file.close()