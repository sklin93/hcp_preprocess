file = open("runall_male.sh", "w+")
with open("male.txt","r") as f:
    line=f.readline()
    while (line!=""):
        file.write("python randWalk.py "+line)
        line=f.readline()
file.close()