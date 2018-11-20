import pandas as pd
import random
import numpy as np
import sys
class white:
    index = 0
    def __init__(self, _index):
        self.index = _index
        self.connections = []
        self.possibility = []

class path:
    logp=0.0
    length=0
    des=0
    def __init__(self, _length, _logp, _des):
        self.logp=_logp
        self.length=_length
        self.des=_des
subj = sys.argv[1]
toWhite={}
toGrey={}
whitematter=[]
with open("grey_connections","r") as f:
    line=f.readline()
    while (line!=""):
        line=line.split()
        index, length=int(line[0]), int(line[1])
        neighbours=[]
        for i in range(length):
            number = int(f.readline())
            toGrey[number]=index
            neighbours.append(number)
        toWhite[index]=neighbours
        line=f.readline()
with open("/share/igert/share/DSI_data/WC_male/"+subj+".txt","r") as f:
    line=f.readline()
    while (line!=""):
        line=line.split()
        index, length=int(line[0]), int(line[1])
        whitematter.append(white(index))
        for i in range(length):
            number = int(f.readline())
            whitematter[index].connections.append(number)
        for i in range(length):
            prob = float(f.readline())
            whitematter[index].possibility.append(prob)
        line=f.readline()

L=1000 #max walk length
correlation=[]
Grey=set(toGrey.values())

for g in range(1,len(Grey)+1):   
    path_list=[]
    for b in range(500):
        loc=random.sample(toWhite[g],1)[0]
        step=0
        logp=0
        while (step<L):
            step+=1
            if (len(whitematter[loc].connections)==0):
                loc=-1
            else:
                rand=np.random.choice(len(whitematter[loc].connections)-1,1,whitematter[loc].possibility)[0]
                logp+=np.log(whitematter[loc].possibility[rand])
                loc=whitematter[loc].connections[rand]
            if (loc == -1): #dead end, set length to max walk length
                path_list.append(path(L,logp,-1))
                break
            elif (loc in toGrey and toGrey[loc]!= g):
                path_list.append(path(step,logp,toGrey[loc]))
                break
            elif (step==L):
                path_list.append(path(L,logp,-1))
    for i in range(1,len(Grey)+1):
        sum=0
        for p in path_list:
            if (p.des==i):
                sum+=np.exp((1-p.length)+p.logp)
            else:
                sum+=np.exp((1-L)+p.logp)
        correlation.append(sum)
correlation=np.array(correlation).reshape((len(Grey),len(Grey)))
#print correlation
df = pd.DataFrame(correlation)
#df.columns=["Grey"+str(i) for i in range(1,len(Grey)+1)]
df.to_pickle("/share/igert/share/DSI_data/Matrix_male/"+subj+".pkl")
