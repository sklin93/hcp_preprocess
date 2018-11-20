import pandas as pd
import numpy as np
import os 
import sys
def dist(a,b):
    #return np.abs(np.sqrt(np.sum(m**2))-np.sqrt(np.sum(mean**2)))
    #return np.sqrt(np.sum((m-mean)**2)) 
    a_mean=a.sum()/a.size
    b_mean=b.sum()/b.size
    return ((a-a_mean)*(b-b_mean)).sum()/(np.sqrt(((a-a_mean)**2)*((b-b_mean)**2))).sum()
def average(cluster):
    sum = np.zeros((129,129))
    for m in cluster:
        sum+=m.matrix
    return sum/len(cluster)
def membership(matrix,mean):
    min=1000000
    cluster = 0
    for i in range(len(mean)):
        distance = dist(matrix,mean[i])
        if distance < min:
            min=distance
            cluster = i
    return cluster
class subject:
    matrix=np.zeros((129,129))
    cluster=0
    gender=''
    def __init__(self,_matrix,_cluster,_gender):
        self.matrix = _matrix
        self.cluster = _cluster
        self.gender=_gender
#k-mean
k=int(sys.argv[1])
mean=[]
new_mean=[]
cluster=[]
for i in range(k):
    cluster.append([])
subj=[]
path="/share/igert/share/DSI_data/Matrix_female/"
for filename in os.listdir(path):
    m=pd.read_pickle(os.path.join(path+filename)).values
    new_subj=subject(m,0,'F')
    subj.append(new_subj)
    cluster[0].append(new_subj)
path="/share/igert/share/DSI_data/Matrix_male/"
for filename in os.listdir(path):
    m=pd.read_pickle(os.path.join(path+filename)).values
    new_subj=subject(m,0,'M')
    subj.append(new_subj)
    cluster[0].append(new_subj)
rand=np.random.choice(len(subj),k)
for i in rand:
    mean.append(subj[i].matrix)
mean=np.array(mean)
c = 0
stable = False
while (not stable): 
    c += 1
    stable = True
    for s in subj:
        new_cluster = membership(s.matrix,mean)
        if (s.cluster != new_cluster):
            cluster[s.cluster].remove(s)
            cluster[new_cluster].append(s)
            s.cluster=new_cluster
            stable=False    
    mean=np.array(list(map(lambda x:average(x), cluster)))  
    print c
M = [0 for i in range(len(cluster))]
F = [0 for i in range(len(cluster))]
for s in subj:
    if s.gender=='M':
        M[s.cluster]+=1
    else:
        F[s.cluster]+=1
for i in range(len(cluster)):
    print "Cluster", i
    print "M:", M[i]," F:",F[i]
num=0
for m in mean:
    df=pd.DataFrame(m)
    df.to_pickle("mean"+str(num)+".pkl")
    num+=1
        
    
    
    
    
    