import sys
import pandas as pd
import gzip
import scipy as sp
import scipy.io
import numpy as np
import nibabel as nib
import random 
import operator 
import math
'''import matplotlib.pyplot as plt'''
df_base = './mda_fibs/%s/warped_mda.h5'
load_df = lambda subj: pd.read_hdf(df_base % subj, str(subj))
# subj = sys.argv[1]
subj = '180129'
scan = load_df(subj)
fibfile_base = './mda_fibs/%s/warped_mda.fib.gz'
fibfile = fibfile_base % subj
with gzip.open(fibfile, 'rb') as fibf:
    m = sp.io.matlab.loadmat(fibf)
    odf_vertices = m['odf_vertices'].T

num_peaks = 4
idxs = ['index%d' % n for n in range(num_peaks)]
fas = ['fa%d' % n for n in range(num_peaks)]

vectors = odf_vertices[scan[idxs].values.astype(np.int),:]  #[#wm,4,3]
magns = scan[fas].values[:,:,np.newaxis]    #[#wm,4,1]
mask_file = './voxel_clustering_MDA/MNI_mask_MC.nii.gz'
mask = nib.load(mask_file).get_data() > 0
idx2voxelIJK = np.array(np.where(mask)).T   ##[#wm,3]

'''
atlas_loc = './ROI_scale60.nii.gz'
atlas = nib.load(atlas_loc).get_data()
w=atlas.shape[0]
l=atlas.shape[1]
h=atlas.shape[2]

class grey:
    index = 0
    def __init__ (self,_index):
        self.index=_index
        self.connections= set([])

#g_list = [grey(x) for x in xrange(0,130)]
g_list=[]
for i in xrange(130):
    g_list.append(grey(i))
'''
whitedict= dict([((loc[0],loc[1],loc[2]),i) for i, loc in enumerate(idx2voxelIJK)])
#greymatter = [(i,j,k) for i in range(w) for j in range(l) for k in range(h) \
              #if atlas[i][j][k] > 0 and atlas[i][j][k] < 130]
#greymatter=set(greymatter)
whitematter=set(whitedict.keys())

'''
for g in greymatter:
    if (g in whitematter):
        g_list[atlas[g[0]][g[1]][g[2]]].connections.add(whitedict[g])
file = open("grey_connections","w")
for g in range(1,len(g_list)):
    print(g_list[g].index)
    file.write("%s %s\n" % (g_list[g].index,len(g_list[g].connections)) )
    for item in g_list[g].connections:
        file.write("%s\n" % item) 
file.close()
'''
class white:
    index = 0
    def  __init__(self, _index):
        self.index = _index
        self.connections = []
        self.possibility = []
        self.vectors = []
        self.magn = []

w_list=[white(i) for i in range(len(scan))]
Limit=100
for w in whitematter:
    index=whitedict[w]
    for m in range(3):
        if (abs(magns[index][m][0])>=4.36984267e-02):
            w_list[index].vectors.append(vectors[index][m])
            w_list[index].vectors.append(-vectors[index][m])
            w_list[index].magn.append(abs(magns[index][m][0])/2)
            w_list[index].magn.append(abs(magns[index][m][0])/2)
    base = sum(w_list[index].magn)
    w_list[index].possibility=[m/base for m in w_list[index].magn]
    for v in range(len(w_list[index].vectors)):
        cur_loc=np.array(w)
        step=0
        found=False
        while (step < Limit):
            step+=1
            cur_loc = cur_loc + w_list[index].vectors[v]*1.5
            next_loc = [math.floor(i) for i in cur_loc]
            if (next_loc !=w): 
                if(tuple(next_loc) in whitematter):
                    found=True
                    break
                else:
                    found=False
                    break     
        if (not found):
            w_list[index].connections.append(-1)
        else:
            w_list[index].connections.append(whitedict[tuple(next_loc)])

import ipdb; ipdb.set_trace()
file = open("/share/igert/share/DSI_data/WC_male/"+subj+".txt","w+")
for w in w_list:
    file.write("%s %s\n" % (w.index,len(w.connections)))
    for item in w.connections:
        file.write("%s\n" % item) 
    for item in w.possibility:
        file.write("%s\n" % item) 
file.close()
    
print("done")
        
#random walk
'''walks=100000
cur_loc=random.sample(g_list[1].connections,1)[0]
print(cur_loc)
found=False
step=1;
while (step<=walks):
    v1,v2,v3,v4=magns[cur_loc][0][0],magns[cur_loc][1][0],magns[cur_loc][2][0],magns[cur_loc][3][0]
    base=v1+v2+v3+v4
    p=[v1/base,v2/base,v3/base,v4/base]
    dir=tuple(vectors[cur_loc][np.random.choice(3,1,p)[0]])
    cur_crdnt=(idx2voxelIJK[cur_loc][0],idx2voxelIJK[cur_loc][1],idx2voxelIJK[cur_loc][2])
    pos=cur_crdnt
    new_crdnt=pos
    while (new_crdnt==cur_crdnt or new_crdnt not in whitematter):    
        pos=tuple(map(operator.add,pos,dir))
        new_crdnt=(math.floor(pos[0]),math.floor(pos[1]),math.floor(pos[2]))
        step+=1
        if (step>walks): break
    cur_loc=whitedict[new_crdnt]
    print(cur_loc)
    if (atlas[new_crdnt[0]][new_crdnt[1]][new_crdnt[2]]>1):
        found=True
        break
if (found): print("I found a path from 1 to "+str(atlas[new_crdnt[0]][new_crdnt[1]][new_crdnt[2]]))
else: print("I didn't find any path")
'''    


'''print 'Dimension of vectors:', vectors.shape, '\n'
print vectors'''
