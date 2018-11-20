"""Calulate effective resistance between grey region pairs"""

import pandas as pd
import random
import numpy as np
import sys
import networkx as nx
from scipy import linalg as LA
import h5py
import pickle

class white:
    index = 0
    def __init__(self, _index):
        self.index = _index
        self.connections = []
        self.possibility = []

# subj = sys.argv[1]
subj = 100206
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
with open("./WC_male/"+str(subj)+".txt","r") as f:
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

Grey=set(toGrey.values())

def test_wm_sym(): 
    """check if every connection in wm is undirected
    Result: False...the connections between voxels are directed."""
    for wm_idx in range(len(whitematter)):
        assert whitematter[wm_idx].index == wm_idx,'wm_idx mismatch '+str(wm_idx)
        for connection_i in range(len(whitematter[wm_idx].connections)):
            cur_connection = whitematter[wm_idx].connections[connection_i]
            if  cur_connection!= -1:
                assert whitematter[cur_connection].index == cur_connection,\
                    'cur_connection index mismatch '+str(wm_idx)+','+str(cur_connection)
                if wm_idx not in whitematter[cur_connection].connections:
                    print(wm_idx)
                    print(cur_connection)
                    return False
    return True

def inverse_dict(di):
    rev_di = {}
    for k in di:
        if di[k] not in rev_di:
            rev_di[di[k]] = []
        rev_di[di[k]].append(k)
    return rev_di

def build_binned_wm(binsize=8):
    with open('ori2binned.pkl', 'rb') as f:
        ori2binned = pickle.load(f)
    binned2ori = inverse_dict(ori2binned)
    bin_num = len(set(ori2binned.values()))
    # transition matrix of binned blocks
    trans_mat = np.zeros((bin_num+1,bin_num+1))
    for wm_idx in range(len(whitematter)):
        assert whitematter[wm_idx].index == wm_idx,'wm_idx mismatch '+str(wm_idx)
        for connection_i in range(len(whitematter[wm_idx].connections)):
            cur_connection = whitematter[wm_idx].connections[connection_i]
            if cur_connection != -1:
                assert whitematter[cur_connection].index == cur_connection,\
                    'cur_connection index mismatch '+str(wm_idx)+','+str(cur_connection)
                trans_mat[ori2binned[wm_idx],ori2binned[cur_connection]] += \
                            whitematter[wm_idx].possibility[connection_i]
            else:
                trans_mat[ori2binned[wm_idx],-1] += \
                        whitematter[wm_idx].possibility[connection_i]
    # normalize each row (there are blocks w/o any connection)
    row_sum = trans_mat.sum(axis=1)
    for i in range(len(trans_mat)):
        if row_sum[i] != 0:
            trans_mat[i] /= row_sum[i]
    return trans_mat

def build_wm_graph():
    trans_mat = build_binned_wm()
    tmp = trans_mat[:-1,:-1]
    tmp_A = (tmp+tmp.T)/2
    print(tmp_A)
    G = nx.Graph()
    G.add_nodes_from(np.arange(len(trans_mat)-1))
    for blk_i in range(len(trans_mat)-1):
        for blk_j in range(len(trans_mat)-1):
            if trans_mat[blk_i,blk_j]!=0 or trans_mat[blk_j,blk_i]!=0:
                G.add_edge(blk_i,blk_j,
                    weight=(trans_mat[blk_i,blk_j]+trans_mat[blk_j,blk_i])/2)
    return G

def er(G):
    # L = nx.laplacian_matrix(G).todense()
    A = nx.adjacency_matrix(G).todense()
    D = np.diag(np.squeeze(np.asarray(np.sum(A,axis=1))))
    L = D-A
    L_ = LA.pinv(L) #laplacian peudo inverse
    d,_  = L_.shape
    er_mat = np.zeros((d,d)) #TODO <OOM>
    for i in range(d):
        for j in range(d):
            if i!=j:
                er_mat[i,j] = er_mat[j,i] = round(L_[i,i]+L_[j,j]-2*L_[i,j],5)
    # save file
    with h5py.File('er.h5', 'w') as hf:
        hf.create_dataset("er", data=er_mat)
    return er_mat

def test_er(test_load=False):
    # series
    G = nx.Graph()
    G.add_edge(0,1)
    G.add_edge(1,2)
    G.add_edge(1,3)
    print(er(G))
    G.clear()
    # parallel
    G.add_edge(0,1)
    G.add_edge(0,2)
    G.add_edge(1,3)
    G.add_edge(2,3)
    print(er(G))
    # add every STOP node as different nodes
    G.add_edge(0,len(G.nodes()))
    G.add_edge(1,len(G.nodes()))
    print(er(G))
    G.clear()
    # add edge weights
    G.add_edge(0,1,weight=0.25)
    G.add_edge(0,2,weight=0.25)
    G.add_edge(1,3,weight=0.25)
    G.add_edge(2,3,weight=0.75)
    # G.add_edge(0,len(G.nodes()),weight=0.5)
    # G.add_edge(1,len(G.nodes()),weight=0.5)
    print(er(G))
    G.clear()
    # add STOP node as 1 node in a directed graph (doesnt work)
    G = nx.DiGraph()
    G.add_edge(0,-1)
    G.add_edge(0,1); G.add_edge(1,0)
    G.add_edge(0,2); G.add_edge(2,0)
    G.add_edge(1,-1)
    G.add_edge(1,3); G.add_edge(3,1)
    G.add_edge(2,3); G.add_edge(3,2)
    print(digrpah_er(G))
    G.clear()

    if test_load:
        with h5py.File('er.h5', 'r') as hf:
            data = hf['er'][:]
        print(data)

# test_er()

G = build_wm_graph()
print('graph built.')
er = er(G)[:len(whitematter),:len(whitematter)]