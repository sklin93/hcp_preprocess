"""Calulate effective resistance between grey region pairs"""

import pandas as pd
import random
import numpy as np
import os,sys
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

subj = sys.argv[1]
# subj = 100206
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
with open('ori2binned.pkl', 'rb') as f:
        ori2binned = pickle.load(f)

Grey = set(toGrey.values())
scale60_to_scale33 = {1: 1, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 7, 10: 7, 11: 8, 12: 8, 13: 8, 14: 8, 15: 9, 16: 10, 17: 10, 18: 10, 19: 10, 20: 11, 21: 12, 22: 13, 23: 14, 24: 15, 25: 16, 26: 16, 27: 16, 28: 17, 29: 17, 30: 18, 31: 18, 32: 18, 33: 19, 34: 19, 35: 20, 36: 20, 37: 21, 38: 22, 39: 23, 40: 23, 41: 24, 42: 24, 43: 25, 44: 25, 45: 26, 46: 27, 47: 28, 48: 29, 49: 29, 50: 30, 51: 30, 52: 31, 53: 32, 54: 32, 55: 33, 56: 34, 57: 34, 58: 35, 59: 36, 60: 37, 61: 38, 62: 39, 63: 40, 64: 41, 65: 42, 66: 42, 67: 43, 68: 44, 69: 45, 70: 46, 71: 47, 72: 48, 73: 48, 74: 48, 75: 49, 76: 49, 77: 49, 78: 49, 79: 50, 80: 51, 81: 51, 82: 51, 83: 51, 84: 52, 85: 53, 86: 54, 87: 55, 88: 56, 89: 57, 90: 57, 91: 57, 92: 58, 93: 58, 94: 59, 95: 59, 96: 59, 97: 60, 98: 60, 99: 61, 100: 61, 101: 62, 102: 63, 103: 64, 104: 64, 105: 65, 106: 65, 107: 66, 108: 66, 109: 67, 110: 68, 111: 69, 112: 70, 113: 70, 114: 71, 115: 71, 116: 72, 117: 73, 118: 73, 119: 74, 120: 75, 121: 75, 122: 76, 123: 77, 124: 78, 125: 79, 126: 80, 127: 81, 128: 82, 129: 83}
Grey33 = set(scale60_to_scale33.values())

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

def build_binned_gm2wm(scale=60):
    if scale==60:
        binned_toGrey = {ori2binned[k]:v for k,v in toGrey.items()}
    if scale==33:
        binned_toGrey = {ori2binned[k]:scale60_to_scale33[v] for k,v in toGrey.items()}
    else:
        print('unsupported scale.')
        return
    binned_toWhite = inverse_dict(binned_toGrey)
    return binned_toGrey, binned_toWhite

def build_wm_graph():
    G = nx.Graph()
    G.add_nodes_from(np.arange(len(trans_mat)-1))
    for blk_i in range(len(trans_mat)-1):
        for blk_j in range(len(trans_mat)-1):
            if trans_mat[blk_i,blk_j]!=0 or trans_mat[blk_j,blk_i]!=0:
                G.add_edge(blk_i,blk_j,
                    weight=(trans_mat[blk_i,blk_j]+trans_mat[blk_j,blk_i])/2)
    return G

def er(X,from_A=False,save_er=False):
    if not from_A:
        A = nx.adjacency_matrix(X).todense()
    else:
        A = X
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
    if save_er:
        with h5py.File('er_structural/'+str(subj)+'_er.h5', 'w') as hf:
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
# G = build_wm_graph()
# er(G)
binned_toGrey, binned_toWhite = build_binned_gm2wm(scale=33)
trans_mat = build_binned_wm()[:-1,:-1]
A = (trans_mat+trans_mat.T)/2
wm_er = er(A,from_A=True)
# use 1/ER as edge strength
strength_wm = np.divide(1,wm_er)
'''
1/ER_gm = sum(1/ER_wm) for all wm paths connecting gm pairs,
thus strength_gm = 1/ER_gm = sum(strength_wm)
'''
strength_gm = np.zeros((len(Grey33),len(Grey33)))
for gm_i in range(len(Grey33)):
    for gm_j in range(len(Grey33)):
        i_toWhite = binned_toWhite[gm_i+1]
        j_toWhite = binned_toWhite[gm_j+1]
        for i_wm in i_toWhite:
            for j_wm in j_toWhite:
                strength_gm[gm_i,gm_j] += strength_wm[i_wm,j_wm]
# normalization
'''current normalization method: 2/(1+e^(-x^0.2))-1'''
normalized_strength = 2/(1+np.exp(-np.power(strength_gm,0.2)))-1
with h5py.File('er_structural/'+str(subj)+'.h5', 'w') as hf:
    hf.create_dataset("S", data=normalized_strength)
print('subject '+str(subj)+'finshed.')