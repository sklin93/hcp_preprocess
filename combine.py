import h5py
import numpy as np
from scipy.io import savemat

subj = [101107, 101410, 105216, 116524, 119833, 126628, 127630, 128127, 131217, 136227, 138231, 138534, 139233, 140925, 141422, 141826, 147030, 150625, 159239, 162733, 163129, 163836, 169444, 179548, 180836, 182739, 185139, 197348, 201414, 210617, 239944, 250932, 255639, 290136, 365343, 366446, 377451, 395958, 415837, 436845, 441939, 541943, 547046, 679568, 683256, 732243, 789373, 861456, 865363, 899885, 917255]
s = []
for i in range(len(subj)):
	fname = 'er_structural/'+str(subj[i])+'.h5'
	with h5py.File(fname,'r') as f: 
		s.append(f['S'][:])
s = np.stack(s)
data = {'subj':np.asarray(subj).reshape(1,51),'X':s.transpose(1,2,0)}
savemat('hcp_S.mat',data)
import ipdb; ipdb.set_trace()

from scipy.io import loadmat
data = loadmat('hcp_S.mat')
subj = data['subj']
s = data['X']