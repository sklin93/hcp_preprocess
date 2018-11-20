import sys
import pandas as pd
import os

file=pd.read_csv('HCP_900_demographics.csv')
subjects=[str(i) for i in file['Subject']]
gender=file['Gender']
hcp=set(subjects)
d=dict(zip(subjects,gender))
path='/share/igert/share/DSI_data/fs125/mda_fibs'
male = open('male.txt','w')
female = open('female.txt','w')
for filename in os.listdir(path):
    if filename in hcp:
        if d[filename]=='M':
            male.write("%s\n" % filename)
        else:
            female.write("%s\n" % filename)
   
male.close()
female.close()