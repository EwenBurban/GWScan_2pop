import pandas as pd
import sys
import random as rd
import re

file = sys.argv[1]
output_file = sys.argv[2]

## load data
data = pd.read_csv(file,sep='\t',comment='#',header=None)
## get the header
header=list()
pattern=re.compile('#')
with open(file,'r') as f:
    line = f.readline()
    while line :
        if pattern.search(line):
            header.append(line)
        line = f.readline()

## pick one haplotype 
def pick_haplo(array):
    choice = rd.randint(0,1)
    y = [ x.split('|')[choice] for x in array]
    return y

test = data.loc[1,].tolist()
i=0
for a in test:
    if '|' not in str(a):
        i+=1

row_start=data.iloc[:,0:i]
diplotypes = data.iloc[:,i:]
haplotypes=diplotypes.apply(pick_haplo)
final_data = pd.concat([row_start,haplotypes],1)
with open(output_file,'w') as o:
    for l in header:
        o.write(l)
final_data.to_csv(output_file,mode='a',index=False,header=False,sep='\t')

