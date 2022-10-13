# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 17:00:05 2022

@author: wilwin
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy as np

a = []
b=[]
c =[]
d =[]
#print(a)
mat_size = 54 #30 d-orbitals of Fe (5 4d-orbitals * 6 atoms) 
              #and  24 orbitals of Sn (4 sp-orbitals *6 atoms)
mat_ind = []
column_wf = []
i=0
j=0
k=0
lines=[]
file = open('out_fit.wandef', 'r')
for line in file:
    lines.append(line)
label=[]
    
for i in range(len(lines)-1):
    #print('LINE = ', line)
    text = lines[i]
    next_text = lines[i+1]
    if 'spin 1: ' in text:
        WF1 = text.split('spin 1: WF(', 1)
        text = WF1[1]
        WF1 = text.split(') -> WF(')
        text = WF1[1]
        WF1 = WF1[0]
        if WF1 not in label:
            label.append(WF1)
        WF2 = text.split(') at relative')
        WF2 = WF2[0]
        mat_ind.append([WF1,WF2])
        i+=1
        #print(i)
        
        if 'hop' in next_text:
            hop = next_text.split('hop=', 1)
            next_text = hop[1]
            hop = next_text.split('+i*', 1)
            real_hop = hop[0]
            im_hop = hop[1]
            j+=1
        else:
            real_hop = 0
            im_hop = 0
        #if abs(float(real_hop))<10:
        a.append((float(real_hop)))
        c.append((float(im_hop)))
        #else:
        #    a.append((float(100)))
        #    c.append((float(100)))            
         
        if len(a)==mat_size:
            b.append(a)
            d.append(c)
            a=[]
            c=[]
        
print(b)
np.diag(b)
idx=np.argsort(np.diag(b))
print(idx)
b_sort=[]
label_sort=[]
for i in idx:
    bs=[]
    label_sort.append(label[i])
    for j in idx:
        bs.append(b[i][j])
    b_sort.append(bs)

plt.rcParams["figure.figsize"] = (50,40)
im = plt.pcolormesh(b_sort,cmap="PiYG")
ax = plt.gca() 
cbar = plt.colorbar(im)
ax.set_xticks(np.arange(len(b))+0.5)
ax.set_yticks(np.arange(len(b))+0.5)
plt.clim(-6, 6)
ax.set_xticklabels(label_sort,rotation=90)
ax.set_yticklabels(label_sort)
cbar.ax.tick_params(labelsize=80) 
#ax.set_aspect('equal')
plt.show()

im = plt.pcolormesh(d, edgecolors='k', linewidth=2)
ax = plt.gca() 
plt.colorbar(im)
ax.set_aspect('equal')
plt.show()