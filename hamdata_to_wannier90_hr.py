# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 17:00:05 2022

@author: wilwin
"""

import numpy as np
import math
def read_hamiltonian(path):

    """
    Read hopping matrix element from the wannier90 output file:wannier90_hr.dat
    wan_num: number of wannier functions
    wsc_num: number of wigner-sitez cells
    wsc_count: variables related to degeneracy and hopping weight
    """

    with open(path,"r") as f:
        lines=f.readlines()

    wan_num=int(lines[1]); wsc_num=int(lines[2])
    ski_row_num=int(np.ceil(wsc_num/15.0))      # skip row numbers
    wsc_count=[]
    for i in range(ski_row_num):
        wsc_count.extend(list(map(int,lines[i+3].split())))

    wsc_tot=np.zeros((wan_num**2*wsc_num,3))
    tem_tot=np.zeros((wan_num**2*wsc_num,2))
    for i in range(wan_num**2*wsc_num):
        wsc_tot[i,:]=list(map(int,lines[3+ski_row_num+i].split()[:3]))
        tem_tot[i,:]=list(map(float,lines[3+ski_row_num+i].split()[5:]))

    wsc_idx=wsc_tot[0:-1:wan_num**2,:]   # the translational vector between wigner-sitez cells
    hop_mat=np.reshape(tem_tot[:,0]+1j*tem_tot[:,1],[wan_num,wan_num,wsc_num],order='F')
    return hop_mat,wsc_tot,wsc_count,wan_num,wsc_num

def write_hamiltonian(hop_mat,wsc_tot,wsc_count,wan_num,wsc_num,path):
    """
    Write hopping matrix element to the wannier90 output file:wannier90_hr.dat
    """
    with open(path,"w") as f:
        f.writelines('Hamiltonian data\n')
        f.writelines(str(wan_num)+'\n'); f.writelines(str(wsc_num)+'\n')
        ski_row_num=int(np.ceil(wsc_num/15.0))      # skip row numbers
        for i in range(ski_row_num-1):
            lists=wsc_count[i*15:(i+1)*15]
            line = '    '.join(str(j) for j in lists)
            f.writelines('    '+line+'\n')
        l=0
        for i in range(wsc_num):
            for j in range(wan_num):
                for k in range(wan_num):
                    line = '    '.join("%.0f" % (n,) for n in wsc_tot[l])
                    line2 ="   %.6f   %.6f" % (hop_mat[k][j][i].real,hop_mat[k][j][i].imag)
                    f.writelines('    '+line+line2+'\n')
                    l+=1
                    print(l)
    return 

def read_hamdata(path):
    """
    Read hopping matrix element from the FPLO +hamdata output file
    """
    with open(path,"r") as f:
        lines=f.readlines()
    translation_list =[]
    re_hopping = []
    im_hopping =[]        
    for i in range(len(lines)-1):
        #print('LINE = ', line)
        text = lines[i]
        next_text = lines[i+1]
        if 'nwan:' in text:
            wan_num = int(next_text)
            Tij = [[ [] for _ in range(wan_num)] for _ in range(wan_num)]
            Hij_re = [[ [] for _ in range(wan_num)] for _ in range(wan_num)]
            Hij_im = [[ [] for _ in range(wan_num)] for _ in range(wan_num)]
            
        if 'Tij, Hij:' in text:
            if '      ' in next_text:
                index = list(map(int, next_text.split( )))
                for j in range(2,10):
                    check = lines[i+j]
                    if 'end Tij, Hij:' in check:
                        maxindex = i+j
                for j in range(maxindex-i-3):
                    if 'E' in lines[j+i+2]:
                        line = list(map(float, lines[j+i+2].split( )))
                        translation = line[:3]
                        real_hop = line[3]
                        im_hop = line[4]
                    else:
                        translation = []
                        real_hop = [0]
                        im_hop = [0]
                    Tij[index[0]-1][index[1]-1].append(translation)
                    Hij_re[index[0]-1][index[1]-1].append(real_hop)
                    Hij_im[index[0]-1][index[1]-1].append(im_hop)

    return wan_num,Tij,Hij_re,Hij_im

def convert_hamdata(wan_num,Tij,Hij_re,Hij_im):

    """
    Convert the format of +hamdata to wannier_hr.dat
    """
    xdiv,ydiv,zdiv = 0,0,0
    for i in range(wan_num):
        for j in range(wan_num):
            for k in range(len(Tij[i][j])):
                if sum([abs(l) for l in Tij[i][j][k]]) != 0:
                    if Tij[i][j][k][0]-round(Tij[i][j][k][0]/4.340)<=0.1:
                        Tij[i][j][k][0] = round(Tij[i][j][k][0]/4.340)
                    else:
                        Tij[i][j][k][0] += 0.14469264133
                        Tij[i][j][k][0] = round(Tij[i][j][k][0]/4.340)
                    if Tij[i][j][k][1]-round(Tij[i][j][k][1]/5.012)<=0.1:
                        Tij[i][j][k][1] = round(Tij[i][j][k][1]/5.012)
                    else:
                        Tij[i][j][k][1] += 0.250615006
                        Tij[i][j][k][1] = round(Tij[i][j][k][1]/5.012)   
                    Tij[i][j][k][2] = round(Tij[i][j][k][2]/4.217)
    lists =[]
    for i in range(-2,3):
        for j in range(-2,3):
            for k in range(-1,2):
                for m in range(wan_num**2):
                    lists.append([i,j,k])
    wsc_tot = lists
    wsc_num = int(len(wsc_tot)/(wan_num**2))
    hop_mat = [[[ [] for _ in range(wsc_num)] for _ in range(wan_num)]for _ in range(wan_num)]
    for i in range(wan_num):
        for j in range(wan_num): 
            for k in range(wsc_num):
                for l in range(len(Tij[i][j])):
                    if Tij[i][j][l] == wsc_tot[k]:
                        hop_mat[i][j][k]=Hij_re[i][j][l]+1j*Hij_im[i][j][l]
                    else: 
                        hop_mat[i][j][k]= 0.+1j*0.
                if hop_mat[i][j][k] == []:
                    hop_mat[i][j][k] = 0.+1j*0.

    #print(hop_mat)
    #wsc_tot[i,:]=list(map(int,1)) 
    #print(wsc_tot[i,:])
    wsc_count=[2 for i in range(30)]
    for i in range(6):
        wsc_count += [1, 2, 2, 2]
                

    return hop_mat,wsc_tot,wsc_count,wan_num,wsc_num

#hop_mat,wsc_tot,wsc_count,wan_num,wsc_num=read_hamiltonian('wannier90_hr_FM.dat')
#print(wsc_tot)
wan_num,Tij,Hij_re,Hij_im = read_hamdata('+hamdata')
hop_mat,wsc_tot,wsc_count,wan_num,wsc_num = convert_hamdata(wan_num,Tij,Hij_re,Hij_im)
write_hamiltonian(hop_mat,wsc_tot,wsc_count,wan_num,wsc_num,'converted.dat')