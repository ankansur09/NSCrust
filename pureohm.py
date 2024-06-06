#!usr/bin/python/env

'''
==========================================================================================================
Description:  Pure Ohmic Solution

Version: 1.0
Date: May 2020

Author: Ankan Sur
Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
==========================================================================================================
'''


import numpy as np
from scipy.special import legendre
import sys, os
from scipy import special
from matplotlib import pyplot as plt
import fnmatch

###############  define input params ######################
Nr = 60
Nth = 60

rmin = 0.75
thtd = 1
l = 1

t_start = 0.0
time_final = 0.02
tNum = 1000000
steps = 100
CFL = 0.1
t = t_start
dt = 1e-6


################## define functions ########################
case = 1
plotting = True
out_folder = 'Output3/'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

def A_init(r,rmin,th):
    if case==1:
        l=1
        n_leg=1
        AA = 2.13404331
        BB = 4.1663394 
        kk = 3
    Pl_n = legendre(n_leg)
    return r*(((np.sin(kk*r)/(kk*r)-np.cos(kk*r))*AA)/(kk*r)+((-np.sin(kk*r)-np.cos(kk*r)/(kk*r))*BB)/(kk*r))*np.sin(th)*np.sin(th)

def solve_A_boundary(A,t):
    for j in range(Nth):
        A[0][j] = 0.0
        A[Nr+1][j] = np.sin(th)*np.sin(th)*np.exp(-t*9)
    for i in range(Nr+2):
        A[i][0] = 0.0
        A[i][Nth-1] = 0.0
    return A
    
    
def load_files(steps):
    for t in range(startnt,endnt):
        a[t] = np.loadtxt(out_folder+'A_%d.txt'%(t*steps))
    return a
    
    


A = np.zeros((Nr+2,Nth))
gsA = np.zeros((Nr+2,Nth))

dr=(1.0-rmin)/(Nr+1)
dth=np.pi/(Nth-1)

#initilize the scalar functions
for i in range(1,Nr+1):
    r = (rmin+(i)*dr)
    for j in range(1,Nth-1):
        th = (j*dth)
        A[i][j]= A_init(r,rmin,th)
        
        
A = solve_A_boundary(A,0)  
np.savetxt(out_folder+'initial_A.txt',A)

###################### start time evolution ##############################


for k in range(tNum+1):

    t+=dt 
    
    for i in range(1,Nr+1):
        for j in range(1,Nth-1):
            r=rmin+(i)*dr
            th = j*dth
            gsA[i][j] = (A[i+1][j]+A[i-1][j]-2*A[i][j])/dr/dr+1/r/r*(A[i][j+1]+A[i][j-1]-2*A[i][j])/dth/dth - (np.cos(th)/np.sin(th))*(A[i][j+1]-A[i][j-1])/2/dth/r**2
    
    

    #dt = 1e-6
    

    for i in range(1,Nr+1):
        for j in range(1,Nth-1):
            
            A[i][j]=A[i][j]+dt*gsA[i][j]

    A = solve_A_boundary(A,t)
    #print output files
    if (k%steps==0):
        np.savetxt(out_folder+'A_%d.txt'%(k), A)

    
    if t>time_final:
        break
    print "dt=", dt, "time=", (t), "step=", k
    
 

################### end of time evolution ###########################


########################### Plotting #################################

radii, thetas = [],[]
for i in range(Nr):
    radii.append(i*dr + rmin)
for j in range(Nth):
    thetas.append(j*dth)

radii = np.array(radii)
thetas = np.array(thetas)


if plotting:
    startnt = 0
    endnt = len(fnmatch.filter(os.listdir(out_folder), '*.txt'))-1
    endtime = time_final
    nt = endnt-startnt
    sct = 1
    r_id = 5
    th_id = 5
    K = 3
    rad = radii[r_id]
    th = thetas[th_id]
    a = np.zeros((nt),object)
    alpha = load_files(steps)
    list_A = []
    for t in range(0,nt):
        AA = alpha[t]
        list_A.append(AA[r_id,th_id])   
    A_rth = np.array(list_A)
    time = np.linspace(0,endtime,nt)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    tau = 1.0/K**2
    alpha_init = np.loadtxt(out_folder+'initial_A.txt')
    alpha_init_rth = alpha_init[r_id,th_id]
    alpha_y = np.zeros(len(time))
    for t in range(len(time)):
        alpha_y[t] = alpha_init_rth*np.exp(-time[t]/tau)
    
    f = plt.figure(figsize=(14,10))
    ax = f.add_subplot(111)
    ax.plot(time,alpha_y, label = r'$\alpha$ analy', lw=3,c='r', alpha=0.5)
    ax.scatter(time[0::sct],A_rth[0::sct],s=50, facecolors='none', edgecolors='r', label="sim")
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'time',size=25)
    ax.set_ylabel(r'$\alpha (t)$',size=25)
    ax.set_title('rad/R='+ str(round(rad,2)) + ",   theta (deg)=" + str(round(th*180/np.pi,2)))
    ax.legend(fontsize=20)
    ax.set_xlim(0,endtime)
    plt.savefig(out_folder+'alpha_time.png', bbox_inches='tight')    
    plt.close(f)