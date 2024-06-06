#!usr/bin/python/env

'''
==========================================================================================================
Description:  Plotting the field components

Version: 1.0
Date: April 2020

Author: Ankan Sur
Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
==========================================================================================================
'''


import numpy as np
from matplotlib import pyplot as plt
import sys, os
from scipy.interpolate import griddata
from scipy import special
from scipy.special import legendre
import math

startnt = 0
endnt = 999
stept = 100
#endtime = 0.00411324786325
#endtime = 0.00822238658777
endtime = 1.0

sct = 25
#endtime = 0.00411119329389


eta = 1.0
nt = endnt-startnt
case=2

a = np.zeros((nt),object)
b = np.zeros((nt),object)

lNum = 1
rmin = 0.75
Nr = 36
Nth = 36
dr=(1.0-rmin)/(Nr+1)
dth=np.pi/(Nth+1)
radii, thetas = [],[]
time = np.linspace(0,endtime,nt)



#out_folder = 'Plots/fieldlines_chuck/'
#out_folder = 'Plots/test2/'
out_folder = 'Output/'
#input_folder = 'output_chuck/'
input_folder = 'Output/'
#input_folder = 'test2/'

if not os.path.exists(out_folder):
  os.makedirs(out_folder)

folder_numbers = []
for i in range(1000):
    folder_numbers.append(100000+i*100+1) 
    
    
def load_files():
    for t in range(startnt,endnt):
    #for t in range(len(folder_numbers)):
        #a[t] = np.loadtxt(input_folder+'A_%d.txt'%(folder_numbers[t]))
        #print t*stept
        a[t] = np.loadtxt(input_folder+'A_%d.txt'%(t*stept+stept))
        b[t] = np.loadtxt(input_folder+'B_%d.txt'%(t*stept+stept))
        #a[t] = np.fromfile(input_folder+'A_%d.bin'%(t*10),'u1')
        #a[t] = aa.reshape(Nr+2,Nth)
        #b[t] = np.fromfile(input_folder+'B_%d.bin'%(t*10),'u1')
        #b[t] = bb.reshape(Nr+2,Nth)
    return a,b


def plotter(Bpol, t, x, z):
    f = plt.figure(figsize=(12,14))
    ax = f.add_subplot(111)
    Bmag = ax.contourf(x, z, Bpol, 50, cmap='hot', extend='both')
    cbar = plt.colorbar(Bmag,fraction=0.046, pad=0.02)
    cbar.set_label(r"B$_{\rm tor}$ / B$_{\rm s}$", size=25)
    ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'x / R$_{\star}$',size=25)
    ax.set_ylabel(r'z / R$_{\star}$',size=25)
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=15)
    cbar.ax.tick_params(labelsize=15) 
    plt.savefig(out_folder+'Bpol_%d.png'%t, bbox_inches='tight')
    f.close()


def spherical_to_cartesian2D(r, theta):
    
    #Given 1D arrays for r and theta, the function makes a spherical (r,theta)
    #grid and then transforms it to cartesian coordinates. It outputs two 2D
    #arrays, one for x and one for z.
    
    theta_matrix, radius_matrix = np.meshgrid(theta, r)
    x = radius_matrix*np.sin(theta_matrix)
    z = radius_matrix*np.cos(theta_matrix)
    return x, z

def spherical_to_cartesian_bfield_2D(N_r, N_theta, theta, B_r, B_theta):
    """
    Given a vector field (B_r,B_theta) in spherical coordinates, the function
    outputs a field in cartesian coordinates.
    """

    B_x = np.zeros((N_r, N_theta))
    B_z = np.zeros((N_r, N_theta))

    for i in range(0, N_theta):
        B_x[:,i] = B_r[:,i]*np.sin(theta[i]) + B_theta[:,i]*np.cos(theta[i])
        B_z[:,i] = B_r[:,i]*np.cos(theta[i]) - B_theta[:,i]*np.sin(theta[i])
    return B_x, B_z

def fieldlines(r, theta, B_r, B_th, t, magB, timenow):
    x, z = spherical_to_cartesian2D(r, theta)
    theta0 = np.linspace(0.0,np.pi,Nth)
    x_r = np.cos(theta0)
    y_r = np.sin(theta0)
    w = 1
    Z_l, X_l = np.mgrid[-w:w:100j, 0:w:100j]
    px_l = x.flatten()
    pz_l = z.flatten()
    px_r = -x.flatten()
    pz_r = z.flatten()
    f = plt.figure(figsize=(12,14))
    ax = f.add_subplot(111)
    B_x_l, B_z_l = spherical_to_cartesian_bfield_2D(Nr, Nth, theta, B_r, B_th)
    pu_l = B_x_l.flatten()
    pv_l = B_z_l.flatten()
    gu_l = griddata(zip(px_l,pz_l), pu_l, (X_l,Z_l))
    gv_l = griddata(zip(px_l,pz_l), pv_l, (X_l,Z_l))
    ax.streamplot(X_l,Z_l,gu_l,gv_l, color='green', linewidth=1.5)
    Bmag = ax.contourf(x, z, magB, 50, cmap='hot', extend='both')
    cbar = plt.colorbar(Bmag,fraction=0.046, pad=0.02)
    cbar.set_label(r"B$_{\rm tor}$ / B$_{\rm 0}$", size=25)
    ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'x / R$_{\star}$',size=25)
    ax.set_ylabel(r'z / R$_{\star}$',size=25)
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=15)
    cbar.ax.tick_params(labelsize=15) 
    ax.set_xlim(0,1)
    ax.set_ylim(-1,1)
    #loc_arr_x = np.array([0.0,0.2,0.4,rmin,1.0])
    #loc_arr_y = np.array([-1.0,-rmin,-0.4,-0.2,0.0,0.2,0.4,rmin,1.0])
    #plt.xticks(loc_arr_x,[0.0,0.2,0.4,0.9,1.0],size=20)
    #plt.yticks(loc_arr_y,[-1.0,-0.9,-0.4,-0.2,0.0,0.2,0.4,0.9,1.0],size=20)
    #print timenow*13.0
    ax.set_title(r'$t$ ='+ str(np.round((timenow*13000.0),3)) +' Thyr',size=25)
    #ax.set_title(r't = %d'%(t),size=25)
    plt.savefig(out_folder+'fig%d.png'%t, bbox_inches='tight')
    plt.close(f)


def fieldlines2(A,B):
    x, z = spherical_to_cartesian2D(radii, thetas)
    B_r = np.zeros((Nr,Nth))
    B_th = np.zeros((Nr,Nth))
    B_phi = np.zeros((Nr,Nth))
    for i in range(Nr-1):
        r = i*dr+ rmin
        for j in range(Nth-1):
            th = j*dth + dth/2.0 
            #Solve each component
            #r component
            B_r[i][j] = 1.0/r/r/np.sin(th)*(A[i][j+1]-A[i][j]+A[i+1][j+1]-A[i+1][j])/2/dth
            #th component
            B_th[i][j] = -1.0/r/np.sin(th)*(A[i+1][j]-A[i][j]+A[i+1][j+1]-A[i][j+1])/2/dr
            #phi component        
            B_phi[i][j] = 1.0/r/np.sin(th)*(B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4

    Bpol = np.sqrt(B_r**2+B_th**2)

    B_r_up = np.zeros((Nr,Nth/2))
    B_th_up = np.zeros((Nr,Nth/2))
    B_r_down = np.zeros((Nr,Nth/2))
    B_th_down = np.zeros((Nr,Nth/2))
    for i in range(Nr-1):
        for j in range(Nth/2-1):
            B_r_up[i,j] = B_r[i,j]
            B_th_up[i,j] = B_th[i,j]
        
    for i in range(Nr-1):
        for j in range(Nth/2-1):
            B_r_down[i,j] = B_r[i,j+Nth/2]
            B_th_down[i,j] = B_th[i,j+Nth/2]



    XAvalues_up=np.zeros((Nr,Nth/2))
    YAvalues_up=np.zeros((Nr,Nth/2))
    theta_up = np.zeros(Nth/2)
    theta_down = np.zeros(Nth/2)
    for i in range(Nr):
        #use twice the value of the radial size of the crust to help visualization
        #r=0
        r= i*dr+rmin
        for j in range(Nth/2):
            th=j*dth
            theta_up[j] = th
            XAvalues_up[i][j]=r*np.sin(th)
            YAvalues_up[i][j]=r*np.cos(th)
        
    XAvalues_down=np.zeros((Nr,Nth/2))
    YAvalues_down=np.zeros((Nr,Nth/2))
    for i in range(Nr):
        #use twice the value of the radial size of the crust to help visualization
        #r=0
        r= i*dr+rmin
        for j in range(Nth/2):
            th=(j+Nth/2)*dth 
            theta_down[j] = th
            XAvalues_down[i][j]=r*np.sin(th)
            YAvalues_down[i][j]=r*np.cos(th)

    px_l_up = XAvalues_up.flatten()
    pz_l_up = YAvalues_up.flatten()
    B_x_l_up, B_z_l_up = spherical_to_cartesian_bfield_2D(Nr, Nth/2, theta_up, B_r_up, B_th_up)
    pu_l_up = B_x_l_up.flatten()
    pv_l_up = B_z_l_up.flatten()
    w = 1
    Z_l_up, X_l_up = np.mgrid[-0.1:1:100j, 0:1:100j]
    gu_l_up = griddata(zip(px_l_up,pz_l_up), pu_l_up, (X_l_up,Z_l_up))
    gv_l_up = griddata(zip(px_l_up,pz_l_up), pv_l_up, (X_l_up,Z_l_up))
    px_l_down = XAvalues_down.flatten()
    pz_l_down = YAvalues_down.flatten()
    B_x_l_down, B_z_l_down = spherical_to_cartesian_bfield_2D(Nr, Nth/2, theta_down, B_r_down, B_th_down)
    pu_l_down = B_x_l_down.flatten()
    pv_l_down = B_z_l_down.flatten()
    #w = 1
    Z_l_down, X_l_down = np.mgrid[-1:0.1:100j, 0:1:100j]
    #print px_l_down
    gu_l_down = griddata(zip(px_l_down,pz_l_down), pu_l_down, (X_l_down,Z_l_down))
    gv_l_down = griddata(zip(px_l_down,pz_l_down), pv_l_down, (X_l_down,Z_l_down))
    f = plt.figure(figsize=(12,14))
    ax = f.add_subplot(111)
    Bmag = ax.contourf(x,z, Bpol, 50, cmap='jet', extend='both', alpha=0.7)
    cbar = plt.colorbar(Bmag,fraction=0.046, pad=0.02)
    cbar.set_label(r"B$_{\rm pol}$ / B$_{\rm 0}$", size=25)
    ax.streamplot(X_l_up,Z_l_up,gu_l_up,gv_l_up,color='black',linewidth=1.0)
    ax.streamplot(X_l_down,Z_l_down,gu_l_down,gv_l_down,color='black',linewidth=1.0)
    ax.set_xlim(0,1)
    ax.set_ylim(-1,1)
    ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'x / R$_{\star}$',size=25)
    ax.set_ylabel(r'z / R$_{\star}$',size=25)
    #ticklabs = cbar.ax.get_yticklabels()
    plt.savefig(out_folder+'fig%d.png'%t, bbox_inches='tight')
    plt.close(f)
    #plt.show()

#def analytical(rad,th,A_rth,B_rth,r_id,th_id):
def analytical(rad,th,A_rth, B_rth, r_id,th_id):
    '''
    if case==1:
        k_a = 7.03266
        k_b = 12.67071
    if case==2:
        k_a = 19.12793
        k_b = 25.18557
    Pl_n = legendre(1)
    #alpha_ana = rad*(Aa*special.spherical_jn(l, k*rad, derivative=False)+Bb*special.spherical_yn(l, k*rad, derivative=False))*Pl_n(np.cos(th))*np.sin(th)
    '''
    l = 1
    k_a = 3
    k_b = 3
    tau_a = 1.0/k_a**2/eta
    tau_b = 1.0/k_b**2/eta
    alpha_init = np.loadtxt(input_folder+'initial_A.txt')
    alpha_init_rth = alpha_init[r_id,th_id]
    beta_init = np.loadtxt(input_folder+'initial_B.txt')
    beta_init_rth = alpha_init[r_id,th_id]
    alpha_y = np.zeros(len(time))
    beta_y = np.zeros(len(time))
    for t in range(len(time)):
        alpha_y[t] = alpha_init_rth*np.exp(-time[t]/tau_a)
        #alpha_y[t] = alpha_ana*np.exp(-time[t]/tau)
        beta_y[t] = beta_init_rth*np.exp(-time[t]/tau_b)
    
    f = plt.figure(figsize=(14,10))
    ax = f.add_subplot(111)
    error_A = np.abs((alpha_y-A_rth)/alpha_y)*100
    error_B = np.abs((beta_y-B_rth)/beta_y)*100
    ax.plot(time,alpha_y, label = r'$\alpha$ analytical', lw=3,c='r', alpha=0.5)
    #ax.plot(time,error_A, label = r'$\alpha$ analy', lw=3,c='r', alpha=0.5)
    ax.plot(time,beta_y, label = r'$\beta$ analytical',lw=3, c='b', alpha=0.5)
    ax.scatter(time[5::sct],A_rth[5::sct],s=80, marker = 'd', facecolors='none', edgecolors='r', label="sim")
    ax.scatter(time[0::sct],B_rth[0::sct],s=80, facecolors='none', edgecolors='b', label="sim")
    #ax.set_aspect('equal',adjustable='box')
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'time',size=25)
    ax.set_ylabel(r'$\alpha (t), \beta(t)$',size=25)
    #ax.set_title('rad/R='+ str(round(rad,2)) + ",   theta (deg)=" + str(round(th*180/np.pi,2)))
    ax.text(0.15, 0.15, 'r/R='+ str(round(rad,2)) + r",  $\theta$ =" + str(round(th*180/np.pi,2)) + r"$^{\circ}$", fontsize=30)
    ax.legend(fontsize=20)
    #ax.set_xlim(0,endtime)
    ax.set_xlim(0,endtime)
    #ax.set_ylim(0,0.5)
    #plt.savefig(out_folder+'alpha_time.png', bbox_inches='tight')    
    plt.savefig(out_folder+'Ohm_AB.pdf', bbox_inches='tight')    
    plt.close(f)
   



def plot_energy(Bpol_e, Btor_e):
    Bpol_e = np.array(Bpol_e)
    Btor_e = np.array(Btor_e)
    f = plt.figure(figsize=(14,10))
    ax = f.add_subplot(111)
    ax.plot(time,Bpol_e/(Bpol_e+Btor_e),label='Bpol', lw=3,c='r', alpha=0.9)
    ax.plot(time,Btor_e/(Bpol_e+Btor_e),label='Btor', lw=3,c='b', alpha=0.9)  
    ax.tick_params(axis="x", labelsize=25) 
    ax.tick_params(axis="y", labelsize=25)
    ax.set_xlabel(r'time',size=25)
    ax.set_ylabel(r'$energy (t)$',size=25)
    ax.legend(fontsize=20)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    #ax.set_xlim(0,0.002)
    #plt.savefig(out_folder+'alpha_time.png', bbox_inches='tight')    
    plt.savefig(out_folder+'energy.png', bbox_inches='tight')    
    plt.close(f)
    

def stream_function(r, theta, B_r, B_theta, order=1):
    """
    Given a spherical (r,theta) grid which can be nonuniform, and a
    divergenceless axisymmetric vector field B(r, theta), this function returns
    a scalar function psi(r, theta)called the Stokes stream function, such that
    B = rot(A), where A is vector potential given by
    A = psi(r,theta)/r*sin(theta) * e_phi. The contours of psi are the field
    lines of the original vector field B.

    Routine originally written for ASH data by ???
    Modified for PLUTO by S. Matt (2 Dec., 2011)
    Translated to Python by F. Bartolic(August, 2016)
    """

    N_r = len(r)
    N_theta = len(theta)

    #calculate derivatives of the stream function
    dpsi_dr = np.zeros((N_r,N_theta))
    dpsi_dtheta = np.zeros((N_r,N_theta))

    for i in range(0, N_theta):
        dpsi_dr[:, i] = -r*np.sin(theta[i])*B_theta[:, i]
        dpsi_dtheta[:, i] = r*r*np.sin(theta[i])*B_r[:,i]

    # start at integrating at pole at inner radius, do all along pole, then do
    # for each theta
    psi = np.zeros((N_r,N_theta))
    psi_2 = np.zeros((N_r, N_theta))
    dtheta = np.zeros(N_theta)
    dr = np.zeros(N_r)
    if order >= 0:
        dr[1:] = r[1:] - r[:-1]
        dtheta[1:] = theta[1:] - theta[:-1]
        psi[1:, 0] = psi[:-1, 0] + dpsi_dr[1:, 0]*dr[1:]

        for i in range(1, N_theta):
             psi[:, i] = psi[:, i-1] + dpsi_dtheta[:, i]*dtheta[i]

    if order <= 0:
        dr[:-1] = r[:-1] - r[1:]
        dtheta[:-1] = theta[:-1] - theta[1:]

        for i in range(N_r-2, -1, -1):
            psi_2[i, N_theta - 1] = psi_2[i + 1, N_theta - 1] +\
                                dpsi_dr[i, N_theta - 1]*dr[i]
        for i in range(N_theta-2, -1, -1):
            psi_2[:, i] = psi_2[:, i + 1] + dpsi_dtheta[:, i]*dtheta[i]

        if order < 0:
            return psi_2
        else:
            psi = 0.5*(psi + psi_2)  # Avg of the two
    return psi
    
alpha,beta = load_files()
#alpha = load_files()


for i in range(Nr):
    radii.append(i*dr + rmin)
for j in range(Nth):
    thetas.append(j*dth)

radii = np.array(radii)
thetas = np.array(thetas)


x, z = spherical_to_cartesian2D(radii, thetas)

#Create array for multipole coefficients
multipoles=np.zeros((lNum));
#store P_l^1(cos(theta)) at the surface, used to solve multipoles. I have an array evaluated in midpoints to perform
#the integrations, and another one solve on mesh points to perform alpha evaluation
plone=np.zeros((lNum,Nth));
plone2=np.zeros((lNum,Nth));
for j in range(Nth):
    alp=special.lpmn(1, lNum+1, math.cos((j+0.5)*dth))
    alp2=special.lpmn(1, lNum+1, math.cos(j*dth))
    for l in range(lNum):
        plone[l][j]=alp[0][1][l+1]
        plone2[l][j]=alp2[0][1][l+1]

#x,z = spherical_to_cartesian2D(radii,thetas)
maxB = []
Bpol_e = []
Btor_e = []
start = 0
intv = 2
clevel = np.zeros((Nr,Nth))


def time_evol():

    for t in range(nt):
        Bpe = 0
        Bte = 0
        time_now = time[t]
        #print time_now
    
        A = alpha[t]
        B = beta[t]
        
        A = A[1:-1]
        B = B[1:-1]

        B_r = np.zeros((Nr,Nth))
        B_th = np.zeros((Nr,Nth))
        B_phi = np.zeros((Nr,Nth))
        for i in range(Nr-1):
            r = i*dr+ rmin + dr/2.0
            #r = radii[i]
            #r = i*dr + dr/2.0 + rmin
            for j in range(Nth-1):
                th = j*dth + dth/2.0 
                #th = j*dth+dth/2.0
            #Solve each component
            #r component
                B_r[i][j] = 1.0/r/r/np.sin(th)*(A[i][j+1]-A[i][j]+A[i+1][j+1]-A[i+1][j])/2/dth
            #th component
                B_th[i][j] = -1.0/r/np.sin(th)*(A[i+1][j]-A[i][j]+A[i+1][j+1]-A[i][j+1])/2/dr
                Bpe += (B_r[i,j]**2+B_th[i,j]**2)
                #phi component        
                B_phi[i][j] = 1.0/r/np.sin(th)*(B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4
                #B_phi[i][j] = B[i][j] 
                Bte += (B_phi[i,j]**2)
        Bpol_e.append(Bpe)
        Btor_e.append(Bte)
        Bpol = np.sqrt(B_r**2+B_th**2)
        Btor = np.sqrt(B_phi**2)
        #B_x, B_z = spherical_to_cartesian_bfield_2D(Nr, Nth, thetas, B_r,B_th)
        Psi = stream_function(radii, thetas, B_r, B_th)
        f = plt.figure(figsize=(14,10))
        ax = f.add_subplot(111)
        Bmag = ax.contourf(x,z, Bpol, 100, cmap='jet', extend='both', alpha=0.7)
        cbar = plt.colorbar(Bmag,fraction=0.046, pad=0.02)
        Psimin = np.amin(Psi)
        #print Psimin
        Bfieldl = ax.contour(x, z, Psi, extend='both', linewidths=2, origin='lower', levels=np.linspace(Psimin,np.amax(Psi),25),colors='black',linestyles='solid') 
        #cbar = plt.colorbar(Bfieldl,fraction=0.046, pad=0.02)
        #ax.quiver(x[start::intv,::intv],z[start::intv,::intv],B_x[start::intv,::intv],B_z[start::intv,::intv],color='black',angles='xy')
        ax.set_xlim(0,1)
        ax.set_ylim(-1,1)
        ax.set_aspect('equal',adjustable='box')
        ax.tick_params(axis="x", labelsize=25) 
        ax.tick_params(axis="y", labelsize=25)
        ax.set_xlabel(r'x / R$_{\star}$',size=25)
        ax.set_ylabel(r'z / R$_{\star}$',size=25)
        plt.savefig(out_folder+'field_new%d.png'%t, bbox_inches='tight')
        plt.close(f)
        
        #fieldlines(radii, thetas, B_r, B_th, t, Bpol, time_now)
        
    return Bpol_e, Btor_e

def store_alpha_beta(r_id, th_id):
    list_A,list_B = [],[]
    for t in range(0,nt):
        AA = alpha[t]
        BB = beta[t]
        list_A.append(AA[r_id,th_id])
        list_B.append(BB[r_id,th_id])

    return np.array(list_A), np.array(list_B)
    #return np.array(list_A)
    
#Bpol_e, Btor_e = time_evol()
#plot_energy(Bpol_e, Btor_e)
r_id = 10
th_id = 10
# these strores alpha, beta at r, th for all times
A_rth, B_rth = store_alpha_beta(r_id, th_id)
#A_rth = store_alpha_beta(r_id, th_id)
analytical(radii[r_id],thetas[th_id], A_rth, B_rth, r_id, th_id)
#analytical(radii[r_id],thetas[th_id], A_rth, B_rth, r_id, th_id)
#for t in range(5):
 #   A = alpha[t]
 #   B = beta[t]
    #fieldlines2(A,B)
    



