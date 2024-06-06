#!usr/bin/python/env

'''
%   ==========================================================================================================
%   Description:  A program to study the evolution of magnetic field through Hall effect and Ohmic dissipation
%              based on the work by Pablo Marchant et. al (2014)
%
%   Version: 6.0
%   Date: May 2020
    Added additional alpha and beta evolution from mass accretion and ambipolar diffusion

%   Author: Ankan Sur
%   Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
%   ==========================================================================================================
'''


import numpy as np
from scipy.special import legendre
import sys, os
from optparse import OptionParser
from scipy import special
import time
import timeit
from scipy import interpolate

PURE_TOROIDAL = False
PUREOHM = False
ACCRETION = True
AMBIPOLAR = False
ohmic_test = False
case=6

start = timeit.default_timer()


#define the initial physical variables

def A_init(r,rmin,th):
    if case==1:
        l = 1
        n_leg = 1
        AA = -0.55882
        BB = -0.52004
        k = 7.03266
    if case==2:
        l = 1
        n_leg = 2
        AA = -0.72288
        BB = -0.20659
        k = 19.12793
    if case==3:
        l = 1
        n_leg = 2
        AA = 1.45340
        BB = 0.41536
        k = 19.12793
    if case==4:
        AA = 2.13404331
        BB = 4.1663394 
        u = np.cos(th)
        kk = 3
        return r*(((np.sin(kk*r)/(kk*r)-np.cos(kk*r))*AA)/(kk*r)+((-np.sin(kk*r)-np.cos(kk*r)/(kk*r))*BB)/(kk*r))*(1-u**2)
     
    if case==5:
        z0 = 0.041
        z = 1.0-r
        s0 = 1.0
        
        if z<z0:
            return (1.0-(z**2/z0**2))*np.sin(th)*np.sin(th)
        else:
            return 0.0 
            
    if case==6:
        x = r
        xmin = rmin
        return (5.0*xmin**3/x - 3*xmin**5/x + 3*x**4 - 5*x**2)*np.sin(th)*np.sin(th)/30.0

           

def B_init(r,rmin,th):
    if case==1:
        l = 1
        n_leg = 1
        CC = -0.55882
        DD = -0.52004
        k = 7.03266
    if case==2:
        l = 1
        n_leg = 2
        CC = -0.72288
        DD = -0.20659
        k = 19.12793
    if case==3:
        l = 1
        n_leg = 2
        CC = 3.55576
        DD = -38.32549
        k = 25.18557
    if case==4:
        CC = 2.13404331
        DD = 4.1663394 
        #AA = 6.51667597 
        #BB = 0.69665458
        u = np.cos(th)
        kk = 3
        return r*(((np.sin(kk*r)/(kk*r)-np.cos(kk*r))*CC)/(kk*r)+((-np.sin(kk*r)-np.cos(kk*r)/(kk*r))*DD)/(kk*r))*(1-u**2)
           
    if case==5:
        return 0
        
    if case==6:
        return 0
     



###  Function n returns electron density at a given radius and theta

def n(r,th):
    return ((1.115-r)/(0.115))**4
    #return 1.0
   
###  Function chi (1/(r^2 sin(th)*sin(th) * n(r,th))) at a given radius and theta

def chi(r,th):
    return 1.0/(pow(r*np.sin(th),2)*n(r,th))

### Function eta returns the value of resistivity at a given radius and theta

def eta(r,th):
    return n(r,th)**(-2.0/3.0)

#all repeated/fixed values

def solve_repeated_values():
    hall_term_A = np.zeros((Nr+2,Nth))
    hall_rflux = np.zeros((Nr+2,Nth))
    hall_thflux = np.zeros((Nr+2,Nth))
    res_term_A = np.zeros((Nr+2,Nth))
    res_rflux = np.zeros((Nr+2,Nth))
    res_thflux = np.zeros((Nr+2,Nth))
    sines = np.zeros(Nth)
    cotans = np.zeros(Nth)
    for j in range(1,Nth-1):
        th = j*dth
        sines[j] = np.sin(th)
        cotans[j] = np.cos(th)/np.sin(th) 
    for i in range(Nr+2):
    #for i in range(1,Nr+1):
        r=rmin+(i-1)*dr
	for j in range(Nth-1):
	    if j==0:
	        th == 0.0001
	    else:
	        th = j*dth
            if PUREOHM == False:
                hall_term_A[i][j]= sines[j]*chi(r,th)/4/dr/dth
                #print hall_term_A[i][j]
	        hall_rflux[i][j] = chi(r+dr/2,th)/8.0/dr/dth
	        hall_thflux[i][j]=-chi(r,th+dth/2)/8.0/dr/dth
	    res_term_A[i][j] = thtd*eta(r,th)
	    res_rflux[i][j]  = thtd*eta(r+dr/2,th)/np.sin(th)/dr/dr
	    res_thflux[i][j] = thtd*eta(r,th+dth/2)/r/r/np.sin(th+dth/2.0)/dth/dth

    boundary_factors1= np.zeros(l)
    boundary_factors2= np.zeros(l)
    legendre_comb = np.zeros((l,Nth))
    legendre_comb_mid = np.zeros((l,Nth))
    if PURE_TOROIDAL == False:
        for n in range(l):
	    boundary_factors1[n]=(n+1.0)/(n+2.0)*np.sqrt(np.pi*(2*n+3))/2*dth
	    boundary_factors2[n]=-2*dr*(n+1)*np.sqrt((2*n+3)/(4*np.pi))
	    for j in range(Nth-1):
	        th=j*dth
                Pl_np1 = legendre(n+1)
                Pl_n = legendre(n)
	        legendre_comb[n][j]=(np.cos(th)*Pl_np1(np.cos(th))- Pl_n(np.cos(th)))
	        legendre_comb_mid[n][j]=(np.cos(th+dth*0.5)*Pl_np1(np.cos(th+dth*0.5))-Pl_n(np.cos(th+dth*0.5)))/np.sin(th+dth*0.5)
	    
    #if PUREOHM == False:
    return hall_term_A, hall_rflux, hall_thflux, res_term_A, res_thflux,res_rflux, sines, cotans, boundary_factors1, boundary_factors2, \
           legendre_comb, legendre_comb_mid
    

def solve_A_boundary(A, a, legendre_comb_mid, legendre_comb, boundary_factors1, boundary_factors2, t):
    
    if ohmic_test==True:
        for j in range(Nth):
            A[0][j] = 0.0
            th = j*dth
            A[Nr+1][j] = np.sin(th)*np.sin(th)*np.exp(-9*t)
        for i in range(Nr+2):
        
            A[i][0] = 0.0
            A[i][Nth-1] = 0.0
            
            
    else:
        for n in range(l):
            a[n]=0
            for j in range(Nth-1):
                a[n]+=(A[Nr][j]+A[Nr][j+1])*legendre_comb_mid[n][j]*boundary_factors1[n]


        #Fix value of A at rmin = 0	
        ###for j in range(1,Nth-1):
        for j in range(Nth):
            ## this BC implies \alpha(rmin, theta) = 0.0 (equation A7)
	    A[1][j]=0
            # Grad Shafranov of A = 0 implies the below condition
            A[0][j]=-A[2][j]
            #### this implies delA/delr = 0 at r=R         
	    A[Nr+1][j]=A[Nr-1][j]
            #### the field is allowed to penetrate outside
	    for n in range(l):
	        A[Nr+1][j]+=boundary_factors2[n]*a[n]*legendre_comb[n][j]   
	        
	        
	for i in range(Nr+2):
	    A[i][0] = 0.0
            A[i][Nth-1] = 0.0
    
    return A


def solve_B_boundary(B, t):
    #Fix beta value just outside the star so the numerical radial  derivative at the surface corresponds
    #to solving it backwards (i.e. using only the point at the surface and the one just below). Also, Apply linear
    #interpolation of beta for the point just below the surface and just above the inner boundary (this last one
    #only in the case of zero boundary conditions at the center).
    if ohmic_test==True:
        for j in range(Nth):
            th = j*dth
            B[0][j] = 0.0
            B[Nr+1][j] = np.sin(th)*np.sin(th)*np.exp(-9*t)
        for i in range(Nr+2):
        
            B[i][0] = 0.0
            B[i][Nth-1] = 0.0
       
    else:
        for j in range(1,Nth-1):
        #for j in range(Nth+1):
	    B[1][j]=0
            #B[0][j] = 0
	    #B[2][j]=B[3][j]/2
            #B[1][j] = B[2][j]/2
	    #B[0][j]=2*dr/sc_factors[j]*(B[1][j]*(B[1][j+1]-B[1][j-1])/(2*dth*pow(rmin,2)*sines[j]))+B[2][j]
	    #B[0][j]=B[2][j]
	    B[Nr][j]=0
	    #B[Nr-1][j]=B[Nr-2][j]/2
	    B[Nr+1][j]=-B[Nr-1][j]
	    
        for i in range(Nr+2):
            B[i][0] = 0.0
            B[i][Nth-1] = 0.0
                     
    return B


#update A

def solve_new_A(A, B, res_term_A, gsA, hall_term_A, dt):
    for i in range(1,Nr+1):
    #for j in range(1,Nth-1):
	r=rmin+(i-1)*dr
	#r = rmin+i*dr
	#th=j*dth
	#Evolve poloidal field function at point
	                    
        Aaux[i][1:Nth-1]=A[i][1:Nth-1]+dt*res_term_A[i][1:Nth-1]*gsA[i][1:Nth-1]
        if PUREOHM == False:
            Aaux[i][1:Nth-1]+= dt*((B[i][2:Nth]-B[i][0:Nth-2])*(A[i+1][1:Nth-1]-A[i-1][1:Nth-1])-(B[i+1][1:Nth-1]-B[i-1][1:Nth-1])*(A[i][2:Nth]-A[i][0:Nth-2]))*hall_term_A[i][1:Nth-1]
		
    return Aaux

#update B

def solve_new_B(A, B, res_rflux, hall_rflux, gsA, res_thflux, hall_thflux, dt, dBr, dBth):
    
    for i in range(Nr+1):
        for j in range(Nth-1):
            if (j!=-1):
                
                dBr[i][j]=dt*res_rflux[i][j]*(B[i+1][j]-B[i][j])
                if PUREOHM == False:
                    dBr[i][j]+=dt*hall_rflux[i][j]*(B[i][j]+B[i+1][j])*(B[i][j+1]+B[i+1][j+1]-B[i][j-1]-B[i+1][j-1])
                    if PURE_TOROIDAL == False:
                        dBr[i][j]+=dt*hall_rflux[i][j]*(gsA[i][j]+gsA[i+1][j])*(A[i][j+1]+A[i+1][j+1]-A[i][j-1]-A[i+1][j-1])
            if (i!=-1):
                
                dBth[i][j]=dt*res_thflux[i][j]*(B[i][j+1]-B[i][j])
		if PUREOHM == False:
                    dBth[i][j]+=dt*hall_thflux[i][j]*(B[i][j]+B[i][j+1])*(B[i+1][j]+B[i+1][j+1]-B[i-1][j]-B[i-1][j+1])
		    if PURE_TOROIDAL == False:
                        dBth[i][j]+=dt*hall_thflux[i][j]*(gsA[i][j]+gsA[i][j+1])*(A[i+1][j]+A[i+1][j+1]-A[i-1][j]-A[i-1][j+1])
								
    return dBr,dBth
    
    
  
    
def solve_new_B_acc(A, B, Mdot, dt, rho):
    B_acc = np.zeros((Nr+2,Nth))
    for i in range(1,Nr+1):
        #for j in range(1,Nth-1):
        r=rmin+(i-1)*dr
	    #B_acc[i][j] = dt*Macc_factor*((B[i+1,j]-B[i-1,j])/2/dr/r**2/rho[i,j] - 2*B[i,j]/r**3/rho[i,j] - B[i,j]*(rho[i+1,j]-rho[i-1,j])/2/dr/rho[i,j]/rho[i,j]/r**2)
        B_acc[i,1:Nth-1] = dt*Macc_factor*((B[i+1,1:Nth-1]/(r+dr)**2/rho[i+1]) - (B[i-1,1:Nth-1]/(r-dr)**2/rho[i-1]))/2.0/dr
    
    return B_acc
  

def solve_new_A_acc(A, B, Mdot, dt, rho):
    '''
    Bint = np.zeros((Nr+2,Nth))
    for j in range(Nth):
        for i in range(Nr+1):
            th = j*dth
            rmin2 = rmin-dr
            r = rmin + i*dr
            B1 = B[:i+1,j]
            radius = np.linspace(rmin2,r,len(B1))
            Bint[i,j] = np.trapz(B1, x=radius)
    A_acc = np.zeros((Nr+2,Nth))
    for i in range(1,Nr+1):
        for j in range(Nth):
	    r=rmin+(i-1)*dr
	    A_acc[i][j] = dt*Macc_factor*(Bint[i,j])/4/np.pi/r**2/rho
    return A_acc
    '''
    A_acc = np.zeros((Nr+2,Nth))
    for i in range(1,Nr+1):
        #for j in range(Nth):
	r=rmin+(i-1)*dr
	A_acc[i][1:Nth-1] = dt*Macc_factor*(A[i+1,1:Nth-1]-A[i-1,1:Nth-1])/2/dr/r**2/rho[i]
	    
    return A_acc
    


############################################## end of all defined functions #######################################

parser = OptionParser()
  
parser.add_option('--Nr', type='int', default=100, help='number of radial points')
parser.add_option('--Nth', type='int', default=100, help='number of angular points')
parser.add_option('--rmin', type='float', default=0.8, help='location of crust-core interface as a multiple of R')
parser.add_option('--thtd', type='float', default=1.0, help='hall timescale divided by dynamica timescale')
parser.add_option('--l', type='int', default=1, help='multipole for outer boundary magnetic field')
parser.add_option('--out', type='string', dest='out_folder', default = 'output/', help='output folder')

(options, args) = parser.parse_args()

Nr = options.Nr
Nth = options.Nth
rmin = options.rmin
thtd = options.thtd
l = options.l
out_folder = options.out_folder

if not os.path.exists(out_folder):
  os.makedirs(out_folder)


#define physical quantities
A = np.zeros((Nr+2,Nth))
Aaux = np.zeros((Nr+2,Nth))
gsA = np.zeros((Nr+2,Nth))
a = np.zeros(l)

B = np.zeros((Nr+2,Nth))
dBr = np.zeros((Nr+2,Nth))
dBth = np.zeros((Nr+2,Nth))
#rho = np.ones((Nr+2,Nth))

dr=(1.0-rmin)/(Nr-1)
dth=np.pi/(Nth-1)

def interpolatef(x,y):
    f = interpolate.interp1d(x, y)
    xnew = np.linspace(8.775*1e5,max(x),Nr+2)
    ynew = f(xnew)
    return xnew,ynew


rad, den = np.genfromtxt('sly10.dat',usecols=(1,4),unpack=True)
ind = np.where(rad>8.0)
rad = rad[ind]
den = den[ind]
r_new, den_new = interpolatef(rad,den)
rho = np.array(den_new)/1e14
'''

rad, mass = np.genfromtxt('sly_m14.all',usecols=(2,6),unpack=True, skip_header=1)
ind = np.where(rad>7.2)
rad = rad[ind]*1e5
mass = mass[ind]*1.98e33
#print rad
density = 3*mass/rad**3/4/np.pi
r_new, den_new = interpolatef(rad,density)
rho = np.array(den_new)/1e14
'''

if ACCRETION == True:
    Mdot = 1.0
    #rho = 1
    Macc_factor = 1.0

rad = []
theta = []
#initilize the scalar functions
for i in range(Nr+2):
    r = (rmin+(i-1)*dr)
    for j in range(Nth):
        th = (j*dth)
        if PURE_TOROIDAL==False:
            A[i][j]= A_init(r,rmin,th)
        B[i][j]= B_init(r,rmin,th)
        theta.append(th)
        rad.append(r)




np.savetxt(out_folder+'initial_A.txt',A)
np.savetxt(out_folder+'initial_B.txt',B)

hall_term_A, hall_rflux, hall_thflux, res_term_A, res_thflux,res_rflux, sines, cotans, boundary_factors1, boundary_factors2,\
    legendre_comb, legendre_comb_mid = solve_repeated_values()


B = solve_B_boundary(B,0)
A = solve_A_boundary(A, a, legendre_comb_mid, legendre_comb, boundary_factors1, boundary_factors2, 0)

tNum = 10000000
steps = 100
time_final = 0.1
factor = 0.1
t = 0
dt = 1.0

theta = np.array(theta)
rad  = np.array(rad)


print_times = np.logspace(-6,-1,100)
pt = 0


################################### start time evolution ##################################################



for k in range(tNum+1):

    newdt=10000000
    localnewdt=10000000
    temp=0

    for i in range(1,Nr+1):
        r=rmin+(i-1)*dr
        gsA[i][1:Nth-1] = (A[i+1][1:Nth-1]+A[i-1][1:Nth-1]-2*A[i][1:Nth-1])/dr/dr+1/r/r*(A[i][2:Nth]+A[i][0:Nth-2]-2*A[i][1:Nth-1])/dth/dth - (np.cos(theta[1:Nth-1])/np.sin(theta[1:Nth-1]))*(A[i][2:Nth]-A[i][0:Nth-2])/2/dth/r**2
    
    
            
    #Obtain critical timestep for Hall drift
    for i in range(1,Nr+1):
        for j in range(1,Nth-1):
	    r=rmin+(i-1)*dr
	    th=j*dth
	    temp=dr/np.sqrt(pow(gsA[i][j]/(r*sines[j]),2)+pow((B[i+1][j]-B[i-1][j])/(2*dr*r*sines[j]),2)+pow((B[i][j+1]-B[i][j-1])/(2*dth*r*r*sines[j]),2))
            #print temp
	    if(temp<localnewdt):
                localnewdt=temp
			
        if(localnewdt<newdt):
            newdt=localnewdt

    #If Ohmic critical step is smaller, use that one
    if(newdt>dr*dr/thtd):
        newdt=dr*dr/thtd
    
    if ACCRETION == True:
    	if (newdt>np.abs(dr/Macc_factor)):
    	    newdt = np.abs(dr/Macc_factor)

    dt=factor*newdt
    if(dt<0.00000001):
        dt=0.00000001
     
    t+=dt
    #print "dt=", dt, "time=", t
    if t>time_final:
        print "finishing"
        break
    
    if PURE_TOROIDAL == False:
        Aaux = solve_new_A(A, B, res_term_A, gsA, hall_term_A, dt)
    
    dBr,dBth = solve_new_B(A, B, res_rflux, hall_rflux, gsA, res_thflux, hall_thflux, dt, dBr, dBth)
    
    if ACCRETION==True:
        B_acc = solve_new_B_acc(A, B, Mdot, dt, rho)
        A_acc = solve_new_A_acc(A, B, Mdot, dt, rho)
 
    for i in range(1,Nr+1):
        #for j in range(1,Nth-1):
            
        #if (np.isinf(dBr[i][1:Nth-2]) or np.isinf(dBth[i,1:Nth-2]) or np.isnan(dBr[i][1:Nth-2]) or np.isnan(dBth[i][1:Nth-2]) or np.isinf(Aaux[i][1:Nth-2]) or np.isnan(Aaux[i][1:Nth-2])):
        #    print "exiting code"
        #    quit()

        B[i][1:Nth-1] = B[i,1:Nth-1] +(dBr[i][1:Nth-1]-dBr[i-1][1:Nth-1]+dBth[i][1:Nth-1]-dBth[i][0:Nth-2])*sines[1:Nth-1]
        
        
        if ACCRETION==True:
            B[i][1:Nth-1]+= B_acc[i][1:Nth-1] 
                 
        if PURE_TOROIDAL == False:
            # Aaux have the A[i][j] term added, so no need to add now! Aaux has the dt term 
            A[i][1:Nth-1]=Aaux[i][1:Nth-1]
            if ACCRETION==True:
                A[i][1:Nth-1]+=A_acc[i][1:Nth-1]
        
               

    B = solve_B_boundary(B, t)

    if PURE_TOROIDAL == False:
        A = solve_A_boundary(A, a, legendre_comb_mid, legendre_comb, boundary_factors1, boundary_factors2, t)

   
    #if (k%steps==0):
        #np.savetxt(out_folder+'A_%d.txt'%k, A)
        #np.savetxt(out_folder+'B_%d.txt'%k, B)
        #np.savetxt(out_folder+'gsA_%d.txt'%k, gsA)
        
    if t>print_times[pt]:
        print "printing %d"%t
        np.savetxt(out_folder+'A_%d.txt'%pt, A)
        np.savetxt(out_folder+'B_%d.txt'%pt, B)
        pt = pt+1
        

    

#################### end of time evolution ################################

stop = timeit.default_timer()

print('Time: ', stop - start) 


    












    
    
