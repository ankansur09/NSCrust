#!usr/bin/python/env

'''
%   ==========================================================================================================
%   Description:  A program to study the evolution of magnetic field through Hall effect and Ohmic dissipation
%              based on the work by Pablo Marchant et. al (2014)
%
%   Version: 2.0
%   Date: May 2020

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
from numba import jit, float32
import faster

PURE_TOROIDAL = False
PUREOHM = True
case=1

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

dr=(1.0-rmin)/(Nr-1)
dth=np.pi/(Nth-1)


#initilize the scalar functions
for i in range(1,Nr+1):
    r = (rmin+(i-1)*dr)
    for j in range(1,Nth-1):
        th = (j*dth)
        if PURE_TOROIDAL==False:
            A[i][j]= faster.A_init(r,case,th)
        B[i][j]= faster.B_init(r,case,th)


np.savetxt(out_folder+'initial_A.txt',A)
np.savetxt(out_folder+'initial_B.txt',B)

#BC at outer radius for B        
for j in range(Nth):
    B[Nr][j] = 0.0
    B[Nr+1][j] = -B[Nr-1][j]
    #A[1][j]=0
    #A[0][j] = A[2][j]
    
#BC at top and bottom for A and B (along the axis)
for i in range(Nr+2):
    B[i][0] = 0.0
    B[i][Nth-1] = 0.0
    if PURE_TOROIDAL==False:
        A[i][0] = 0.0
        A[i][Nth-1] = 0.0


hall_term_A, hall_rflux, hall_thflux, res_term_A, res_thflux,res_rflux, sines, cotans, boundary_factors1, boundary_factors2,\
    legendre_comb, legendre_comb_mid = faster.solve_repeated_values(rmin, Nr, Nth, dr, dth, thtd, l)


B = faster.solve_B_boundary(B, Nr, Nth)
A = faster.solve_A_boundary(A, a, Nr, Nth, l, legendre_comb_mid, legendre_comb, boundary_factors1, boundary_factors2)

tNum = 100
steps = 10
factor = 0.1
t = 0
dt = 1.0

faster.time_evolve(A, B, rmin, thtd, Nr, Nth, res_term_A, gsA, hall_term_A, res_rflux, hall_rflux, res_thflux, hall_thflux, dBr, dBth,  tNum, steps, factor, t, dt, dr, dth, cotans, sines)
