#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 09:56:17 2019

@author: camila
"""

#%%

import healpy as hp
import math
import numpy as np
from scipy import special
#import sph_harm_byHand
import time

#import range_multipole_Map as rm
import matplotlib.pyplot as plt
from pylab import cm
cmap = cm.jet #stores old healpy colloring scheme (cmap='jet')
cmap.set_under('w')

#%%
#.....................................
# Definitions:
#.....................................

def norms(nmax,L,lz): #generates a list with the norms up to a given cutoff
    L2 = L/lz # or: float("%.3f" % (L/lz))
    l = []
    for nx in np.arange(0,nmax+1,1):
        for ny in np.arange(0,nmax+1,1):
            for nz in range(1,nmax+1):
                norm = (nx)**2 + (ny)**2 + (L2*nz)**2 # <-- nx, ny, nz correspond to (L/lx)*nx, (L/ly)*ny, (L/lz)*nz
#                print(norm,nx,ny,nz)
                l.append(norm) 
    l = list(set(l)) #we remove duplicates with set() 
    l.pop(0) #remove zero norm
    return sorted(l)

###############################

def vector_n(nmax,L,lz):
    L2 = L/lz
    vec_n = {}
#    print('norms=',len(norms(nmax,L,lz)))
    for norm in norms(nmax,L,lz): #norm2 := norm*norm
#        print('   >>norm=',norm)        
        vec_n[math.sqrt(norm)] = list()
        nx_count = math.sqrt(norm)
        for nx in np.arange(0,nx_count+2,1): #we do not include components at and bellow the nz=0 plane
            ny2_max = norm - (nx)**2
            if ny2_max > 0.:
                ny_count = math.sqrt(ny2_max)
                for ny in np.arange(0,ny_count+2,1):
                    if norm - (nx)**2 - (ny)**2 > 0.:
                        L2nz = math.sqrt(norm - (nx)**2 - (ny)**2)
#                        print(L2nz/L2,'-', round(L2nz/L2), '/', (nx)**2 + (ny)**2 + (L2nz)**2,'-',norm)
                        if L2nz !=0 and float(format(L2nz/L2, '0.8f')) == round(L2nz/L2) and (nx)**2 + (ny)**2 + (L2nz)**2 == norm:
#                            print('***')                            
                            vec_n[math.sqrt(norm)].append([nx,ny,L2nz])
#                            print(nx,ny,L2nz)
    return vec_n

###############################

def vectors_to_angles(n_dict):
    """
    :param n_dict: { norm=sqrt(((L/lx)*nx)^2 + ((L/ly)*ny)^2 + ((L/lz)*nz)^2), 
                            [ (L/lx)*nx, (L/ly)*ny, (L/lz)*nz ] }
    :return: the conversion of ((L/lx)*nx, (L/ly)*ny, (L/lz)*nz), cartesian coordinates, 
                            to (theta_k,phi_k), spherical coordinates. 
    """
    new_n = dict()
    for k in n_dict.keys(): # <- for in 'norm'
        angles = []
        for vector in n_dict[k]:
            # calculate the azimuthal angle
            if vector[0] != 0: # <- (L/lx)*nx
                phi = math.atan(float(vector[1])/float(vector[0])) 
                # <- arctan( (L/ly)*ny / (L/lx)*nx) = tan( (Lx*ny) / (Ly*nx) )
            else:
                phi = math.pi/2  # tan^-1(inf) = pi/2
                if vector[1] == 0:  # tan^-1(0/0) = 0.0
                    phi = 0.0
            # calculate the polar angle
            if vector[2] != 0: # <- (L/lz)*nz
                # sqrt returns float, float/int = float
                theta = math.atan((math.sqrt(vector[0]**2 + vector[1]**2))/vector[2])
            else:
                theta = math.pi/2  # tan^-1(inf) = pi/2
            angles.append([theta, phi]) # <- theta_k, phi_k
        new_n.update({k: angles})
    return new_n

#%%
#.....................................
# Calculation of norms, vectors and angles:
#.....................................

nmax = 4
Lmax = 4    # Size of the fundamental domain: Lx = Ly = Lmax and Lz.
Lz   = 0.1

partial_time = time.time()

normvec_dict = vector_n(nmax,Lmax,Lz)   # <- {sqrt(((L/lx)*nx)^2 + ((L/ly)*ny)^2 + ((L/lz)*nz)^2), [(L/lx)*nx, (L/ly)*ny, (L/lz)*nz]}
print(len(norms(nmax,Lmax,Lz)))
normang_dict = vectors_to_angles(normvec_dict) #initialize norms and angles
print(len(normang_dict))


#.....................................
# Tranforming the dict into a matrix:
#.....................................

n = len(normang_dict.keys())

print(list(normang_dict.keys())[0])
#print(normang_dict[list(normang_dict.keys())[0]][1][0])

norm_ang = []
for i in range(n):
  for j in range(len(normang_dict[list(normang_dict.keys())[i]])):
    norm_ang.append((list(normang_dict.keys())[i], 
                     normang_dict[list(normang_dict.keys())[i]][j][0], 
                     normang_dict[list(normang_dict.keys())[i]][j][1]))

norm_ang = np.array(norm_ang)
print(len(norm_ang))

print('Finished calculating norms, vectors and angles in {time} s'.format(time = time.time() - partial_time))

# theta = norm_ang[:,1] -> latitude
# phi   = norm_ang[:,2] -> longitude

#%%

#.....................................
# Defining parameters:
#.....................................

lmax = 35 #int(2.*64.-1.) # Maximum multipole
Nside = 64

l, m  = hp.Alm.getlm(lmax)
index = hp.Alm.getidx(lmax, l, m)

A = 1 # np.sqrt(2./np.pi)
R = 1 # Mpc ==>> R = R_H = \chi_rec = 45.7 x 10^9 anos luz (1 pc = 3.26 anos-luz) => R = 14.0184 Mpc
k = 2.0*np.pi*norm_ang[:,0]/(Lmax*R) # norm of Fourier vector
sqrt_Pspec = np.sqrt(A*(k**(-3.0)))
Cte = (1./3.) * np.sqrt(2./np.pi)

#####


kR = k * R

alm = np.complex128(np.zeros(len(index)))

np.random.seed(seed=654321)#(idx**2)*round(2.*Lz))
rand_num1 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num2 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num3 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
rand_num4 = ( np.random.randn(len(norm_ang)) + np.random.randn(len(norm_ang))*1j )
        
for idx in index:
    if l[idx] > 1:
        jl_kR = special.spherical_jn(l[idx], kR)
        Ylm1  = special.sph_harm( np.float64(m[idx]), np.float64(l[idx]), norm_ang[:,2], norm_ang[:,1]) #( m, l, phi, theta)
        #print(m[idx], '/', l[idx], '/', '/', Ylm1[0])
        psi1 = (rand_num1 + ((-1.)**l[idx])*np.conjugate(rand_num1))
        psi2 = (rand_num2 + ((-1.)**l[idx])*np.conjugate(rand_num2))
        psi3 = (rand_num3 + ((-1.)**l[idx])*np.conjugate(rand_num3))
        psi4 = (rand_num4 + ((-1.)**l[idx])*np.conjugate(rand_num4))
        sum1 = ( ((psi1 + ((-1.)**m[idx])*psi3) * np.conjugate(Ylm1)) 
              +  ((psi2 + ((-1.)**m[idx])*psi4) * Ylm1) )
        sum2 = sqrt_Pspec * jl_kR
        Total_sum = sum2 * sum1
        alm[idx] = ((1j)**l[idx] * Cte) * np.sum( Total_sum ) # para um par l,m

mapa = hp.alm2map(alm, Nside)
#hp.mollview(mapa, title='Slab, Lz=0.1, Lx=Ly=4, nmax=4', cmap=cmap, rot=[0.,210.,70.])
#plt.savefig('Slab_SWonly_seed654321_Lz0p1_LxLy4_nmax4.png', dpi=100)
hp.mollview(mapa, title='Torus, L = 0.1, nmax=4', cmap=cmap, rot=[0.,210.,70.])
plt.savefig('Torus_SWonly_seed654321_L0p1_nmax4.png', dpi=100)


