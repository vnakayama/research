#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:08:16 2019

@author: camila
"""

#%%

import healpy as hp
import numpy as np
from scipy import special
import time
from astropy.io import fits
import sys

#%%
#.....................................
## Reading the norms and angles:
#.....................................

params = (sys.argv)[1:] # -> nameFile, nmaps, nmax, Lmax, Lz

if len(params) == 5:
    nameFile = params[0]
    nmaps    = int(params[1])
    nmax     = float(params[2])
    Lmax     = float(params[3])
    Lz       = float(params[4])
    typeTop  = 'Slab'; folder = 'SlabTop'
    file_log = open("run_SWonly_"+typeTop+"_"+str(Lz)+".log","a")
    print(typeTop, 'topology:', nmaps, 'maps with', 'Lx=Ly=', Lmax, ', Lz=', Lz, ', nmax=', nmax)
    file_log.write('Slab topology: '+str(nmaps)+' maps with Lx=Ly= '+str(Lmax)+
                   ' , Lz= '+str(Lz)+', nmax= '+str(nmax))

if len(params) == 4: 
    nameFile = params[0]
    nmaps    = int(params[1])
    nmax     = float(params[2])
    Lmax     = float(params[3])
    Lz       = Lmax
    typeTop  = 'Torus'; folder = 'TorusTop'
    file_log = open("run_SWonly_"+typeTop+"_"+str(Lz)+".log","a")
    print(typeTop, 'topology:', nmaps, 'maps with', 'Lx=Ly=Lz=', Lmax, ', nmax=', nmax)
    file_log.write('Torus topology: '+str(nmaps)+' maps with Lx=Ly=Lz='+str(Lmax)+', nmax='+str(nmax))
    file_log.flush()

if len(params) < 4:    
    print('The input parameters should be, in this order: nameFile, nmaps, nmax, Lmax, Lz!')
    sys.exit()

hdu0 = fits.open(nameFile+'.fits')
norm_ang = hdu0[0].data 
print('Norms and angles - size:', len(norm_ang))

norm_ang = norm_ang[np.where(norm_ang[:,0] <= nmax)[0]]
print('Norms and angles - size:', len(norm_ang))
file_log.write('\n Norms and angles - size: %i'% len(norm_ang))
file_log.flush()


#%%
#.....................................
# Defining parameters:
#.....................................

# NOTE:
# The choice of a R value is important here, since different value of R lead to 
# differente values of the k and, therefore, different values of the transfer function.
# The A value on the other hand happens to be only a multiplicative value and can 
#corrected by the posterior nomalization of the simulations.


lmax = 20 # Maximum multipole
Nside = 16

l, m  = hp.Alm.getlm(lmax)
index = hp.Alm.getidx(lmax, l, m)

####
# Parameters from Planck 2018:
# pivot_scalar (Mpc^-1)
ks = 0.05
# scalar_amp          
As = 2.100549e-9
# scalar_spectral_index 
ns = 0.9660499
####

R = 13869.653929413116 # Mpc at zstar = 1089.9096 .. calculated by CAMB
#A = 1 # np.sqrt(2./np.pi)
k = 2.0*np.pi*norm_ang[:,0]/(Lmax*R) # norm of Fourier vector
Pk = ((2.*(np.pi**2.)*As)/(ks**(ns-1.))) * k**(ns-4.)
sqrt_Pspec = np.sqrt(Pk)
V = (Lmax*Lmax*Lz*(R**3.)) # Volume of the fundamental domain.
Cte = ((2.*np.pi)**3. / V )



#%%

#.....................................
# Finally a loop for each l,m pair:
#.....................................

kR = k * R

partial_time = time.time()

for nmp in range(1,nmaps+1):
    np.random.seed(seed=(nmp**2)*round(2.*Lz))
    rand_num1 = ( np.random.randn(len(k)) + np.random.randn(len(k))*1j )
    rand_num2 = ( np.random.randn(len(k)) + np.random.randn(len(k))*1j )
    rand_num3 = ( np.random.randn(len(k)) + np.random.randn(len(k))*1j )
    rand_num4 = ( np.random.randn(len(k)) + np.random.randn(len(k))*1j )
    alm = np.complex128(np.zeros(len(index)))
    for idx in index:
        if l[idx] > 1:
            jl_kR = special.spherical_jn(l[idx], kR)
            Ylm1  = special.sph_harm( np.float64(m[idx]), np.float64(l[idx]), norm_ang[:,2], norm_ang[:,1]) #( m, l, phi, theta)
            psi1 = (rand_num1 + ((-1.)**l[idx])*np.conjugate(rand_num1))
            psi2 = (rand_num2 + ((-1.)**l[idx])*np.conjugate(rand_num2))
            psi3 = (rand_num3 + ((-1.)**l[idx])*np.conjugate(rand_num3))
            psi4 = (rand_num4 + ((-1.)**l[idx])*np.conjugate(rand_num4))
            sum1 = ( ((psi1 + ((-1.)**m[idx])*psi3) * np.conjugate(Ylm1)) 
                  +  ((psi2 + ((-1.)**m[idx])*psi4) * Ylm1) )
            sum2 = sqrt_Pspec * jl_kR
            Total_sum = sum2 * sum1
            alm[idx] = ((1j)**l[idx] * Cte) * np.sum( Total_sum ) # para um par l,m
    if typeTop == 'Slab':
        hp.write_alm(folder+'/almSlab_SWonly_Ns16_nmax'+str(round(nmax))+'_LxLy'+str(round(Lmax))+
                     '_Lz'+str(Lz)+'_mc'+str(nmp)+'.fits', alm, overwrite=True)  
    if typeTop == 'Torus':
        hp.write_alm(folder+'/almTorus_SWonly_Ns16_nmax'+str(round(nmax))+'_L'+str(Lmax)+
                     '_mc'+str(nmp)+'.fits', alm, overwrite=True)    
#    print('Map:', nmp, 'Seed=', (nmp**2)*round(2.*Lz), 'in %f s' % (time.time() - partial_time))
    file_log.write('\n Map:'+str(nmp)+', Seed='+str((nmp**2)*round(2.*Lz))+', in %f s' 
                   % (time.time() - partial_time))
    file_log.flush()

#print('Loop finished in {time} s'.format(time = time.time() - partial_time))
file_log.write('\n Loop finished in %f s' % (time.time() - partial_time))
file_log.flush()

file_log.close()
