#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Jul 1st 2019

@author: camila
@editor: valentina
"""

from context import *

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