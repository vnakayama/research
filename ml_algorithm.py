"""
Computation of Cosmic Microwave Background Radiation (CMB) anisitropies simulated-maps 
with HEALpy, Astropy, Shogun, Scikit-Learn, Numpy and Pandas Python-libraries 
"""

import healpy as hp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import sympy as sym
#import sklearn as s
import sys
import time
#import seaborn as sns
#import camb
#import astropy
#import shogun

"""
@author: "Valentina Nakayama"
@data: "2019-09-18"
@email: "nakayama@ufrj.br"
"""

# Path to FITS Image files of four CMB missions, from the Plank Legacy Archive
observed_maps_dict = {"SMICA":"~/Documents/Research/Codes/Maps/COM_CMB_IQU-smica_2048_R3.00_full.fits",
            "COMMANDER":"~/Documents/Research/Codes/Maps/COM_CMB_IQU-commander_2048_R3.00_full.fits",
            "NILC":"~/Documents/Research/Codes/Maps/COM_CMB_IQU-nilc_2048_R3.00_full.fits",
            "SEVEM":"~/Documents/Research/Codes/Maps/COM_CMB_IQU-sevem_2048_R3.00_full.fits"
}

# Path to FITS Image files of four CMB missions , from the Plank Legacy Archive
observed_maps_masks_dict = {"SMICA":"",
                            "COMMANDER":"",
                            "NILC":"",
                            "SEVEM":""
}

# Define parameters for the given topology as [n_max,L,L_z,L/Lz]

topology_restraints = {"SLAB":[[4.0,4.0,1.0,4/1.0],
                               [4.0,4.0,1.2,4/1.2],
                               [4.0,4,1.4,4/1.4],
                               [4.0,4,1.6,4/1.6],
                               [4.0,4,1.8,4/1.8],
                               [4.0,4,2.0,4/1.2],
                               [4.0,4.0,4.0,4.0/4.0]],
                       "TORUS":[[4,1.2,1.2,4/1.2],
                               [4,1.4,1.4,4/1.4],
                               [4,1.6,1.6,4/1.6],
                               [4,1.8,1.8,4/1.8]]
} # don't use TORUS yet

# Restraints calculus

def norms(topology): # pass argument as string "topology" to match dict.keys()
    n_max = topology_restraints[topology][0]
    L = topology_restraints[topology][1]
    L_z = topology_restraints[topology][2]
    L_2 = topology_restraints[topology][3]
    norms = list()
    for n_x in np.arange(0,n_max+1,1):
        for n_y in np.arange(0,n_max+1,1):
            for n_z in range(1,n_max+1): # maybe change to np.arange
                N = (n_x)**2 + (n_y)**2 + (L_2*n_z)**2
                norms.append(N)
    norms = list(set(norms)) # remove duplicates with set() 
    norms.pop(0) #remove zero norm
    return sorted(norms)



def vectors(topology): 
    L_2 = topology_restraints[topology][3]
    vectors = {}
    for norm in norms(topology): #norm2 := norm*norm   
        vectores[np.sqrt(norm)] = list()
        n_x_count = np.sqrt(norm)
        for n_x in np.arange(0,n_x_count+2,1): # we do not include components at and bellow the nz=0 plane
            n_y2_max = norm - (n_x)**2
            if n_y2_max > 0.:
                n_y_count = np.sqrt(n_y2_max)
                for n_y in np.arange(0,n_y_count+2,1):
                    if norm - (n_x)**2 - (n_y)**2 > 0.:
                        L_2n_z = math.sqrt(norm - (nx)**2 - (ny)**2)
                        if L_2n_z !=0 and float(format(L2nz/L2, '0.8f')) == round(L_2n_z/L_2) and (n_x)**2 + (n_y)**2 + (L_2n_z)**2 == norm:                         
                            vectors[np.sqrt(norm)].append([n_x,n_y,L_2n_z])
    return vectors

# Simulation
def simulated_maps(topology,l_max,Nside):
    angles = list() # angle matrix
    for vector in vectors(topology).keys():
        angles.append(list(hp.vec2ang(np.array(vectors(topology)[vector][0]))))
    l, m  = hp.Alm.getlm(l_max)
    index = hp.Alm.getidx(l_max, l, m)
    # more calculus
    # mapa = hp.alm2map()
    # Cl_s = hp.anafast(mapa)
    # plt.plot()
    # plt.title("%s " % map)
    # plt.xtitle()
    # plt.ytitle()
    # plt.legend()
    # hp.mollview()

# Observed map import and analysis
def observed_maps(map):
    # apply mask observed_maps_masks_dict[map]
    observed_map = hp.read_map(observed_maps_dict[map])
    observed_cls = hp.anafast(observed_map)
    observed_alms = hp.map2alm(observed_map)
    # plt.plot()
    # plt.title("%s " % map)
    # plt.xtitle()
    # plt.ytitle()
    # plt.legend()
    # hp.mollview()

# Classification algorithm using the simulations as score

# Plot observed, simulated topology and best fit