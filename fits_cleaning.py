"""
Computation of Cosmic Microwave Background Radiation (CMB) anisitropies simulated-maps 
with HEALpy, Astropy, Shogun, Scikit-Learn, Numpy and Pandas Python-libraries

Cleaning FITS Images
"""

import healpy as hp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sklearn as sk
import sys
import time
#import seaborn as sns
#import camb
#import astropy
#import shogun

"""
@author: "Valentina Nakayama"
@data: "2019-10-1"
@email: "nakayama@ufrj.br"
"""

# Path to FITS Image files of four CMB missions, from the Plank Legacy Archive
observed_maps_dict = {"SMICA":"/Users/gui/Documents/Research/Codes/Maps/COM_CMB_IQU-smica_2048_R3.00_full.fits",
            "COMMANDER":"/Users/gui/Documents/Research/Codes/Maps/COM_CMB_IQU-commander_2048_R3.00_full.fits",
            "NILC":"/Users/gui/Documents/Research/Codes/Maps/COM_CMB_IQU-nilc_2048_R3.00_full.fits",
            "SEVEM":"/Users/gui/Documents/Research/Codes/Maps/COM_CMB_IQU-sevem_2048_R3.00_full.fits"
}

# Path to FITS Image files of four CMB missions confidence masks, from the Plank Legacy Archive
# Find the pipeline specific one
observed_maps_masks_dict = {"SMICA":"",
                            "COMMANDER":"",
                            "NILC":"",
                            "SEVEM":""
}

def apply_mask(mission):
    observed_map = fits.open(observed_maps_dict[mission], memap=True)
    observed_map_mask = fits.open(observed_maps_masks_dict[mission], memap=True)
    # ...

