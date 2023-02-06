#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:38:05 2022

@author: pickle
"""

import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

font = {'size'   : 14}
plt.matplotlib.rc('font', **font)
import scipy.constants
E_0 = scipy.constants.physical_constants["proton mass energy equivalent in MeV"][0]*0.001 # (938.2720813, 'MeV', 5.8e-06) in GeV

sys.path.append("/Users/pickle/Documents/Diablo_2/")
from beam import beam
import madx_utils.madx_utils as mxu

import interactive_beam_envelop as ibe

distribution_file = "/Users/pickle/ADAM/L-HEBT/L-HEBT_Matching/13MKS_Phase_Offset/CCL_Output_71.7.dat"
lattice_file = "lattice_example/optics_71.7_MeV.twiss"

interactive_plot = ibe.Interactive(distribution_file, lattice_file, "H01")

interactive_plot.k1q11 = -10.79
interactive_plot.k1q12 = 4.585246875
interactive_plot.k1q13 = -6.76480708
interactive_plot.k1q14 = 8.216087706
interactive_plot.k1q15 = -5.623878856

interactive_plot.plotInteractive()


#interactive_plot.plotSimple()

print(interactive_plot.sigma)

print(interactive_plot.k1q11)
print(interactive_plot.k1q12)
print(interactive_plot.k1q13)
print(interactive_plot.k1q14)
print(interactive_plot.k1q15)