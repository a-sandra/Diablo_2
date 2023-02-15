import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys,os

sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")
#from beam import beam
import madx_utils.madx_utils as mxu
import interactive_tool.sequence_2_python as seq
import optics_utils.transfer_matrices as opu

print(os.getcwd())

lattice_file = "lattice_example\optics_71.7_MeV.twiss"
df_optics = mxu.get_twiss(lattice_file)
#print(df_optics[:10])

# test with the sequence file
bl_with_input = seq.Beamline(lattice_file)

#