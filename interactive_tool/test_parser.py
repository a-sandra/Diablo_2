import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys,os

sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")
from beam import beam
import madx_utils.madx_utils as mxu
import interactive_tool.sequence_2_python as seq
import optics_utils.transfer_matrices as opu

print(os.getcwd())

lattice_file = "lattice_example\optics_71.7_MeV.twiss"
df_optics = mxu.get_twiss(lattice_file)
input_distribution = "distribution_example/CCL_Output_71.7.dat"

# test with the sequence file
#bl_with_input = seq.Beamline(lattice_file)
#print(bl_with_input.quadrupole_list)

#updated_HEBT_quadrupole = ["H01_EMQ_01", "H01_EMQ_02",  "H01_EMQ_04", "H01_EMQ_05"]
#updated_strength = np.array([5.3, -10.4, 3.0, -1.34,  -8.5])

#bl_with_input.update_magnet_parameter("QUADRUPOLE", updated_HEBT_quadrupole, "K1L",updated_strength)

#bm = beam.Beam(input_distribution)
#print(bm.distribution)

print(df_optics.iloc[2]["KEYWORD"], df_optics.iloc[2]["S"])