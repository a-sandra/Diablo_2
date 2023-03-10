import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys,os

sys.path.append("/Users/pickle/Documents/GitHub/Diablo_2")
#sys.path.append(r"\Python_Packages\Diablo_2")
from beam import beam
import madx_utils.madx_utils as mxu
import interactive_tool.build_beamline as seq
import optics_utils.transfer_matrices as opu

print(os.getcwd())

lattice_file = "lattice_example/optics_71.7_MeV.twiss"
df_optics = mxu.get_twiss(lattice_file)
input_distribution = "distribution_example/CCL_Output_71.7.dat"

# test with the sequence file
bl_with_input = seq.Beamline(lattice_file, input_distribution)
print(bl_with_input.df_beam_size_along_s)

plt.plot(bl_with_input.df_beam_size_along_s["S"], bl_with_input.df_beam_size_along_s["X"], "o")
plt.plot(df_optics["S"], df_optics["SIGMA_X"], "o")
plt.show()

#updated_HEBT_quadrupole = ["H01_EMQ_01", "H01_EMQ_02",  "H01_EMQ_04", "H01_EMQ_05"]
#updated_strength = np.array([5.3, -10.4, 3.0, -1.34,  -8.5])

#bl_with_input.update_magnet_parameter("QUADRUPOLE", updated_HEBT_quadrupole, "K1L",updated_strength)
