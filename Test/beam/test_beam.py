import sys, os
import scipy

sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")

#sys.path
from beam import beam

# No input beam distribution
bm = beam.Beam()
print(bm.distribution, bm.E_ptc)

# with a beam distribution
path_to_distribution = r'C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2\Test\distributions\CCL_Output_99.0.dat'
input_distribution = r'..\distributions\CCL_Output_99.0.dat'
bm_distribution = beam.Beam(path_to_distribution)

print(os.getcwd())
print(os.listdir())

print(bm_distribution.distribution)

