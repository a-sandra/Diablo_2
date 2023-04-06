import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys,os

matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

sys.path.append("/Users/pickle/Documents/GitHub/Diablo_2")
#sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")
from beam import beam
import madx_utils.madx_utils as mxu
import interactive_tool.build_beamline as seq
import optics_utils.transfer_matrices as opu
from dashboard.dashboard import MainWindow


print(os.getcwd())

lattice_file = "lattice_example/optics_71.7_MeV.twiss"
df_optics = mxu.get_twiss(lattice_file)
input_distribution = "distribution_example/CCL_Output_71.7.dat"

# test with the sequence file
## create the beam line
bl_with_input = seq.Beamline(lattice_file, input_distribution)

# change quadrupole strength
## before
#print(bl_with_input.dataframe_madx_sequence[bl_with_input.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])

##after
#list_quad = ["H01_EMQ_01", "H01_EMQ_02", "H01_EMQ_03"]
#k1 = [1.5, -8.5, 5.0]
#bl_with_input.update_magnet_parameter("QUADRUPOLE", list_quad, "K1L", k1)
#print(bl_with_input.dataframe_madx_sequence[bl_with_input.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])

#------ Dashboard ----------
app = QtWidgets.QApplication(sys.argv)
window = QtWidgets.QWidget()

window = MainWindow(bl_with_input)
#
window.show()
app.exec_()