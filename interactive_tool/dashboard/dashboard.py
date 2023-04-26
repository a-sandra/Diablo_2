import sys, copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

from PyQt5.QtWidgets import QMainWindow, QApplication, QDockWidget, QWidget, QGridLayout, QSlider, QLabel, QDoubleSpinBox
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

matplotlib.rcParams.update({'font.size': 9})

class MainWindow(QMainWindow):
    def __init__(self, tranferline,*args, **kwargs):
        QMainWindow.__init__(self, *args, **kwargs)
        self.beamline = tranferline
        self.set_h01emq001(0)
        self.set_h01emq002(0)
        self.set_h01emq003(0)
        self.set_h01emq004(0)
        self.set_h01emq005(0)

        dock = QDockWidget ("Values")
        self.addDockWidget (Qt.RightDockWidgetArea, dock)
        sliders = QWidget ()
        sliders_grid = QGridLayout (sliders)

        def add_slider(col):
            sld = QSlider(Qt.Vertical, sliders)
            sld.setPageStep(0.5)
            sld.setMaximum(50.0)
            sld.setMinimum(-50.0)
            sld.setFocusPolicy(Qt.NoFocus)
            sliders_grid.addWidget (sld, 1, col)
            return sld

        def add_spinbox(foo, col, slider):
            spinbox = QDoubleSpinBox()
            spinbox.setMaximum(50.0)
            spinbox.setMinimum(-50.0)
            spinbox.valueChanged[float].connect(foo) #When the slider's value has changed
            spinbox.valueChanged.connect(self.plot)
            spinbox.valueChanged[float].connect(slider.setValue)
            sliders_grid.addWidget (spinbox, 0, col)
            slider.sliderReleased.connect(lambda:spinbox.setValue(slider.value()))


        sld1 = add_slider (col = 0)
        sld2 = add_slider (col = 1)
        sld3 = add_slider (col = 2)
        sld4 = add_slider (col = 3)
        sld5 = add_slider (col = 4)
        add_spinbox(foo = self.set_h01emq001, col = 0, slider = sld1)
        add_spinbox(foo = self.set_h01emq002, col = 1, slider = sld2)
        add_spinbox(foo = self.set_h01emq003, col = 2, slider = sld3)
        add_spinbox(foo = self.set_h01emq004, col = 3, slider = sld4)
        add_spinbox(foo = self.set_h01emq005, col = 4, slider = sld5)

        dock.setWidget(sliders)
        self.plot()
    
#    def convert_to_int(self, val):


    def set_h01emq001(self, val):
        self.k1_h01emq001 = val
        
    def set_h01emq002(self, val):
        self.k1_h01emq002 = val

    def set_h01emq003(self, val):
        self.k1_h01emq003 = val
        
    def set_h01emq004(self, val):
        self.k1_h01emq004 = val

    def set_h01emq005(self, val):
        self.k1_h01emq005 = val

    def plot(self):
        print("bonjour")
        print(self.beamline.dataframe_madx_sequence[self.beamline.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])
        list_quad = ["H01_EMQ_01", "H01_EMQ_02", "H01_EMQ_03", "H01_EMQ_04", "H01_EMQ_05"]
        k1l = [self.k1_h01emq001, self.k1_h01emq002, self.k1_h01emq003, self.k1_h01emq004, self.k1_h01emq005]
        #print(self.beamline.df_beam_size_along_s["X"])
        #print(k1l)
        self.beamline.update_magnet_parameter("QUADRUPOLE", list_quad, "K1L", k1l)
        print(self.beamline.dataframe_madx_sequence[self.beamline.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])

        ##tracking
        self.beamline.data = []
        #self.beamline.df_beam_size_along_s = []
        for index in np.arange(len(self.beamline.dataframe_madx_sequence)):
            self.beamline.check_element_type(self.beamline.dataframe_madx_sequence.iloc[index])
        self.beamline.main_dataframe_sequence = pd.DataFrame(self.beamline.data, columns=["NAME", "KEYWORD", "S","L", "K1","TRANSFER_MATRIX"])
        self.beamline.tm_as_function_s_array = np.asarray(self.beamline.build_beamline_transfer_matrix(), dtype=object)
        self.beamline.beam_distribution_at_every_element = self.beamline.cumulative_tracking()

        print(self.beamline.df_beam_size_along_s["X"])
        print(self.beamline.main_dataframe_sequence)

        s = self.beamline.dataframe_madx_sequence["S"]-self.beamline.dataframe_madx_sequence["L"]/2
        self.beamline.aperture_beamline()
        self.beamline.build_dataframe_element()
        self.maximum_s = max(self.beamline.dataframe_madx_sequence["S"])

        list_emq = [Rectangle(( self.beamline.df_emq["S"][i]-self.beamline.df_emq["L"][i]/2, -0.2/2), self.beamline.df_emq["L"][i], 0.2,color ='skyblue') for i in self.beamline.df_emq.index]
        list_stm = [Rectangle(( self.beamline.df_stm["S"][i]-self.beamline.df_stm["L"][i]/2, -0.2/2), self.beamline.df_stm["L"][i], 0.2,color ='limegreen') for i in self.beamline.df_stm.index]
        list_bpm = [Rectangle(( self.beamline.df_bpm["S"][i]-self.beamline.df_bpm["L"][i]/2, -0.2/2), self.beamline.df_bpm["L"][i], 0.2,color ='orange') for i in self.beamline.df_bpm.index]
        list_bend = [Rectangle(( self.beamline.df_bend["S"][i]-self.beamline.df_bend["L"][i]/2, -0.2/2), self.beamline.df_bend["L"][i], 0.2,color ='crimson') for i in self.beamline.df_bend.index]
        list_fsm = [Rectangle(( self.beamline.df_fsm["S"][i]-self.beamline.df_fsm["L"][i]/2, -0.2/2), self.beamline.df_fsm["L"][i], 0.2,color ='purple') for i in self.beamline.df_fsm.index]

        self.fi, ((self.ax1,self.ax3),(self.ax2, self.ax4))= plt.subplots(2, 2, figsize = (20,15), gridspec_kw={'height_ratios':[1,4]})
        for i in list_emq:
            nq=copy.copy(i)
            self.ax1.add_patch(nq)
        for i in list_stm:
            ns=copy.copy(i)
            self.ax1.add_patch(ns)
        for i in list_bpm:
            nb=copy.copy(i)
            self.ax1.add_patch(nb)
        for i in list_bend:
            nbend=copy.copy(i)
            self.ax1.add_patch(nbend)
        for i in list_fsm:
            nfsm=copy.copy(i)
            self.ax1.add_patch(nfsm)

        self.ax1.plot([0, self.maximum_s],[0,0], "k", linewidth=0.5)
        self.ax1.set_ylim(-0.2,0.6)
        self.ax1.set_xlim(0,6)
        self.ax1.grid()
        self.ax1.axes.yaxis.set_visible(False)
        self.ax1.axes.xaxis.set_visible(False)
        self.ax1.spines['top'].set_visible(False)
        self.ax1.spines['bottom'].set_visible(False)
        self.ax1.spines['left'].set_visible(False)
        self.ax1.spines['right'].set_visible(False)

        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], 3.0*self.beamline.df_beam_size_along_s["X"]*1000, "ob")
        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], -3.0*self.beamline.df_beam_size_along_s["X"]*1000, "ob")

        self.ax2.step(s, self.beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax2.step(s, -self.beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")

        self.ax2.step(s, self.beamline.aperture_y, where = 'pre', color = 'r', label = "_nolabel_")
        self.ax2.step(s, -self.beamline.aperture_y, where = 'pre', color = 'r', label = "_nolabel_")
        self.ax2.set_xlabel("s(m)")
        self.ax2.set_ylabel("Beam evelop (mm)")
        self.ax2.set_ylim(-30, 30)

        self.ax3.hist(self.beamline.bm.distribution["X(mm)"], bins=50)
        self.ax3.set_xlabel("X(mm)")
        self.ax3.set_title("initial distribution -  OPT/FSM in the future")

        self.ax4.scatter(self.beamline.bm.distribution["X(mm)"], self.beamline.bm.distribution["Y(mm)"], color = 'purple', s=0.5)
        self.ax4.set_xlabel("X(mm)")
        self.ax4.set_ylabel("Y(mm)")
 

        self.canvas = FigureCanvas(self.fi)
        self.setCentralWidget(self.canvas)
        self.canvas.draw()