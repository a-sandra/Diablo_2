import sys, copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

from scipy import stats
import random

from PyQt5.QtWidgets import QMainWindow, QApplication, QDockWidget, QWidget, QGridLayout, QSlider, QLabel, QDoubleSpinBox
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

matplotlib.rcParams.update({'font.size': 9})

class MainWindow(QMainWindow):
    def __init__(self, tranferline,*args, **kwargs):
        QMainWindow.__init__(self, *args, **kwargs)
        self.beamline = tranferline

        #self.beamline.get_quadrupole_dataframe()[["NAME", "K1"]]
        #for index_quadrupole in self.beamline.get_quadrupole_dataframe()[["NAME", "K1"]]
        #self.beamline.get_quadrupole_dataframe()

        self.set_h01emq001(self.beamline.quadrupole_df.iloc[0]["K1"]*self.beamline.quadrupole_df.iloc[0]["L"])
        self.set_h01emq002(self.beamline.quadrupole_df.iloc[1]["K1"]*self.beamline.quadrupole_df.iloc[1]["L"])
        self.set_h01emq003(self.beamline.quadrupole_df.iloc[2]["K1"]*self.beamline.quadrupole_df.iloc[2]["L"])
        self.set_h01emq004(self.beamline.quadrupole_df.iloc[3]["K1"]*self.beamline.quadrupole_df.iloc[3]["L"])
        self.set_h01emq005(self.beamline.quadrupole_df.iloc[4]["K1"]*self.beamline.quadrupole_df.iloc[4]["L"])
        self.set_h02emq001(self.beamline.quadrupole_df.iloc[5]["K1"]*self.beamline.quadrupole_df.iloc[5]["L"])
        self.set_h03emq001(self.beamline.quadrupole_df.iloc[6]["K1"]*self.beamline.quadrupole_df.iloc[6]["L"])
        self.set_h03emq002(self.beamline.quadrupole_df.iloc[7]["K1"]*self.beamline.quadrupole_df.iloc[7]["L"])
        self.set_h03emq003(self.beamline.quadrupole_df.iloc[8]["K1"]*self.beamline.quadrupole_df.iloc[8]["L"])
        self.set_h03emq004(self.beamline.quadrupole_df.iloc[9]["K1"]*self.beamline.quadrupole_df.iloc[9]["L"])
        self.set_h03emq005(self.beamline.quadrupole_df.iloc[10]["K1"]*self.beamline.quadrupole_df.iloc[10]["L"])

        dock = QDockWidget ("Values")
        self.addDockWidget (Qt.BottomDockWidgetArea, dock)
        sliders = QWidget ()
        sliders_grid = QGridLayout (sliders)

        def add_slider(col):
            sld = QSlider(Qt.Horizontal, sliders)
            sld.setPageStep(0.5)
            sld.setMaximum(50.0)
            sld.setMinimum(-50.0)
            sld.setFocusPolicy(Qt.NoFocus)
            sliders_grid.addWidget (sld, 1, col)
            return sld

        def add_spinbox(foo, col, slider, name, unit):
            spinbox = QDoubleSpinBox()
            spinbox.setMaximum(50.0)
            spinbox.setMinimum(-50.0)
            spinbox.valueChanged[float].connect(foo) #When the slider's value has changed
            spinbox.valueChanged.connect(self.plot)
            spinbox.valueChanged[float].connect(slider.setValue)
            spinbox.setPrefix(name)
            spinbox.setSuffix(unit)
            sliders_grid.addWidget (spinbox, 0, col)
            slider.sliderReleased.connect(lambda:spinbox.setValue(slider.value()))

        sld1 = add_slider (col = 0)
        sld2 = add_slider (col = 1)
        sld3 = add_slider (col = 2)
        sld4 = add_slider (col = 3)
        sld5 = add_slider (col = 4)

        sld6 = add_slider (col = 5)
        sld7 = add_slider (col = 6)
        sld8 = add_slider (col = 7)
        sld9 = add_slider (col = 8)
        sld10 = add_slider (col = 9)
        sld11 = add_slider (col = 10)
        add_spinbox(foo = self.set_h01emq001, col = 0, slider = sld1, name = "H01-EMQ-01: ", unit = " m^-2")
        add_spinbox(foo = self.set_h01emq002, col = 1, slider = sld2, name = "H01-EMQ-02: ", unit = " m^-2")
        add_spinbox(foo = self.set_h01emq003, col = 2, slider = sld3, name = "H01-EMQ-03: ", unit = " m^-2")
        add_spinbox(foo = self.set_h01emq004, col = 3, slider = sld4, name = "H01-EMQ-04: ", unit = " m^-2")
        add_spinbox(foo = self.set_h01emq005, col = 4, slider = sld5, name = "H01-EMQ-05: ", unit = " m^-2")
        #add_spinbox(foo = self.set_h02emq001, col = 5, slider = sld6, name = "H02-EMQ-01: ", unit = " m^-2")
        add_spinbox(foo = self.set_h03emq001, col = 6, slider = sld7, name = "H03-EMQ-01: ", unit = " m^-2")
        add_spinbox(foo = self.set_h03emq002, col = 7, slider = sld8, name = "H03-EMQ-02: ", unit = " m^-2")
        add_spinbox(foo = self.set_h03emq003, col = 8, slider = sld9, name = "H03-EMQ-03: ", unit = " m^-2")
        add_spinbox(foo = self.set_h03emq004, col = 9, slider = sld10, name = "H03-EMQ-04: ", unit = " m^-2")
        add_spinbox(foo = self.set_h03emq005, col = 10, slider = sld11, name = "H03-EMQ-05: ", unit = " m^-2")

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
    
    def set_h02emq001(self, val):
        self.k1_h02emq001 = val

    def set_h03emq001(self, val):
        self.k1_h03emq001 = val
    
    def set_h03emq002(self, val):
        self.k1_h03emq002 = val

    def set_h03emq003(self, val):
        self.k1_h03emq003 = val

    def set_h03emq004(self, val):
        self.k1_h03emq004 = val

    def set_h03emq005(self, val):
        self.k1_h03emq005 = val

    def plot(self):
        #print(self.beamline.dataframe_madx_sequence[self.beamline.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])
        list_quad = ["H01_EMQ_01", "H01_EMQ_02", "H01_EMQ_03", "H01_EMQ_04", "H01_EMQ_05"]
        k1l = [self.k1_h01emq001, self.k1_h01emq002, self.k1_h01emq003, self.k1_h01emq004, self.k1_h01emq005]
        #print(self.beamline.df_beam_size_along_s["X"])
        #print(k1l)
        self.beamline.update_magnet_parameter("QUADRUPOLE", list_quad, "K1L", k1l)
        #print(self.beamline.dataframe_madx_sequence[self.beamline.dataframe_madx_sequence["KEYWORD"]=="QUADRUPOLE"])

        ##tracking
        self.beamline.data = []
        #self.beamline.df_beam_size_along_s = []
        for index in np.arange(len(self.beamline.dataframe_madx_sequence)):
            self.beamline.check_element_type(self.beamline.dataframe_madx_sequence.iloc[index])
        self.beamline.main_dataframe_sequence = pd.DataFrame(self.beamline.data, columns=["NAME", "KEYWORD", "S","L", "K1","TRANSFER_MATRIX"])
        self.beamline.tm_as_function_s_array = np.asarray(self.beamline.build_beamline_transfer_matrix(), dtype=object)
        self.beamline.beam_distribution_at_every_element = self.beamline.cumulative_tracking()

        #print(self.beamline.df_beam_size_along_s["X"])
        #print(self.beamline.main_dataframe_sequence)

        s = self.beamline.dataframe_madx_sequence["S"]-self.beamline.dataframe_madx_sequence["L"]/2
        self.beamline.aperture_beamline()
        self.beamline.build_dataframe_element()
        self.maximum_s = max(self.beamline.dataframe_madx_sequence["S"])

        list_emq = [Rectangle(( self.beamline.df_emq["S"][i]-self.beamline.df_emq["L"][i]/2, -0.2/2), self.beamline.df_emq["L"][i], 0.2,color ='skyblue') for i in self.beamline.df_emq.index]
        list_stm = [Rectangle(( self.beamline.df_stm["S"][i]-self.beamline.df_stm["L"][i]/2, -0.2/2), self.beamline.df_stm["L"][i], 0.2,color ='limegreen') for i in self.beamline.df_stm.index]
        list_bpm = [Rectangle(( self.beamline.df_bpm["S"][i]-self.beamline.df_bpm["L"][i]/2, -0.2/2), self.beamline.df_bpm["L"][i], 0.2,color ='orange') for i in self.beamline.df_bpm.index]
        list_bend = [Rectangle(( self.beamline.df_bend["S"][i]-self.beamline.df_bend["L"][i]/2, -0.2/2), self.beamline.df_bend["L"][i], 0.2,color ='crimson') for i in self.beamline.df_bend.index]
        list_fsm = [Rectangle(( self.beamline.df_fsm["S"][i]-self.beamline.df_fsm["L"][i]/2, -0.2/2), self.beamline.df_fsm["L"][i], 0.2,color ='purple') for i in self.beamline.df_fsm.index]

        self.fi = plt.figure(figsize = (20,20)) #, tight_layout = True
        gs = gridspec.GridSpec(3,4, width_ratios=[1,1,1,3],height_ratios = [3,0.5,3], hspace = 0.3)

        self.ax3 = self.fi.add_subplot(gs[0, 0])
        self.ax4 = self.fi.add_subplot(gs[0, 2])
        self.ax5 = self.fi.add_subplot(gs[0, 1])
        self.ax6 = self.fi.add_subplot(gs[0, 3])
        self.ax1 = self.fi.add_subplot(gs[1, :])
        self.ax2 = self.fi.add_subplot(gs[2, :])

##-----------------------Layout of the line-----------------------##
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

##-----------------------Beam envelop-----------------------##

        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], 3.0*self.beamline.df_beam_size_along_s["X"]*1000, "--ob")
        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], -3.0*self.beamline.df_beam_size_along_s["X"]*1000, "--ob")
        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], 3.0*self.beamline.df_beam_size_along_s["Y"]*1000, "--og")
        self.ax2.plot(self.beamline.df_beam_size_along_s["S"], -3.0*self.beamline.df_beam_size_along_s["Y"]*1000, "--og")

        self.ax2.plot(self.beamline.dataframe_madx_sequence["S"], 3.0*self.beamline.dataframe_madx_sequence["SIGMA_X"]*1000, "--y")
        self.ax2.plot(self.beamline.dataframe_madx_sequence["S"], -3.0*self.beamline.dataframe_madx_sequence["SIGMA_X"]*1000, "--y")
        self.ax2.plot(self.beamline.dataframe_madx_sequence["S"], 3.0*self.beamline.dataframe_madx_sequence["SIGMA_Y"]*1000, "--k")
        self.ax2.plot(self.beamline.dataframe_madx_sequence["S"], -3.0*self.beamline.dataframe_madx_sequence["SIGMA_Y"]*1000, "--k")

        self.ax2.step(s, self.beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax2.step(s, -self.beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax2.step(s, self.beamline.aperture_y, where = 'pre', color = 'r', label = "_nolabel_")
        self.ax2.step(s, -self.beamline.aperture_y, where = 'pre', color = 'r', label = "_nolabel_")
        self.ax2.set_xlabel("s(m)")
        self.ax2.set_ylabel("Beam evelop (mm)")
        self.ax2.set_ylim(-30, 30)

##-----------------------Histogram-----------------------##
        # gaussian fit for beam distribution
        x_mean = np.mean(self.beamline.distribution_fsm_h01[:,0]*1000)
        x_rms_size = np.std(self.beamline.distribution_fsm_h01[:,0]*1000)
        horizontal_gaussian_fit = stats.norm.pdf(self.beamline.distribution_fsm_h01[:,0]*1000, x_mean, x_rms_size)

        y_mean = np.mean(self.beamline.distribution_fsm_h01[:,2]*1000)
        y_rms_size = np.std(self.beamline.distribution_fsm_h01[:,2]*1000)
        vertical_gaussian_fit = stats.norm.pdf(self.beamline.distribution_fsm_h01[:,2]*1000, y_mean, y_rms_size)

        self.ax3.hist(self.beamline.distribution_fsm_h01[:,0]*1000, bins=50, density = True, histtype = 'step')
        self.ax3.scatter(self.beamline.distribution_fsm_h01[:,0]*1000, horizontal_gaussian_fit, s=0.5)
        self.ax3.set_xlabel("X(mm)")
        self.ax3.set_xlim(-16,16)
        self.ax3.set_title("H01-FSM-001")


        x_fsm = np.arange(-16, 16, 0.125)
        simu_pdf = stats.norm.pdf(x_fsm,loc=0.01, scale=0.16)
        size_fsm_simu = len(simu_pdf)
        noise = [random.random()/4 for i in range(size_fsm_simu)]
        fsm_signal = simu_pdf + noise

        self.ax5.hist(self.beamline.distribution_fsm_h01[:,2]*1000, bins=16, density = True, histtype = 'step', label = "model")
        self.ax5.scatter(self.beamline.distribution_fsm_h01[:,2]*1000, vertical_gaussian_fit, s=0.5, label = "gaussian fit")
        self.ax5.bar(x_fsm, fsm_signal, width = 0.125, fill = False, label = "Measured with noise")
        self.ax5.set_xlabel("Y(mm)")
        self.ax5.set_xlim(-16,16)
        self.ax5.set_title("H01-FSM-001")
        self.ax5.legend()

        self.ax4.scatter(self.beamline.distribution_opt_h01[:,0]*1000, self.beamline.distribution_opt_h01[:,2]*1000, color = 'purple', s=0.5)
        self.ax4.set_xlabel("X(mm)")
        self.ax4.set_ylabel("Y(mm)")
        self.ax4.set_xlim(-16,16)
        self.ax4.set_ylim(-16,16)
        self.ax4.set_aspect('equal')
        self.ax4.grid(linestyle='--')
        self.ax4.locator_params(tight=True, nbins=10)
        self.ax4.set_title("H01-OPT-001")

# Phase space
        self.ax6.scatter(self.beamline.distribution_fsm_h01[:,0]*1000, self.beamline.distribution_fsm_h01[:,1]*1000, color = 'purple', s=0.5, label = "Hor.X/X'")
        self.ax6.scatter(self.beamline.distribution_fsm_h01[:,2]*1000, self.beamline.distribution_fsm_h01[:,3]*1000, color = 'blue', s=0.5, label = "Vert.Y/Y'")
        self.ax6.set_xlim(-16,16)
        self.ax6.set_ylim(-16,16)
        self.ax6.legend()
        self.ax6.set_title("Phase Space FSM")
        self.ax6.grid(linestyle='--')
        self.ax6.locator_params(tight=True, nbins=10)

        self.canvas = FigureCanvas(self.fi)
        self.setCentralWidget(self.canvas)
        self.canvas.draw()