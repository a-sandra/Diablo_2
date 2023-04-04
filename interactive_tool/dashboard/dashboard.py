import sys, copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, beamline,*args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setGeometry(100, 100, 1000, 700) 
        self.setWindowTitle('L-HEBT envelop')
        grid = QtWidgets.QGridLayout()
        self.setLayout(grid)

        self.sp = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        #self.sp.setGeometry(QtCore.QRect(190, 100, 160, 16))
        #grid.addWidget(self.sp, 5, 0)

        """
        s = beamline.dataframe_madx_sequence["S"]-beamline.dataframe_madx_sequence["L"]/2
        beamline.aperture_beamline()
        beamline.build_dataframe_element()
        self.maximum_s = max(beamline.dataframe_madx_sequence["S"])

        list_emq = [Rectangle(( beamline.df_emq["S"][i]-beamline.df_emq["L"][i]/2, -0.2/2), beamline.df_emq["L"][i], 0.2,color ='skyblue') for i in beamline.df_emq.index]
        list_stm = [Rectangle(( beamline.df_stm["S"][i]-beamline.df_stm["L"][i]/2, -0.2/2), beamline.df_stm["L"][i], 0.2,color ='limegreen') for i in beamline.df_stm.index]
        list_bpm = [Rectangle(( beamline.df_bpm["S"][i]-beamline.df_bpm["L"][i]/2, -0.2/2), beamline.df_bpm["L"][i], 0.2,color ='orange') for i in beamline.df_bpm.index]
        list_bend = [Rectangle(( beamline.df_bend["S"][i]-beamline.df_bend["L"][i]/2, -0.2/2), beamline.df_bend["L"][i], 0.2,color ='crimson') for i in beamline.df_bend.index]
        list_fsm = [Rectangle(( beamline.df_fsm["S"][i]-beamline.df_fsm["L"][i]/2, -0.2/2), beamline.df_fsm["L"][i], 0.2,color ='purple') for i in beamline.df_fsm.index]

        fi, ((self.ax1, self.ax3),(self.ax2, self.ax4))= plt.subplots(2, 2, figsize = (20,15), gridspec_kw={'height_ratios':[1,4]})

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

        #self.ax1 = fi.add_subplot(221)
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
        #self.ax2 = fi.add_subplot(223)
        self.ax2.plot(beamline.df_beam_size_along_s["S"], 3.0*beamline.df_beam_size_along_s["X"]*1000, "ob")
        self.ax2.plot(beamline.df_beam_size_along_s["S"], -3.0*beamline.df_beam_size_along_s["X"]*1000, "ob")
        self.ax2.step(s, beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax2.step(s, -beamline.aperture_x, where = 'pre', color = 'k', label = "_nolabel_")

        for i in list_emq:
            nq=copy.copy(i)
            self.ax3.add_patch(nq)
        for i in list_stm:
            ns=copy.copy(i)
            self.ax3.add_patch(ns)
        for i in list_bpm:
            nb=copy.copy(i)
            self.ax3.add_patch(nb)
        for i in list_bend:
            nbend=copy.copy(i)
            self.ax3.add_patch(nbend)
        for i in list_fsm:
            nfsm=copy.copy(i)
            self.ax3.add_patch(nfsm)

        self.ax3.plot([0, self.maximum_s],[0,0], "k", linewidth=0.5)
        self.ax3.set_ylim(-0.2,0.6)
        self.ax3.grid()
        self.ax3.axes.yaxis.set_visible(False)
        self.ax3.axes.xaxis.set_visible(False)
        self.ax3.spines['top'].set_visible(False)
        self.ax3.spines['bottom'].set_visible(False)
        self.ax3.spines['left'].set_visible(False)
        self.ax3.spines['right'].set_visible(False)
        self.ax3.set_xlim(0,6)
        self.ax4.step(s, beamline.aperture_y, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax4.step(s, -beamline.aperture_y, where = 'pre', color = 'k', label = "_nolabel_")
        self.ax4.set_xlabel("s(m)")
        self.ax4.set_xlim(0,6)
        self.ax4.set_ylim(-20,20)
        self.ax4.set_ylabel("Vertical (mm)")

        central_widget = FigureCanvas(fi)
        self.setCentralWidget(central_widget)
        """