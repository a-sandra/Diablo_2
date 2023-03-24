import sys, copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, beamline,*args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setGeometry(100, 100, 680, 500)

        s = beamline.dataframe_madx_sequence["S"]-beamline.dataframe_madx_sequence["L"]/2
        beamline.aperture_beamline()
        beamline.build_dataframe_element()

        list_emq = [Rectangle(( self.df_emq["S"][i]-self.df_emq["L"][i]/2, -0.2/2), self.df_emq["L"][i], 0.2,color ='skyblue') for i in self.df_emq.index]
        list_stm = [Rectangle(( self.df_stm["S"][i]-self.df_stm["L"][i]/2, -0.2/2), self.df_stm["L"][i], 0.2,color ='limegreen') for i in self.df_stm.index]
        list_bpm = [Rectangle(( self.df_bpm["S"][i]-self.df_bpm["L"][i]/2, -0.2/2), self.df_bpm["L"][i], 0.2,color ='orange') for i in self.df_bpm.index]
        list_bend = [Rectangle(( self.df_bend["S"][i]-self.df_bend["L"][i]/2, -0.2/2), self.df_bend["L"][i], 0.2,color ='crimson') for i in self.df_bend.index]
        list_fsm = [Rectangle(( self.df_fsm["S"][i]-self.df_fsm["L"][i]/2, -0.2/2), self.df_fsm["L"][i], 0.2,color ='purple') for i in self.df_fsm.index]

        fi, ((self.ax1, self.ax3),(self.ax2, self.ax4))= plt.subplots(2, 2, figsize = (20,15), gridspec_kw={'height_ratios':[1,4]})

        for i in self.list_emq:
            nq=copy.copy(i)
            self.ax1.add_patch(nq)
        for i in self.list_stm:
            ns=copy.copy(i)
            self.ax1.add_patch(ns)
        for i in self.list_bpm:
            nb=copy.copy(i)
            self.ax1.add_patch(nb)
        for i in self.list_bend:
            nbend=copy.copy(i)
            self.ax1.add_patch(nbend)
        for i in self.list_fsm:
            nfsm=copy.copy(i)
            self.ax1.add_patch(nfsm)

        self.ax1.plot([0, beamline.maximum_s],[0,0], "k", linewidth=0.5)
        self.ax1.set_ylim(-0.2,0.6)
        self.ax1.set_xlim(0,6)
        self.ax1.grid()
        self.ax1.axes.yaxis.set_visible(False)
        self.ax1.axes.xaxis.set_visible(False)
        self.ax1.spines['top'].set_visible(False)
        self.ax1.spines['bottom'].set_visible(False)
        self.ax1.spines['left'].set_visible(False)
        self.ax1.spines['right'].set_visible(False)
        self.ax2 = fi.add_subplot(111)
        self.ax2.plot(beamline.df_beam_size_along_s["S"], 3.0*beamline.df_beam_size_along_s["X"]*1000, "ob")
        self.ax2.plot(beamline.df_beam_size_along_s["S"], -3.0*beamline.df_beam_size_along_s["X"]*1000, "ob")


        central_widget = FigureCanvas(fi)
        self.setCentralWidget(central_widget)

        self.show()

