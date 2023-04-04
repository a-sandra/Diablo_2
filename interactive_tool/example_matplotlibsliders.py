import sys

from PyQt5.QtWidgets import QMainWindow, QApplication, QDockWidget, QWidget, QGridLayout, QSlider

from PyQt5.QtCore import Qt

import numpy as np
import matplotlib.pyplot
import matplotlib.backends.backend_qt5agg

class MainWindow(QMainWindow):

    x = np.arange(0,10, 0.1)
    cos = 0
    sin = 0

    def __init__(self):
        QMainWindow.__init__(self)

        self.figure = matplotlib.pyplot.figure()
        self.drawing = self.figure.add_subplot(111)
        self.canvas = matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg(self.figure)

        self.setCentralWidget(self.canvas)

        dock = QDockWidget ("Values")
        self.addDockWidget (Qt.RightDockWidgetArea, dock)

        sliders = QWidget ()
        sliders_grid = QGridLayout (sliders)

        def add_slider(foo, col):
           sld = QSlider(Qt.Horizontal, sliders)
           sld.setFocusPolicy(Qt.NoFocus)
           sld.valueChanged[int].connect(foo)
           sld.valueChanged.connect(self.plot)
           sliders_grid.addWidget (sld, 0, col)
        
        add_slider (foo = self.setcos, col = 0)
        add_slider (foo = self.setsin, col = 1)

        dock.setWidget (sliders)
        self.plot ()

    def setcos(self, v):
        self.cos = v / float (100)
        
    def setsin(self, v):
        self.sin = v / float (100)
        
    def plot(self):
        #self.drawing.hold(False)
        s = np.sin(self.x + self.sin)
        c = np.cos(self.x + self.cos)
        self.drawing.plot(self.x, s, 'r', self.x, c, 'r', self.x, s +c, 'b')
        self.drawing.set_ylim(-2, 2)
        self.canvas.draw()


if __name__ == "__main__":
  app = QApplication (sys.argv)
  main = MainWindow ()
  main.show ()
  sys.exit (app.exec_ ())