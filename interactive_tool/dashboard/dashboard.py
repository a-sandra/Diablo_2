import sys
import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
#from ..build_beamline import Beamline

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, myfigure, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes_env = fig.add_subplot(111)

        #self.axes_profile = fig.add_subplot(122)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, beamline, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        widget = QtWidgets.QWidget()
        self.setCentralWidget(widget)
        grid = QtWidgets.QGridLayout()

        # Build the plot with th aperture
        beamline.aperture_beamline()

        self._sc = MplCanvas(self, width=5, height=4, dpi=100)

        #plt.show()


        #self._sc.axes_env.plot([0, beamline.maximum_s], [10,1,20,3,40])

        #grid.addWidget(QtWidgets.QLabel("param1"), 0, 0)
        #grid.addWidget(QtWidgets.QLabel("param2"), 0, 0)

        #self._sl = QtWidgets.QSlider()
        #self._sl.setMinimum(-10)
        #self._sl.setMaximum(10)
        #self._sl.setTickPosition(QtWidgets.QSlider.TicksBelow)
        #self._sl.setTickInterval(1)
        #self._sl.valueChanged.connect(self.valuechange)
        #grid.addWidget(self._sl, 0, 1)
        #grid.addWidget(self._sc, 0, 2)
        widget.setLayout(grid)
        self.show()

    def valuechange(self):
        print("value changed")
        param = self._sl.value()
        x = [i for i in range(0,100)]
        y = [i**2 + i*param - 2*param for i in x]
        self._sc.axes_env.cla()
        self._sc.axes_env.plot(x, y)
        self._sc.draw()