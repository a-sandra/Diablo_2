import sys
import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        widget = QtWidgets.QWidget()
        self.setCentralWidget(widget)
        grid = QtWidgets.QGridLayout()
        self._sc = MplCanvas(self, width=5, height=4, dpi=100)
        self._sc.axes.plot([0,1,2,3,4], [10,1,20,3,40])
        grid.addWidget(QtWidgets.QLabel("param1"), 0, 0)
        self._sl = QtWidgets.QSlider()
        self._sl.setMinimum(-10)
        self._sl.setMaximum(10)
        self._sl.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self._sl.setTickInterval(1)
        self._sl.valueChanged.connect(self.valuechange)
        grid.addWidget(self._sl, 0, 1)
        grid.addWidget(self._sc, 0, 2)
        widget.setLayout(grid)
        self.show()

    def valuechange(self):
        print("value changed")
        param = self._sl.value()
        x = [i for i in range(0,100)]
        y = [i**2 + i*param - 2*param for i in x]
        self._sc.axes.cla()
        self._sc.axes.plot(x, y)
        self._sc.draw()