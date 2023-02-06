#from _dependencies import _plt
import matplotlib.pyplot as _plt
font = {'size'   : 14}
_plt.matplotlib.rc('font', **font)

from _dependencies import beam as bm

def hist(column, xlab, ylab, legend, place="upper left"):
    n_occ, x_bin, list_obj = _plt.hist(column, bins='auto', label=legend)
    _plt.grid(linestyle="--")
    _plt.xlabel(xlab)
    _plt.ylabel(ylab)
    _plt.legend(loc = place)
