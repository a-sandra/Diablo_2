from _dependencies import pandas as _pd
from _dependencies import numpy as _np

def mqd(k1x, length):
    return _np.array([[ _np.cosh(_np.sqrt(_np.abs(k1x))*length), 1/_np.sqrt(_np.abs(k1x))*_np.sinh(_np.sqrt(_np.abs(k1x))*length)], [_np.sqrt(_np.abs(k1x))*_np.sinh(_np.sqrt(_np.abs(k1x))*length), _np.cosh(_np.sqrt(_np.abs(k1x))*length)]])
    #return _np.array([[ 1, 0], [k1x*length, 1]])

def md(length):
    return _np.array([[1, length], [0, 1]])

def mqf(k1x, length):
    return _np.array([[ _np.cos(_np.sqrt(k1x)*length), 1/_np.sqrt(k1x)*_np.sin(_np.sqrt(k1x)*length)], [-_np.sqrt(k1x)*_np.sin(_np.sqrt(k1x)*length), _np.cos(_np.sqrt(k1x)*length)]])
    #return _np.array([[ 1, 0], [k1x*length, 1]])

