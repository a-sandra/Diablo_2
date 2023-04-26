from _dependencies import pandas as _pd
from _dependencies import numpy as _np

#def mqd(k1x, length):
#    return _np.array([[ _np.cosh(_np.sqrt(_np.abs(k1x))*length), 1/_np.sqrt(_np.abs(k1x))*_np.sinh(_np.sqrt(_np.abs(k1x))*length)], [_np.sqrt(_np.abs(k1x))*_np.sinh(_np.sqrt(_np.abs(k1x))*length), _np.cosh(_np.sqrt(_np.abs(k1x))*length)]])
#    #return _np.array([[ 1, 0], [k1x*length, 1]])

def mqd(k1x, length):
    k1l = _np.sqrt(_np.abs(k1x))*length
    sqrt_k1 = _np.sqrt(_np.abs(k1x))
    #return _np.array([[ _np.cosh(k1l), 1/sqrt_k1*_np.sinh(k1l)], [sqrt_k1*_np.sinh(k1l), _np.cosh(k1l)]])

    return _np.array([[_np.cosh(k1l), 1/sqrt_k1*_np.sinh(k1l), 0, 0],
                  [sqrt_k1*_np.sinh(k1l), _np.cosh(k1l), 0, 0],
                  [0, 0, _np.cos(k1l), 1/sqrt_k1*_np.sin(k1l)],
                  [0, 0, -sqrt_k1*_np.sin(k1l), _np.cos(k1l)]])

def md(length):
    #return _np.array([[1, length], [0, 1]])
    return _np.array([[1, length, 0, 0], [0, 1, 0, 0], [0, 0, 1, length], [0, 0, 0, 1]])

def mqf(k1x, length):
    k1l = _np.sqrt(_np.abs(k1x))*length
    sqrt_k1 = _np.sqrt(_np.abs(k1x))
    #return _np.array([[ _np.cos(k1l), 1/sqrt_k1*_np.sin(k1l)], [ -sqrt_k1*_np.sin(k1l), _np.cos(k1l)]])

    return _np.array([[ _np.cos(k1l), 1/sqrt_k1*_np.sin(k1l), 0, 0],
                      [ -sqrt_k1*_np.sin(k1l), _np.cos(k1l), 0, 0],
                      [0, 0, _np.cosh(k1l), 1/sqrt_k1*_np.sinh(k1l)],
                      [0, 0, sqrt_k1*_np.sinh(k1l), _np.cosh(k1l)]])

def sector_dipole(length, rho):
    theta = length/rho
    return _np.array([[ _np.cos(theta), _np.sin(theta)*rho, 0, 0],
                     [ -1/rho*_np.sin(theta), _np.cos(theta), 0, 0],
                     [0, 0, 1, length],
                     [0, 0, 0, 1]])

def unity():
    #return _np.array([[1, 0], [0, 1]])
    return _np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0,], [0, 0, 0, 1]])


