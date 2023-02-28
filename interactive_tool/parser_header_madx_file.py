#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:48:20 2023

@author: s.aumon
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")
#sys.path.append("/home/s.aumon/Python_Packages/Diablo_2/")
from beam import beam
import madx_utils.madx_utils as mxu
import re

madx_fname = "lattice_example\optics_71.7_MeV.twiss"
string = "@ PC               %le        0.3737500416"

match = re.search(r'^@\sPC.', string)
if match:
    print(match)
else:
    print("tagueule")


with open(madx_fname, "r") as fp:
    lines = fp.readlines()
    for line in lines:
        if re.search(r'^@\sMASS.', line):
            print(line)
            string_split = line.split()
            print(string_split[-1])

