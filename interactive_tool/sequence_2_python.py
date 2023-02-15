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
import optics_utils.transfer_matrices as opu

class Beamline(object):
    def __init__(self, sequence_file={}, kin_energy = 70.0):
        self.kinetic_energy = kin_energy
        self.sequence_file = sequence_file
        self.full_transfert_matrice = []
        print(bool(self.sequence_file))
        if bool(self.sequence_file) is True:
            print("beamline created")
            self.full_transfert_matrice = self.build_beamline_from_optics_file(self.sequence_file)
        else:
            print("No Twiss file has been provided, the beam line is empty")
        print(self.full_transfert_matrice)
        
    
    def build_beamline_from_optics_file(self, sequence_file = ""):
        self.sequence_file = sequence_file
        if not self.sequence_file:
            print("No Twiss file")
        else:
            print("sequence file given") # method to construct the beam line
            self.madx_sequence = mxu.get_twiss(self.sequence_file)
            #print(self.madx_sequence[:10])
            new_df = self.madx_sequence[:3]
            #print(new_df)
            self.full_sequence = []
            for index in np.arange(len(new_df)):
                #print(index)
                #print(new_df.iloc[index]["KEYWORD"])
                self.check_element_type(new_df.iloc[index])
        return self.full_sequence # build a dataframe Name, Keyword, l, s, k1l, angle
    
    def build_beamline_transfer_matrix(self):
        self.full_transfert_matrice = []
        for index in range(1, len(self.full_transfert_matrice)):
            ftm = self.full_transfert_matrice[index].dot(self.full_transfert_matrice[index-1])
        return ftm

    def check_element_type(self, element):
        if element["KEYWORD"]=="QUADRUPOLE":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a quad") 
            quad = Quadrupole(self, element["K1L"],  element["NAME"], element["L"])
            print(element["NAME"])
            self.full_transfert_matrice.append(quad.quadrupole_transfert_matrix)

        elif element["KEYWORD"]=="SBEND":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is sbend")
            sbend = Dipole(self, element["L"], 0.785398, element["NAME"])
            print(element["NAME"])
            self.full_transfert_matrice.append(sbend.dipole_transfert_matrix)

        elif element["KEYWORD"]=="DRIFT":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a drift")
            drift = Drift(self,  element["L"], element["NAME"])
            print(element["NAME"])
            self.full_transfert_matrice.append(drift.drift_transfert_matrix)

        elif element["KEYWORD"]=="INSTRUMENT":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is an instrument, a drift")
            drift =Drift(self,  element["L"], element["NAME"])
            print(element["NAME"])
            self.full_transfert_matrice.append(drift.drift_transfert_matrix)

        elif element["KEYWORD"]=="MARKER":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a MARKER, a drift")
            pass

        elif element["KEYWORD"]=="MONITOR":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a monitor, a drift")
            drift = Drift(self,  element["L"], element["NAME"])
            print(element["NAME"])
            self.full_transfert_matrice.append(drift.drift_transfert_matrix)

        elif element["KEYWORD"]=="KICKER": # let's consider the kicker as a dritt at the moment
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a steerer")
            drift = Drift(self,  element["L"], element["NAME"])
            print(element["NAME"])
            self.full_transfert_matrice.append(drift.drift_transfert_matrix)     


class Quadrupole(object):
    def __init__(self, beamline, k1l, name = "noname_quadrupole", magnetic_length = 1):
        self.beamline = beamline
        self.kinetic_energy = beamline.kinetic_energy
        self.magnetic_quadrupole_length = magnetic_length
        self.quadrupole_name = name
        self.quadrupole_strength = k1l/magnetic_length
        if self.quadrupole_strength >= 0 :
            self.quadrupole_transfert_matrix = opu.mqf(self.magnetic_quadrupole_length, self.magnetic_quadrupole_length)
        else:
            self.quadrupole_transfert_matrix = opu.mqd(self.magnetic_quadrupole_length, self.magnetic_quadrupole_length)

class Drift(object):
    def __init__(self, beamline, dr_length = 1, name = "noname_drift"):
        self.drift_length = dr_length
        self.drift_name = name
        self.drift_transfert_matrix = np.array([[1.0, 0.0], [0, self.drift_length]])
        
class Dipole(object):
    def __init__(self, beamline, dip_length = 1, angle = 0.785398, name = "noname_dipole"):
        self.beamline = beamline
        self.kinetic_energy = beamline.kinetic_energy
        self.bending_angle = angle
        self.dipole_magnetic_length = dip_length
        self.rho = self.dipole_magnetic_length/self.bending_angle
        self.dipole_name = name
        self.dipole_transfert_matrix = np.array([[np.cos(self.bending_angle), self.rho*np.sin(self.bending_angle)], [-1/self.bending_angle*np.sin(self.bending_angle), np.cos(self.bending_angle)]])

class Steerer(object):
    def __init__(self, beamline, name = "noname_steerer", st_length = 1):
        self.steerer_length = st_length
        self.steerer_name = name
        self.drift_transfert_matrix = np.array([[1.0, 0.0], [0, dr_length]])




    
    
    
    
    


        


        
        
        
        