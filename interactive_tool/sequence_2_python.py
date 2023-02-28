#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:48:20 2023

@author: s.aumon
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys, re

sys.path.append(r"C:\Users\s.aumon\OneDrive - Advanced Oncotherapy plc\Python_Packages\Diablo_2")
#sys.path.append("/home/s.aumon/Python_Packages/Diablo_2/")
from beam import beam
import madx_utils.madx_utils as mxu
import optics_utils.transfer_matrices as opu

class Beamline(object):
    def __init__(self, sequence_file={}, kin_energy = 70.0):
        self.kinetic_energy = kin_energy
        self.brho = 0
        self.mass = 0
        self.sequence_file = sequence_file
        self.full_transfert_matrice = []
        self.quadrupole_list = []
        self.sbend_list = []
        self.steerer_list = []
        self.non_magnetic_element = []

        # check if an input file is given.
        # this part should just build the beam line.
        if bool(self.sequence_file) is True:
            print("beamline created")
            self.build_beamline_from_sequence_file(self.sequence_file)
        else:
            print("No Twiss file has been provided, the beam line is empty")

    # get the kinetic energy of the beam
    def get_kinetic_energy(self):
        return self.kinetic_energy
    
    # this method is meant to build the beam line from the sequence file
    def build_beamline_from_sequence_file(self, sequence_file = ""):
        self.sequence_file = sequence_file
        if not self.sequence_file:
            print("No Twiss file")
        else:
            print("sequence file given")
            self.dataframe_madx_sequence = mxu.get_twiss(self.sequence_file) # Here, we supposed that the sequence file is a MADX sequence file
            # this could be done as an external function.
            with open(self.sequence_file, "r") as fp:
                lines = fp.readlines()
                for line in lines:
                    if re.search(r'^@\sPC.', line):
                        string_split = line.split()
                        self.pc = float(string_split[-1])
                        self.brho = 3.33*self.pc # only valid for proton
                    elif re.search(r'^@\sMASS.', line):
                        string_split0 = line.split()
                        self.mass = float(string_split0[-1])
            self.kinetic_energy = 1000*(np.sqrt(self.mass**2 + self.pc**2)-self.mass)

            #check element by element, in order to build a dataframe with NAME, KEYWORD, etc.. And Transfer matrix.
            self.main_dataframe_sequence = pd.DataFrame({"NAME": pd.Series(dtype="str"), "KEYWORD": pd.Series(dtype="str"), "S": pd.Series(dtype="float"),
                                                        "L": pd.Series(dtype="float"), "K1": pd.Series(dtype="float"), "ANGLE": pd.Series(dtype="float"),
                                                        "TRANSFER_MATRIX": pd.Series(dtype="float")})
            self.data = []
            #for index in np.arange(len(self.dataframe_madx_sequence)):
            #    self.check_element_type(self.dataframe_madx_sequence.iloc[index])
    
    def build_beamline_transfer_matrix(self):
        self.full_transfert_matrice = []
        for index in range(1, len(self.full_transfert_matrice)):
            ftm = self.full_transfert_matrice[index].dot(self.full_transfert_matrice[index-1])
        return ftm

    def check_element_type(self, element):
        if element["KEYWORD"]=="QUADRUPOLE":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a quad") 
            quad = Quadrupole(self, element["K1L"],  element["NAME"], element["L"])
            self.quadrupole_list.append(quad)
            if (element["K1L"]>=0):
                k1 = element["K1L"]/element["L"]
                self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"]/element["L"], 
                                element["ANGLE"], opu.mqf(k1, element["L"])])
            else:
                k1 = element["K1L"]/element["L"]
                self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"]/element["L"], 
                                element["ANGLE"], opu.mqd(k1, element["L"])])


        elif element["KEYWORD"]=="SBEND":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is sbend")
            sbend = Dipole(self, element["L"], 0.785398, element["NAME"])
            self.sbend_list.append(sbend)
            self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"], element["ANGLE"], opu.md(element["L"])])

        elif element["KEYWORD"]=="DRIFT":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a drift")
            drift = Drift(self,  element["L"], element["NAME"])
            self.non_magnetic_element.append(drift)
            self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"], element["ANGLE"], opu.md(element["L"])])

        elif element["KEYWORD"]=="INSTRUMENT":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is an instrument, a drift")
            drift =Drift(self,  element["L"], element["NAME"])
            self.non_magnetic_element.append(drift)
            self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"], element["ANGLE"], opu.md(element["L"])])

        elif element["KEYWORD"]=="MARKER":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a MARKER, a drift")
            pass

        elif element["KEYWORD"]=="MONITOR":
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a monitor, a drift")
            drift = Drift(self,  element["L"], element["NAME"])
            self.non_magnetic_element.append(drift)
            self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"], element["ANGLE"], opu.md(element["L"])])

        elif element["KEYWORD"]=="KICKER": # let's consider the kicker as a drift at the moment
            #print(element["NAME"] +" "+ element["KEYWORD"]+", this element is a steerer")
            steerer = Drift(self,  element["L"], element["NAME"])
            self.steerer_list.append(steerer)
            self.data.append([element["NAME"], element["KEYWORD"], element["S"], element["L"], element["K1L"], element["ANGLE"], opu.md(element["L"])])

    def update_magnet_parameter(self, keyword, list_of_magnet_to_update, keyword_parameter_to_update, parameter_to_update):
        dataframe_as_keyword = self.dataframe_madx_sequence[self.dataframe_madx_sequence["KEYWORD"]==keyword]
        #print(dataframe_as_keyword)
        for magnet, field  in zip(list_of_magnet_to_update, parameter_to_update):
            index_magnet = self.dataframe_madx_sequence[self.dataframe_madx_sequence['NAME']==magnet].index.values.astype(int)[0]
            self.dataframe_madx_sequence.at[index_magnet, keyword_parameter_to_update] = field
        #print(self.dataframe_madx_sequence[self.dataframe_madx_sequence["KEYWORD"]==keyword])

"""
print(df_optics[df_optics["KEYWORD"]=="QUADRUPOLE"])

for quadrupole_name, strength in zip(updated_HEBT_quadrupole, updated_strength):
    index_quadrupole = df_optics[df_optics['NAME']==quadrupole_name].index.values.astype(int)[0]
    df_optics.at[index_quadrupole, "K1L"] = strength

print(df_optics[df_optics["KEYWORD"]=="QUADRUPOLE"])

"""

class Quadrupole(object):
    def __init__(self, beamline, k1l, name = "noname_quadrupole", magnetic_length = 1): # MADX gives K1L, strength*length, thus it has to be converted into strength, gradient
        self.beamline = beamline
        self.kinetic_energy = beamline.kinetic_energy
        self.magnetic_quadrupole_length = magnetic_length
        self.quadrupole_name = name
        self.integrated_strength = k1l
        self.quadrupole_strength = k1l/magnetic_length
        self.gradient = self.quadrupole_strength*beamline.brho

# def mq(k1x, length), strength, length

    def quadrupole_transfer_matrix(self):
        if self.quadrupole_strength >= 0 :
            self.quadrupole_transfert_matrix = opu.mqf(self.quadrupole_strength, self.magnetic_quadrupole_length)
        else:
            self.quadrupole_transfert_matrix = opu.mqd(self.quadrupole_strength, self.magnetic_quadrupole_length)

## setter
    def set_gradient(self, gradient):
        self.gradient = gradient
        return self.gradient

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




    
    
    
    
    


        


        
        
        
        