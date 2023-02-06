

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append("/Users/pickle/Documents/Diablo_2/")
from beam import beam
import madx_utils.madx_utils as mxu
import optics_utils.transfer_matrices as opu
import pickle as pkl
import copy
from matplotlib.patches import Rectangle

#from plot_aperture import fi


from ipywidgets import interact
import ipywidgets as widgets
from matplotlib.widgets import Slider, Button

class Interactive(object):
    def __init__(self,  file_input_distribution, input_lattice, lattice = "H01"):
        
        self.bm = 0
        self._k1q11 = 0.0
        self._k1q12 = 5.0
        self._k1q13 = -5.0
        self._k1q14 = 5.0
        self._k1q15 = -5.0
        
        self._k1q21 = 0.0
        
        self._k1q31 = 0.0
        self._k1q32 = 0.0
        self._k1q33 = 0.0
        self._k1q34 = 0.0
        self._k1q35 = 0.0

        self.input_beamfile = file_input_distribution
        self.bm = beam.Beam()
        self.bm.read_distribution(self.input_beamfile)
        self._sigma = 0.0
        
        if lattice == "H01":
            df_optics = mxu.get_twiss(input_lattice);
            df_emq = df_optics[df_optics["KEYWORD"]=="QUADRUPOLE"]
            df_fsm = df_optics[df_optics["NAME"]=="H01_FSM_01"]
            self._length_quad = 0.14
            
            #Defining the distance between the quadrupoles
            length_drift_q1 = df_optics[df_optics["NAME"]=="H01_MECHAL_011"]["S"].values[0]
            length_q1_q2 = df_optics[df_optics["NAME"]=="H01_MECHAL_021"]["S"].values[0] - df_optics[df_optics["NAME"]=="H01_MECHAL_012"]["S"].values[0]
            length_q2_q3 = df_optics[df_optics["NAME"]=="H01_MECHAL_031"]["S"].values[0] - df_optics[df_optics["NAME"]=="H01_MECHAL_022"]["S"].values[0]
            length_q3_q4 = df_optics[df_optics["NAME"]=="H01_MECHAL_041"]["S"].values[0] - df_optics[df_optics["NAME"]=="H01_MECHAL_032"]["S"].values[0]
            length_q4_q5 = df_optics[df_optics["NAME"]=="H01_MECHAL_051"]["S"].values[0] - df_optics[df_optics["NAME"]=="H01_MECHAL_042"]["S"].values[0]
            length_q4_fsm = df_optics[df_optics["NAME"]=="H01_FSM_01"]["S"].values[0] -  df_optics[df_optics["NAME"]=="H01_MECHAL_042"]["S"].values[0]
            
            #Defining the drift matrices between the quadrupoles
            self.drift_q1 = opu.md(length_drift_q1)
            self.drift_q1_q2 = opu.md(length_q1_q2)
            self.drift_q2_q3 = opu.md(length_q2_q3)
            self.drift_q3_q4 = opu.md(length_q3_q4)
            self.drift_q4_q5 = opu.md(length_q4_q5)
            self.drift_q4_fsm = opu.md(length_q4_fsm)
            
            self.s_q1 = df_emq[df_emq["NAME"]=="H01_EMQ_01"]["S"].values[0]
            self.s_q2 = df_emq[df_emq["NAME"]=="H01_EMQ_02"]["S"].values[0]
            self.s_q3 = df_emq[df_emq["NAME"]=="H01_EMQ_03"]["S"].values[0]
            self.s_q4 = df_emq[df_emq["NAME"]=="H01_EMQ_04"]["S"].values[0]
            self.s_q5 = df_emq[df_emq["NAME"]=="H01_EMQ_05"]["S"].values[0]
            self.s_fms = df_fsm[df_fsm["NAME"]=="H01_FSM_01"]["S"].values[0]
            #print(self.s_fms)
            
            hor_distribution = self.bm.distribution[["X(mm)", "XP(mrad)"]]*0.001
            vert_distribution = self.bm.distribution[["Y(mm)", "YP(mrad)"]]*0.001
            
            self.transposed_hor = hor_distribution.T
            self.transposed_vert = vert_distribution.T

            self.start_std_x = np.std(hor_distribution["X(mm)"])
            self.start_std_y = np.std(vert_distribution["Y(mm)"])
            
        else:
            pass

#%%    
    @property
    def k1q11(self):
        return self._k1q11
    @k1q11.setter
    def k1q11(self, value):
        self._k1q11 = value
        
    @property
    def k1q12(self):
        return self._k1q12
    @k1q12.setter
    def k1q12(self, value):
        self._k1q12 = value
    
    @property
    def k1q13(self):
        return self._k1q13
    # a setter function
    @k1q13.setter
    def k1q13(self, a):
        self._k1q13 = a
         
    @property
    def k1q14(self):
        return self._k1q14
    # a setter function
    @k1q14.setter
    def k1q14(self, a):
        self._k1q14 = a
         
    @property
    def k1q15(self):
        return self._k1q15
    # a setter function
    @k1q15.setter
    def k1q15(self, a):
        self._k1q15 = a
      
    @property
    def k1q21(self):
        return self._k1q21
    # a setter function
    @k1q21.setter
    def k1q21(self, a):
        self._k1q21 = a

    @property
    def k1q31(self):
        return self._k1q31
    # a setter function
    @k1q31.setter
    def k1q31(self, a):
        self._k1q31 = a
    
    @property
    def k1q32(self):
        return self._k1q32
    # a setter function
    @k1q32.setter
    def k1q32(self, a):
        self._k1q32 = a
         
    @property
    def k1q33(self):
        return self._k1q33
    # a setter function
    @k1q33.setter
    def k1q33(self, a):
        self._k1q33 = a
         
    @property
    def k1q34(self):
        return self._k1q34
    # a setter function
    @k1q34.setter
    def k1q34(self, a):
        self._k1q34 = a
         
    @property
    def k1q35(self):
        return self._k1q35
    # a setter function
    @k1q35.setter
    def k1q35(self, a):
        self._k1q35 = a
#%%
        
   # def computeHorizontalOptics(self):
   #     #Horizontal
   #     mid_hor_q1 = opu.mqd(self._k1q11, self._length_quad).dot(self.drift_q1).dot(self.transposed_hor)
   #     mid_hor_q2 = opu.mqf(self._k1q12, self._length_quad).dot(self.drift_q1_q2).dot(mid_hor_q1)
   #     mid_hor_q3 = opu.mqd(self._k1q13, self._length_quad).dot(self.drift_q2_q3).dot(mid_hor_q2)
   #     mid_hor_q4 = opu.mqf(self._k1q14, self._length_quad).dot(self.drift_q3_q4).dot(mid_hor_q3)
   #     mid_hor_q5 = opu.mqd(self._k1q15, self._length_quad).dot(self.drift_q4_q5).dot(mid_hor_q4)
   #     self.fsm_hor =self.drift_q4_fsm.dot(mid_hor_q4)

    #    self.sx = np.array([[0, self.start_std_x],[self.s_q1, np.std(mid_hor_q1[0,:])],
    #                        [self.s_q2, np.std(mid_hor_q2[0,:])],[self.s_q3, np.std(mid_hor_q3[0,:])], 
    #                        [self.s_q4, np.std(mid_hor_q4[0,:])],[self.s_q5, np.std(mid_hor_q5[0,:])]])
        
    def computeHorizontalOptics(self):
        #Horizontal
        mid_hor_q1 = self.quadrupoleMatrix(self._k1q11, self._length_quad).dot(self.drift_q1).dot(self.transposed_hor)
        mid_hor_q2 = self.quadrupoleMatrix(self._k1q12, self._length_quad).dot(self.drift_q1_q2).dot(mid_hor_q1)
        mid_hor_q3 = self.quadrupoleMatrix(self._k1q13, self._length_quad).dot(self.drift_q2_q3).dot(mid_hor_q2)
        mid_hor_q4 = self.quadrupoleMatrix(self._k1q14, self._length_quad).dot(self.drift_q3_q4).dot(mid_hor_q3)
        mid_hor_q5 = self.quadrupoleMatrix(self._k1q15, self._length_quad).dot(self.drift_q4_q5).dot(mid_hor_q4)
        self.fsm_hor =self.drift_q4_fsm.dot(mid_hor_q4)

        self.sx = np.array([[0, self.start_std_x],[self.s_q1, np.std(mid_hor_q1[0,:])],
                            [self.s_q2, np.std(mid_hor_q2[0,:])],[self.s_q3, np.std(mid_hor_q3[0,:])], 
                            [self.s_q4, np.std(mid_hor_q4[0,:])],[self.s_q5, np.std(mid_hor_q5[0,:])]])
        
    def quadrupoleMatrix(self, k1, length):
        if k1>= 0 :
            return opu.mqf(k1, length)
        else:
            return opu.mqd(k1, length)
        
        
    def computeVerticalOptics(self):
        #Vertical 
        mid_vert_q1 = self.quadrupoleMatrix(-self._k1q11, self._length_quad).dot(self.drift_q1).dot(self.transposed_vert)
        mid_vert_q2 = self.quadrupoleMatrix(-self._k1q12, self._length_quad).dot(self.drift_q1_q2).dot(mid_vert_q1)
        mid_vert_q3 = self.quadrupoleMatrix(-self._k1q13, self._length_quad).dot(self.drift_q2_q3).dot(mid_vert_q2)
        mid_vert_q4 = self.quadrupoleMatrix(-self._k1q14, self._length_quad).dot(self.drift_q3_q4).dot(mid_vert_q3)
        mid_vert_q5 = self.quadrupoleMatrix(-self._k1q15, self._length_quad).dot(self.drift_q4_q5).dot(mid_vert_q4)
        self.fsm_ver = self.drift_q4_fsm.dot(mid_vert_q4)
            
                
        self.sy = np.array([[0, self.start_std_y],[self.s_q1, np.std(mid_vert_q1[0,:])],
                            [self.s_q2, np.std(mid_vert_q2[0,:])],[self.s_q3,np.std(mid_vert_q3[0,:])], 
                            [self.s_q4, np.std(mid_vert_q4[0,:])],[self.s_q5, np.std(mid_vert_q5[0,:])]])
    
    def sigmaIntoDataframe(self):
        self.sigma = pd.DataFrame(data = {"S": self.sx[:,0], "SIGMA_X": self.sx[:,1], "SIGMA_Y": self.sy[:,1]})


    def plotSimple(self):
        optics_file = "aperture_H01.twiss"
        aperture_file = mxu.get_twiss(optics_file)
        s = aperture_file["S"]-aperture_file["L"]/2
        aperture_x = 1000*aperture_file["APER_1"].replace(0, 0.016)
        aperture_y = 1000*aperture_file["APER_2"].replace(0, 0.016)

        tw = mxu.get_twiss(optics_file)
        maximum_s = max(tw["S"])
        df_emq = tw[tw["KEYWORD"]=="QUADRUPOLE"]
        df_stm = tw[tw["KEYWORD"]=="KICKER"]
        df_bpm = tw[tw["NAME"].str.contains(r"BPM")]
        df_bend = tw[tw["NAME"].str.contains(r"EMD")]
        df_fsm = tw[tw["NAME"].str.contains(r"FSM")]
        
        list_emq = [Rectangle(( df_emq["S"][i]-df_emq["L"][i]/2, -0.2/2),df_emq["L"][i], 0.2,color ='skyblue') for i in df_emq.index]
        list_stm = [Rectangle(( df_stm["S"][i]-df_stm["L"][i]/2, -0.2/2),df_stm["L"][i], 0.2,color ='limegreen') for i in df_stm.index]
        list_bpm = [Rectangle(( df_bpm["S"][i]-df_bpm["L"][i]/2, -0.2/2),df_bpm["L"][i], 0.2,color ='orange') for i in df_bpm.index]
        list_bend = [Rectangle(( df_bend["S"][i]-df_bend["L"][i]/2, -0.2/2),df_bend["L"][i], 0.2,color ='crimson') for i in df_bend.index]
        list_fsm = [Rectangle(( df_fsm["S"][i]-df_fsm["L"][i]/2, -0.2/2),df_fsm["L"][i], 0.2,color ='purple') for i in df_fsm.index]

        fi, ((ax1, ax3),(ax2, ax4),(ax5, ax6))= plt.subplots(3, 2, figsize = (20,15), gridspec_kw={'height_ratios':[1,4,5]})
        for i in list_emq:
            nq=copy.copy(i)
            ax1.add_patch(nq) 
        for i in list_stm:
            ns=copy.copy(i)
            ax1.add_patch(ns)    
        for i in list_bpm:
            nb=copy.copy(i)
            ax1.add_patch(nb)
        for i in list_bend:
            nbend=copy.copy(i)
            ax1.add_patch(nbend)
        for i in list_fsm:
            nfsm=copy.copy(i)
            ax1.add_patch(nfsm)
            
        ax1.plot([0, maximum_s],[0,0], "k", linewidth=0.5)
        ax1.set_ylim(-0.2,0.6)
        ax1.set_xlim(0,6)
        ax1.grid()
        ax1.axes.yaxis.set_visible(False)
        ax1.axes.xaxis.set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        
        ax2.step(s, aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        ax2.step(s, -aperture_x, where = 'pre', color = 'k', label = "_nolabel_")
        ax2.set_xlabel("s(m)")
        ax2.set_xlim(0,6)
        ax2.set_ylim(-30,30)
        ax2.set_ylabel("Horizontal (mm)")

        for i in list_emq:
            nq=copy.copy(i)
            ax3.add_patch(nq) 
        for i in list_stm:
            ns=copy.copy(i)
            ax3.add_patch(ns)    
        for i in list_bpm:
            nb=copy.copy(i)
            ax3.add_patch(nb)
        for i in list_bend:
            nbend=copy.copy(i)
            ax3.add_patch(nbend)
        for i in list_fsm:
            nfsm=copy.copy(i)
            ax3.add_patch(nfsm)
                
        ax3.plot([0, maximum_s],[0,0], "k", linewidth=0.5)
        ax3.set_ylim(-0.2,0.6)
        ax3.grid()
        ax3.axes.yaxis.set_visible(False)
        ax3.axes.xaxis.set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.spines['bottom'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.set_xlim(0,6)
        
        ax4.step(s, aperture_y, where = 'pre', color = 'k', label = "_nolabel_")
        ax4.step(s, -aperture_y, where = 'pre', color = 'k', label = "_nolabel_")
        ax4.set_xlabel("s(m)")
        ax4.set_xlim(0,6)
        ax4.set_ylim(-20,20)
        ax4.set_ylabel("Vertical (mm)")
        
        ax2.plot(self.sigma["S"], 3.0*self.sigma["SIGMA_X"]*1000, "-ro", label = "X")
        ax2.plot(self.sigma["S"], -3.0*self.sigma["SIGMA_X"]*1000, "-ro")
        
        ax4.plot(self.sigma["S"], 3.0*self.sigma["SIGMA_Y"]*1000, "-bo", label = "Y")
        ax4.plot(self.sigma["S"], -3.0*self.sigma["SIGMA_Y"]*1000, "-bo")
        
        ax5.scatter(self.fsm_hor[0,:]*1000, self.fsm_hor[1,:]*1000,label = "Hor. H01-FSM-001")
        ax5.scatter(self.fsm_ver[0,:]*1000, self.fsm_ver[1,:]*1000, label = "Vert. H01-FSM-001")
        ax5.set_xlabel("X,Y(mm)")
        ax5.set_ylabel("XP,YP (mrad)")
        ax5.grid(linestyle = "--")
        ax5.legend()
        
        ax6.hist(self.fsm_hor[0,:]*1000, bins=50, histtype = "step", label = "X at FMS")
        ax6.hist(self.fsm_ver[0,:]*1000, bins=50, histtype = "step", label = "Y at FSM")
        ax6.set_xlabel("X, Y(mm)")
        ax6.set_ylabel("#")
        ax6.grid(linestyle = "--")
        ax6.legend()
        
        fi.subplots_adjust(wspace=0.2, hspace=0.15)

    def combineComputeOpticsplotSimple(self, k1, k2, k3, k4, k5):
        self._k1q11 = k1
        self._k1q12 = k2
        self._k1q13 = k3
        self._k1q14 = k4
        self._k1q15 = k5
        self.computeHorizontalOptics()
        self.computeVerticalOptics()
        self.sigmaIntoDataframe()
        self.plotSimple()

        
    def plotInteractive(self):
        interact(self.combineComputeOpticsplotSimple, 
                 k1 = widgets.FloatSlider(value=-5,description = "H01-EMQ1", min=-50, max=50, step=0.02),
                 k2 = widgets.FloatSlider(value=5, description = "H01-EMQ2", min=-50, max=50, step=0.02),
                 k3 = widgets.FloatSlider(value=-5,description = "H01-EMQ3", min=-50, max=50, step=0.02),
                 k4 = widgets.FloatSlider(value=5, description = "H01-EMQ4", min=-50, max=50, step=0.02),
                 k5 = widgets.FloatSlider(value=-5,description = "H01-EMQ5", min=-50, max=50, step=0.02))
        
        
        
        

