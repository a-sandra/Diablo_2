import scipy.constants
from scipy import interpolate
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time


# Class Beam 
## 10/02/2023: I added an argument in the constructor, the input beam distrubution. 
class Beam(object):
    def __init__(self, distribution_file_name="", type_particle="PROTON" ):

        # A particle type
        if type_particle == "ELECTRON":
            self.E_0 = scipy.constants.physical_constants["electron mass energy equivalent in MeV"][0]
            self.q = -1.0
        elif type_particle == "PROTON":
            self.E_0 = scipy.constants.physical_constants["proton mass energy equivalent in MeV"][0]
            self.q = 1.0

        # 10/02/2023 I added a condition: if the input file is empty, then a message notifying the user is printing. 
        ## if the file is giving, then the method self.read_distribution is used.
        if not distribution_file_name: 
            self.distribution = []
            self.n_particle = 0
            self._p_mom_av = 0
            self._e_total_av = 0
            self.ptc_input_distribution = []
            self.ptc_input_distribution_cut = [] # ?
            self.t = 0
            self.E_ptc = 0
            self.ptc_coordinates = []
            self.travel_distribution = []
            print("beam created without any input beam distribution. Use read_distribution(filename) et voila")
        else: 
            self.read_distribution(distribution_file_name)
            

    @property
    def p_mom_av(self):
        return self._p_mom_av

    def p_mom_av(self):
        self._p_mom_av = math.sqrt((self.w_av + self.E_0) ** 2 - self.E_0 ** 2)
        return self._p_mom_av

    def e_total_av(self):
        self._e_total_av = self.E_0 + self.w_av

    def lorentz_relativist_factor(self):
        self.gamma_r_av = self._e_total_av / self.E_0
        self.beta_r_av = math.sqrt(1 - 1 / self.gamma_r_av ** 2)
        pass

    def get_energy_parameters(self):
        self.p_mom_av()
        self.e_total_av()
        self.lorentz_relativist_factor()
        pass

    # This method reads a beam distribution with a standard format ["X(mm)", "XP(mrad)", "Y(mm)", "YP(mrad)", "Z(mm)", "W(MeV)"] in mm, mrad, MeV
    def read_distribution(self, filename):
        print("Hello")
        try:
            print("je suis la")
            lheader = ["X(mm)", "XP(mrad)", "Y(mm)", "YP(mrad)", "Z(mm)", "W(MeV)"]
            self.distribution = pd.read_csv(filename, delim_whitespace=True, skiprows=1, names=lheader)
            self.x = self.distribution["X(mm)"]
            self.y = self.distribution["Y(mm)"]
            self.xp = self.distribution["XP(mrad)"]
            self.yp = self.distribution["YP(mrad)"]
            self.z = self.distribution["Z(mm)"]
            self.w = self.distribution["W(MeV)"]
            self.x_av = np.mean(self.x)
            self.y_av = np.mean(self.y)
            self.xp_av = np.mean(self.xp)
            self.yp_av = np.mean(self.yp)
            self.w_av = np.mean(self.w)
            self.w_rms = np.std(self.w)
            self.w_med = np.median(self.w)
            self.z_av = np.mean(self.z)
            self.get_energy_parameters()
            self.dpp = (self.w + self.E_0 - self._e_total_av) / self._e_total_av / (
                        self.beta_r_av ** 2)
            self.n_particle = len(self.w)
            return "The beam distribution has been loaded successfully , You should check the distribution is correctly formatted into the dataframe - type thenameofyourobject.distribution"
        except IOError:
            print("Ooops, the file is not found ¯\_(ツ)_/¯")
        finally:
            pass

    def get_distribution(self, df):
        self.distribution = df
        self.x = self.distribution["X(mm)"]
        self.y = self.distribution["Y(mm)"]
        self.xp = self.distribution["XP(mrad)"]
        self.yp = self.distribution["YP(mrad)"]
        self.w = self.distribution["W(MeV)"]
        self.z = self.distribution["Z(mm)"]
        self.x_av = np.mean(self.x)
        self.y_av = np.mean(self.y)
        self.xp_av = np.mean(self.xp)
        self.yp_av = np.mean(self.yp)
        self.w_av = np.mean(self.w)
        self.w_med = np.median(self.w)
        self.z_av = np.mean(self.z)
        self.n_particle = len(self.w)
        self.get_energy_parameters()
        return self.distribution

    def build_ptc_distribution(self):
        self.E_ptc = (self.w - self.w_av) / self._p_mom_av
        self.t = -self.z / self.beta_r_av
        self.ptc_input_distribution = pd.DataFrame(
            {"X": self.distribution["X(mm)"] * 0.001, "XP": self.distribution["XP(mrad)"] * 0.001,
             "Y": self.distribution["Y(mm)"] * 0.001, "YP": self.distribution["YP(mrad)"] * 0.001, "T": self.t * 0.001,
             "PT": self.E_ptc})
        return self.ptc_input_distribution

    def write_distribution(self, prefix="distribution", suffix="madx"):
        print_coord = open(prefix+"_" + str(round(self.w_av, 0)) + "." + str(suffix), "w")
        for line in self.ptc_coordinates:
            print_coord.write(line)
            print_coord.write("\n")
        print_coord.close()
        return print("print_coord acquired")

    def standard2ptc(self, prefix, suffix, n=20000):
        self.build_ptc_distribution()
        if n <= len(self.ptc_input_distribution["X"]):
            self.ptc_input_distribution_cut = self.ptc_input_distribution.sample(n).reset_index()
            self.ptc_coordinates = ["ptc_start, x=" + str(self.ptc_input_distribution_cut["X"][i]) + ", px=" + str(
                self.ptc_input_distribution_cut["XP"][i]) + ", y=" + str(
                self.ptc_input_distribution_cut["Y"][i]) + ", py=" + str(
                self.ptc_input_distribution_cut["YP"][i]) + ", t=" + str(
                self.ptc_input_distribution_cut["T"][i]) + ", pt=" + str(self.ptc_input_distribution_cut["PT"][i]) + ";"
                                    for i in np.arange(n)]
            self.write_distribution(prefix, str(suffix))
        else:
            print("The number of particles to sample is larger than the length of the distribution. By default n =20k")
        return "Translated"  # "Save the input PTC distribution at " + str(round(self.w_av, 0))

    def get_ptc_distribution(self, filename):
        self.PTC_distribution = pd.read_csv(filename, skiprows=1, delim_whitespace=True, names=["NUMBER", "TURN", "X", "PX", "Y", "PY", "T", "PT", "S", "E"])
        self.x = self.PTC_distribution["X"]*1000
        self.y = self.PTC_distribution["Y"]*1000
        self.xp = self.PTC_distribution["PX"]*1000
        self.yp = self.PTC_distribution["PY"]*1000
        self.t = -self.PTC_distribution["T"]
        self.E_ptc = self.PTC_distribution["PT"]
        self.E = self.PTC_distribution["E"]
        self._e_total_0 = np.mean(self.E)*1000 # in MeV
        self.w_0 = self._e_total_0 - self.E_0
        self.p_0 = math.sqrt(self._e_total_0**2 - self.E_0**2)
        self.w_spread = self.E_ptc*self.p_0 # in MeV
        self.w = self.E_ptc * self.p_0 + self.w_0
        self.gamma_r_0 = self._e_total_0 / self.E_0
        self.beta_r_0 = math.sqrt(1 - 1 / self.gamma_r_0 ** 2)
        self.z = -self.PTC_distribution["T"] * self.beta_r_0*1000
        self.distribution = pd.DataFrame({"X(mm)": self.x, "XP(mrad)": self.xp,
                                          'Y(mm)': self.y, "YP(mrad)": self.yp,
                                          "Z(mm)": self.z, "W(MeV)": self.w})
        self.x_av = np.mean(self.x)
        self.y_av = np.mean(self.y)
        self.xp_av = np.mean(self.xp)
        self.yp_av = np.mean(self.yp)
        self.w_av = np.mean(self.w)
        self.z_av = np.mean(self.z)
        return self.distribution

    def standard2travel(self, filename , phi_ref_rad = 0, rf_frequency=3.0e9):
        len_distribution = len(self.x)
        f = open(filename, "w")
        f.write("TRAVEL BEAM DATA FILE\n")
        f.write("Date created: " + time.strftime("%x") + "\n")
        f.write(str(self._p_mom_av*0.001) + "    !REFERENCE MOMENTUM [GeV/c]\n")
        f.write(str(phi_ref_rad) + "  !REFERENCE PHASE [rad]\n")
        f.write(str(rf_frequency) + "    !FREQUENCY [Hz]\n")
        f.write(str(self.E_0*0.001) + "    !REFERENCE MASS [GeV/c2]\n")
        f.write(str(self.q) + "       !REFERENCE CHARGE STATE\n")
        f.write(str(len_distribution) + "    !NUMBER OF PARTICLES\n")

        self.travel_distribution = np.zeros((len_distribution, 10))
        self.travel_distribution[:, 0] = np.arange(0, len_distribution) + 1
        self.travel_distribution[:, 1] = self.x * 0.001
        self.travel_distribution[:, 2] = self.xp * 0.001
        self.travel_distribution[:, 3] = self.y * 0.001
        self.travel_distribution[:, 4] = self.yp * 0.001
        self.travel_distribution[:, 5] = self.z * 0.001 # not correct
        self.travel_distribution[:, 6] = (self.w + self.E_0 - self._e_total_av) / self._e_total_av / (self.beta_r_av ** 2)
        self.travel_distribution[:, 8] = np.ones(len_distribution)
        self.travel_distribution[:, 9] = np.ones(len_distribution) * self.E_0*0.001  # convert from MeV to GeV
        for line in range(0, len_distribution):
            f.write("\t".join(str(i) for i in self.travel_distribution[line, :]) + "\n")
        f.close()

    #def read_TRAVEL_distribution(filename):
        #df = pd.read_csv(filename, skiprows=8, sep="    ",
         #                names=["NB", "X", "XP", "Y", "YP", "PHASE", "DP/P", "Q", "LOSS", "E0"], engine='python')
        #return df
    # ---------------------------------------------------------
    # EMITTANCE
    def emittance_calculation(self):
        x_var = np.mean(self.distribution["X(mm)"] * self.distribution["X(mm)"])
        xxp_var = np.mean(self.distribution["X(mm)"] * self.distribution["XP(mrad)"])
        xp_var = np.mean(self.distribution["XP(mrad)"] * self.distribution["XP(mrad)"])
        y_var = np.mean(self.distribution["Y(mm)"] * self.distribution["Y(mm)"])
        yyp_var = np.mean(self.distribution["Y(mm)"] * self.distribution["YP(mrad)"])
        yp_var = np.mean(self.distribution["YP(mrad)"] * self.distribution["YP(mrad)"])

        xy_var = np.mean(self.distribution["X(mm)"] * self.distribution["Y(mm)"])
        xyp_var = np.mean(self.distribution["X(mm)"] * self.distribution["YP(mrad)"])
        xpy_var = np.mean(self.distribution["XP(mrad)"] * self.distribution["Y(mm)"])
        xpyp_var = np.mean(self.distribution["YP(mrad)"] * self.distribution["YP(mrad)"])

        self.e_x = math.sqrt(x_var * xp_var - xxp_var * xxp_var)
        self.e_y = math.sqrt(y_var * yp_var - yyp_var * yyp_var)
        self.matrix_4d = np.array(
            [[x_var, xxp_var, xy_var, xyp_var], [xxp_var, xp_var, xpy_var, xpyp_var], [xy_var, xpy_var, y_var, yyp_var],
             [xyp_var, xpyp_var, yyp_var, yp_var]])

        self.bx = x_var / self.e_x
        self.ax = -xxp_var / self.e_x
        self.gx = xp_var / self.e_x
        self.by = y_var / self.e_y
        self.ay = -yyp_var / self.e_y
        self.gy = yp_var / self.e_y

    def emittance_ellipse(self):
        j = np.arange(100)
        thdeg = j * 13
        teta = thdeg * 2.0 * 3.14 / 360.0
        self.xell = np.sqrt(self.e_x * self.bx) * np.cos(teta)
        self.xpell = -np.sqrt(self.e_x / self.bx) * (self.ax * np.cos(teta) + np.sin(teta))
        self.yell = np.sqrt(self.e_y * self.by) * np.cos(teta)
        self.ypell = -np.sqrt(self.e_y / self.by) * (self.ay * np.cos(teta) + np.sin(teta))

    # ---------------------------------------------------------
    # PLOTTING

    def hist(self, column, xlab, ylab, legend, place="upper left"):
        n_occ, x_bin, list_obj = plt.hist(column, bins='auto', label=legend)
        plt.grid(linestyle="--")
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.legend(loc=place)

    def phase_space_hist(self, size=(12, 5), legend="upper left"):
        plt.figure(figsize=size)
        plt.subplot(2, 2, 1)
        self.hist(self.x, "X(mm)", "# Occurence", "X", legend)
        plt.subplot(2, 2, 2)
        self.hist(self.xp, "XP(mrad)", "# Occurence", "XP", legend)
        plt.subplot(2, 2, 3)
        self.hist(self.y, "Y(mm)", "# Occurence", "Y", legend)
        plt.subplot(2, 2, 4)
        self.hist(self.yp, "YP(mrad)", "# Occurence", "YP", legend)

    # PLOTTING

    def single_2D_phase_plot(self, convert_in_mm, xdata, ydata, xmin, xmax, ymin, ymax, lx, ly):
        nbinsx = 250
        nbinsy = 250
        H_xy, xy_xedges, xy_yedges = np.histogram2d(xdata* convert_in_mm, ydata * convert_in_mm, bins=[nbinsx, nbinsy])
        H_xy = np.rot90(H_xy)
        H_xy = np.flipud(H_xy)
        Hmasked_xy = np.ma.masked_where(H_xy == 0, H_xy)

        xy_xcenters = xy_xedges[:-1] + 0.5 * (xy_xedges[1:] - xy_xedges[:-1])
        xy_ycenters = xy_yedges[:-1] + 0.5 * (xy_yedges[1:] - xy_yedges[:-1])

        f_xy = interpolate.interp2d(xy_xcenters, xy_ycenters, Hmasked_xy, kind='cubic')
        axes = []

        dens_xy = [f_xy(i, j) for i, j in zip(xdata * convert_in_mm, ydata * convert_in_mm)]

        fig = plt.figure(figsize=(30, 30))

        grid = plt.GridSpec(4, 4)
        main_ax_1 = fig.add_subplot(grid[1, 1])

        x_hist = fig.add_subplot(grid[1, 1], sharex=main_ax_1)
        x_hist.axis('off')  # no axis for x histogram
        y_hist = fig.add_subplot(grid[1, 1], sharey=main_ax_1)
        y_hist.axis('off')  # no axis for y histogram

        main_ax_1.scatter(xdata * convert_in_mm, ydata * convert_in_mm, c=np.array(dens_xy).reshape(len(xdata)), s=0.2, cmap=plt.cm.jet, alpha=0.4)
        main_ax_1.grid(alpha=0.5)

        main_ax_1.set_xlim(xmin, xmax)
        main_ax_1.set_ylim(ymin, ymax)
        main_ax_1.set_xlabel(lx)
        main_ax_1.set_ylabel(ly)
        # main_ax_1.text(-14,13,"sig_x="+str('%.3f' %self.sigma_x*convert_in_mm)+"mm",fontsize=12)# main_ax_1.text(-14,11, "sig_y="+str('%.3f' %self.sigma_y*convert_in_mm)+"mm", fontsize=12)

        x_hist.hist(xdata * convert_in_mm, nbinsx, histtype='step', orientation='vertical', alpha=0.3, color='b', edgecolor='k');
        y_hist.hist(ydata * convert_in_mm, nbinsy, histtype='step', orientation='horizontal', alpha=0.3, color='b', edgecolor='k');
        # y_hist.set_xlim(0,2000)
        # x_hist.set_ylim(0,2000)