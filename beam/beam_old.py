import scipy.constants
import math
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt

# This class
class Beam(object):
    def __init__(self, type_particle="PROTON"):

        # A particle type
        if type_particle == "ELECTRON":
            self.E_0 = scipy.constants.physical_constants["electron mass energy equivalent in MeV"][0]
            self.q = -1.0
        elif type_particle == "PROTON":
            self.E_0 = scipy.constants.physical_constants["proton mass energy equivalent in MeV"][0]
            self.q = 1.0

        # Attributes

        self.E_total_av = 0 # Average total energy
        self.w_av = 0 # average kinetic energy of the distribution
        self.rel_gamma_ref = 0 # Reference relativistic gamma
        self.rel_beta_ref = 0  # Reference relativistic beta
        self.p = 0 # Momentum of the beam       np.sqrt(self.E_total * self.E_total - self.E_0 * self.E_0)

        self.distribution = [] # Standard Beam distribution
        self.n_particles = 0

        self.x = 0.0
        self.y = 0.0
        self.xp = 0.0
        self.yp = 0.0
        self.phi_rad = 0.0
        self.phi_deg = 0
        self.delta = 0.0
        self.w = 0.0
        self.z = 0.0
        self.t = 0.0
        self.Eptc = 0.0

        # not clean

        #self.E_total = (self.kinetic + self.E_0)
        #self.gamma_r = self.E_total / self.E_0
        #self.beta_r = math.sqrt(1.0 - 1 / (self.gamma_r * self.gamma_r))
        #self.pc = np.sqrt(self.E_total * self.E_total - self.E_0 * self.E_0)

        self.PTC_distribution = []



        self.sigma_x = 0.0
        self.sigma_y = 0.0
        self.sigma_xp = 0.0
        self.sigma_yp = 0.0
        self.sigma_phi = 0.0
        self.sigma_z = 0.0
        self.sigma_deltap = 0.0
        #
        self.x_av = 0
        self.y_av = 0
        self.xp_av = 0
        self.yp_av = 0
        self.phi_av = 0
        self.deltap_av = 0
        self.w_av = 0
        self.z_av = 0
        #
        self.e_x = 0.0
        self.e_y = 0.0
        # self.el = 0.0
        #
        # # Twiss parameters at that location
        self.matrix_4d = []
        self.bx = 0.0
        self.ax = 0.0
        self.gx = 0.0
        self.by = 0.0
        self.ay = 0.0
        self.gy = 0.0
        #
        self.xell = []
        self.xpell = []
        self.yell = []
        self.ypell = []

        self.PTC_input_distribution = []

    # Method to read the standard distribution
    def read_standard_distribution(self, filename):
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
        self.z_av = np.mean(self.z)
        self.E_total_av = self.w_av + self.E_0
        self.rel_gamma_ref = self.E_total_av/self.E_0
        self.delta = 1.0 / (1.0 + 1.0 / self.rel_gamma_ref) * (self.w - self.w_av) / self.w_av
        return self.distribution
    #
    def read_PTC_file_distribution(self, filename):
        self.PTC_distribution = pd.read_csv(filename, skiprows=1, delim_whitespace=True, names=["NUMBER", "TURN", "X", "PX", "Y", "PY", "T", "PT", "S", "E"])
        self.x = self.PTC_distribution["X"]
        self.y = self.PTC_distribution["Y"]
        self.xp = self.PTC_distribution["PX"]
        self.yp = self.PTC_distribution["PY"]
        self.z = -self.PTC_distribution["T"] * self.beta_r
        self.deltap = self.PTC_distribution["PT"] * self.pc / self.E_total / (self.beta_r * self.beta_r)
        self.w = self.deltap * (1 + 1.0 / self.gamma_r) * self.kinetic + self.kinetic
        self.distribution = pd.DataFrame({"X": self.PTC_distribution["X"], "XP": self.PTC_distribution["PX"], "Y": self.PTC_distribution["Y"], "YP": self.PTC_distribution["PY"], "Z": self.z, "DELTAP": self.deltap})
        self.x_av = np.mean(self.x)
        self.y_av = np.mean(self.y)
        self.xp_av = np.mean(self.xp)
        self.yp_av = np.mean(self.yp)
        self.w_av = np.mean(self.w)
        self.z_av = np.mean(self.z)
        return self.distribution


    def get_beam_distribution(self, filename):
        lheader = ["X(mm)", "XP(mrad)", "Y(mm)", "YP(mrad)", "Z(mm)", "W(MeV)"]
        df = pd.read_csv(filename, dtype=float, delim_whitespace=True, header=None, names=lheader, skiprows=1)
        self.distribution = self.get_distribution(df)
        return "Distribution assigned to beam object"

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
        self.z_av = np.mean(self.z)
        return self.distribution

    def build_PTC_distribution(self):
        self.w_av = np.mean(self.w)
        self.Eptc = (self.w - self.w_av) / self.pc
        self.t = self.z / self.beta_r
        self.PTC_input_distribution = pd.DataFrame({"X": self.distribution["X(mm)"] * 0.001, "XP": self.distribution["XP(mrad)"] * 0.001, "Y": self.distribution["Y(mm)"] * 0.001, "YP": self.distribution["YP(mrad)"] * 0.001, "T": self.distribution["Z(mm)"] / self.beta_r * 0.001, "PT": self.Eptc})
        return self.PTC_input_distribution

    def emittance_calculation(self):
        x_var = 0.0
        xxp_var = 0.0
        xp_var = 0.0
        y_var = 0.0
        yyp_var = 0.0
        yp_var = 0.0
        xy_var = 0.0
        xyp_var = 0.0
        xyp_var = 0.0
        xpyp_var = 0.0

        x_var = np.mean(self.distribution["X(mm)"] * self.distribution["X(mm)"])
        xxp_var = np.mean(self.distribution["X"] * self.distribution["XP"])
        xp_var = np.mean(self.distribution["XP"] * self.distribution["XP"])
        y_var = np.mean(self.distribution["Y"] * self.distribution["Y"])
        yyp_var = np.mean(self.distribution["Y"] * self.distribution["YP"])
        yp_var = np.mean(self.distribution["YP"] * self.distribution["YP"])

        xy_var = np.mean(self.distribution["X"] * self.distribution["Y"])
        xyp_var = np.mean(self.distribution["X"] * self.distribution["YP"])
        xpy_var = np.mean(self.distribution["XP"] * self.distribution["Y"])
        xpyp_var = np.mean(self.distribution["YP"] * self.distribution["YP"])

        self.e_x = math.sqrt(x_var * xp_var - xxp_var * xxp_var)
        self.e_y = math.sqrt(y_var * yp_var - yyp_var * yyp_var)
        self.matrix_4d = np.array([[x_var, xxp_var, xy_var, xyp_var], [xxp_var, xp_var, xpy_var, xpyp_var], [xy_var, xpy_var, y_var, yyp_var], [xyp_var, xpyp_var, yyp_var, yp_var]])

        self.bx = x_var / self.e_x
        self.ax = -xxp_var / self.e_x
        self.gx = xp_var / self.e_x
        self.by = y_var / self.e_y
        self.ay = -yyp_var / self.e_y
        self.gy = yp_var / self.e_y

    def compute_standard_deviation_distribution(self):
        self.sigma_x = np.std(self.x)
        self.sigma_y = np.std(self.y)
        self.sigma_xp = np.std(self.xp)
        self.sigma_yp = np.std(self.yp)
        self.sigma_z = np.std(self.z)
        self.sigma_deltap = np.std(self.deltap)

    def compute_average_distribution(self):
        self.x_av = np.mean(self.x)
        self.y_av = np.mean(self.y)
        self.xp_av = np.mean(self.xp)
        self.yp_av = np.mean(self.yp)
        self.z_av = np.mean(self.z)
        self.deltap_av = np.mean(self.deltap)

    def emittance_calculation(self):
        x_var = 0.0
        x_xp_var = 0.0
        xp_var = 0.0
        y_var = 0.0
        y_yp_var = 0.0
        yp_var = 0.0
        x_y_var = 0.0
        x_yp_var = 0.0
        xp_y_var = 0.0
        xp_yp_var = 0.0

        x_var = np.mean(self.distribution["X(mm)"] * self.distribution["X(mm)"])
        x_xp_var = np.mean(self.distribution["X(mm)"] * self.distribution["XP(mrad)"])
        xp_var = np.mean(self.distribution["XP(mrad)"] * self.distribution["XP(mrad)"])
        y_var = np.mean(self.distribution["Y(mm)"] * self.distribution["Y(mm)"])
        y_yp_var = np.mean(self.distribution["Y(mm)"] * self.distribution["YP(mrad)"])
        yp_var = np.mean(self.distribution["YP(mrad)"] * self.distribution["YP(mrad)"])

        x_y_var = np.mean(self.distribution["X(mm)"] * self.distribution["Y(mm)"])
        x_yp_var = np.mean(self.distribution["X(mm)"] * self.distribution["YP(mrad)"])
        xp_y_var = np.mean(self.distribution["XP(mrad)"] * self.distribution["Y(mm)"])
        xp_yp_var = np.mean(self.distribution["XP(mrad)"] * self.distribution["YP(mrad)"])

        self.e_x = math.sqrt(x_var * xp_var - x_xp_var * x_xp_var)
        self.e_y = math.sqrt(y_var * yp_var - y_yp_var * y_yp_var)
        self.matrix_4d = np.array([[x_var, x_xp_var, x_y_var, x_yp_var], [x_xp_var, xp_var, x_yp_var, xp_yp_var], [x_y_var, xp_y_var, y_var, y_yp_var], [x_yp_var, xp_yp_var, y_yp_var, yp_var]])

        self.bx = x_var / self.e_x
        self.ax = -x_xp_var / self.e_x
        self.gx = xp_var / self.e_x
        self.by = y_var / self.e_y
        self.ay = -y_yp_var / self.e_y
        self.gy = yp_var / self.e_y

    def emittance_ellipse(self):
        j = np.arange(100)
        thdeg = j * 13
        teta = thdeg * 2.0 * 3.14 / 360.0
        self.xell = np.sqrt(self.e_x * self.bx) * np.cos(teta)
        self.xpell = -np.sqrt(self.e_x / self.bx) * (self.ax * np.cos(teta) + np.sin(teta))
        self.yell = np.sqrt(self.e_y * self.by) * np.cos(teta)
        self.ypell = -np.sqrt(self.e_y / self.by) * (self.ay * np.cos(teta) + np.sin(teta))

    def select_particle_circle(self, cut):
        return self.distribution[np.sqrt(self.distribution["X(mm)"] * self.distribution["X(mm)"] + self.distribution["Y(mm)"] * self.distribution["Y(mm)"]) <= cut]

        #   def hist(self, data, b):
        #       fig = plt.figure()
        #       plt.hist(data, bins = b)
        #       plt.grid()
        #       return fig

    def phase_space_plot(self, convert_in_mm, xmin, xmax, ymin, ymax, xpmin, xpmax, ypmin, ypmax):
        nbinsx = 500
        nbinsy = 500
        nbinspx = 500
        nbinspy = 1000
        # nbinsz = 100
        # nbinsdp = 100

        # return len(self.x), len(self.y)

        H_xy, xy_xedges, xy_yedges = np.histogram2d(self.x * convert_in_mm, self.y * convert_in_mm, bins=[nbinsx, nbinsy])
        H_xpx, xpx_xedges, xpx_pxedges = np.histogram2d(self.x * convert_in_mm, self.xp * convert_in_mm, bins=[nbinsx, nbinspx])
        H_ypy, ypy_yedges, ypy_pyedges = np.histogram2d(self.y * convert_in_mm, self.yp * convert_in_mm, bins=[nbinsy, nbinspy])
        # H_zdp, zdp_zedges, zdp_dpedges = np.histogram2d(self.z, self.deltap, bins=[nbinsz, nbinsdp])

        H_xy = np.rot90(H_xy)
        H_xy = np.flipud(H_xy)
        H_xpx = np.rot90(H_xpx)
        H_xpx = np.flipud(H_xpx)
        H_ypy = np.rot90(H_ypy)
        H_ypy = np.flipud(H_ypy)
        # H_zdp = np.rot90(H_zdp)
        # H_zdp = np.flipud(H_zdp)
        Hmasked_xy = np.ma.masked_where(H_xy == 0, H_xy)
        Hmasked_xpx = np.ma.masked_where(H_xpx == 0, H_xpx)
        Hmasked_ypy = np.ma.masked_where(H_ypy == 0, H_ypy)
        # Hmasked_zdp = np.ma.masked_where(H_zdp == 0, H_zdp)

        xy_xcenters = xy_xedges[:-1] + 0.5 * (xy_xedges[1:] - xy_xedges[:-1])
        xy_ycenters = xy_yedges[:-1] + 0.5 * (xy_yedges[1:] - xy_yedges[:-1])

        xpx_xcenters = xpx_xedges[:-1] + 0.5 * (xpx_xedges[1:] - xpx_xedges[:-1])
        xpx_pxcenters = xpx_pxedges[:-1] + 0.5 * (xpx_pxedges[1:] - xpx_pxedges[:-1])

        ypy_ycenters = ypy_yedges[:-1] + 0.5 * (ypy_yedges[1:] - ypy_yedges[:-1])
        ypy_pycenters = ypy_pyedges[:-1] + 0.5 * (ypy_pyedges[1:] - ypy_pyedges[:-1])

        # zdp_zcenters = zdp_zedges[:-1] + 0.5 * (zdp_zedges[1:] - zdp_zedges[:-1])
        # zdp_dpcenters = zdp_dpedges[:-1] + 0.5 * (zdp_dpedges[1:] - zdp_dpedges[:-1])

        f_xy = interpolate.interp2d(xy_xcenters, xy_ycenters, Hmasked_xy, kind='cubic')
        f_xpx = interpolate.interp2d(xpx_xcenters, xpx_pxcenters, Hmasked_xpx, kind='cubic')
        f_ypy = interpolate.interp2d(ypy_ycenters, ypy_pycenters, Hmasked_ypy, kind='cubic')
        # f_zdp = interpolate.interp2d(zdp_zcenters, zdp_dpcenters, Hmasked_zdp, kind='cubic')

        axes = []

        dens_xy = [f_xy(i, j) for i, j in zip(self.x * convert_in_mm, self.y * convert_in_mm)]
        dens_xpx = [f_xpx(i, j) for i, j in zip(self.x * convert_in_mm, self.xp * convert_in_mm)]
        dens_ypy = [f_ypy(i, j) for i, j in zip(self.y * convert_in_mm, self.yp * convert_in_mm)]
        # dens_zdp = [f_zdp(i, j) for i, j in zip(self.z, self.deltap)]

        fig = plt.figure(figsize=(30, 30))

        grid = plt.GridSpec(4, 4)
        main_ax_1 = fig.add_subplot(grid[1, 1])

        x_hist = fig.add_subplot(grid[1, 1], sharex=main_ax_1)
        x_hist.axis('off')  # no axis for x histogram
        y_hist = fig.add_subplot(grid[1, 1], sharey=main_ax_1)
        y_hist.axis('off')  # no axis for y histogram

        main_ax_1.scatter(self.x * convert_in_mm, self.y * convert_in_mm, c=np.array(dens_xy).reshape(len(self.x)), s=0.2, cmap=plt.cm.jet, alpha=0.4)
        main_ax_1.grid(alpha=0.5)

        main_ax_1.set_xlim(xmin, xmax)
        main_ax_1.set_ylim(ymin, ymax)
        main_ax_1.set_xlabel("x(mm)")
        main_ax_1.set_ylabel("y(mm)")
        # main_ax_1.text(-14,13,"sig_x="+str('%.3f' %self.sigma_x*convert_in_mm)+"mm",fontsize=12)
        # main_ax_1.text(-14,11, "sig_y="+str('%.3f' %self.sigma_y*convert_in_mm)+"mm", fontsize=12)

        x_hist.hist(self.x * convert_in_mm, nbinsx, histtype='step', orientation='vertical', alpha=0.3, color='b', edgecolor='k');
        y_hist.hist(self.y * convert_in_mm, nbinsy, histtype='step', orientation='horizontal', alpha=0.3, color='b', edgecolor='k');
        # y_hist.set_xlim(0,2000)
        # x_hist.set_ylim(0,2000)

        # --------------------------------------------------------------------------------------------------------------
        main_ax_2 = fig.add_subplot(grid[1, 2])

        xpx_xhist = fig.add_subplot(grid[1, 2], sharex=main_ax_2)
        xpx_xhist.axis('off')  # no axis for x histogram
        xpx_pxhist = fig.add_subplot(grid[1, 2], sharey=main_ax_2)
        xpx_pxhist.axis('off')  # no axis for y histogram

        main_ax_2.scatter(self.x * convert_in_mm, self.xp * convert_in_mm, c=np.array(dens_xpx).reshape(len(self.x)), s=0.2, cmap=plt.cm.jet, alpha=0.4)
        main_ax_2.grid(alpha=0.5)
        # main_ax_2.plot(self.xell, self.xpell, "-r")

        main_ax_2.set_xlim(xmin, xmax)
        main_ax_2.set_ylim(xpmin, xpmax)
        main_ax_2.set_xlabel("x(mm)")
        main_ax_2.set_ylabel("x'(mrad)")

        xpx_xhist.hist(self.x * convert_in_mm, nbinsx, histtype='step', orientation='vertical', alpha=0.3, color='b', edgecolor='k');
        xpx_pxhist.hist(self.xp * convert_in_mm, nbinspx, histtype='step', orientation='horizontal', alpha=0.3, color='b', edgecolor='k');

        # ----------------------------------------------------------------------------------------------------------------
        main_ax_3 = fig.add_subplot(grid[2, 1])
        ypy_yhist = fig.add_subplot(grid[2, 1], sharex=main_ax_3)
        ypy_yhist.axis('off')  # no axis for x histogram
        ypy_pyhist = fig.add_subplot(grid[2, 1], sharey=main_ax_3)
        ypy_pyhist.axis('off')  # no axis for y historgram

        main_ax_3.scatter(self.y * convert_in_mm, self.yp * convert_in_mm, c=np.array(dens_ypy).reshape(len(self.y)), s=0.2, cmap=plt.cm.jet, alpha=0.4)
        main_ax_3.grid(alpha=0.5)
        # main_ax_3.plot(self.yell , self.ypell, "-r")

        main_ax_3.set_xlim(ymin, ymax)
        main_ax_3.set_ylim(ypmin, ypmax)
        main_ax_3.set_xlabel("y(mm)")
        main_ax_3.set_ylabel("y'(mrad)")

        ypy_yhist.hist(self.y * convert_in_mm, nbinsy, histtype='step', orientation='vertical', alpha=0.3, color='b', edgecolor='k');
        ypy_pyhist.hist(self.yp * convert_in_mm, nbinspy, histtype='step', orientation='horizontal', alpha=0.3, color='b', edgecolor='k');
        # xpx_xhist.set_xlim(0,2000)
        # xpx_pxhist.set_ylim(0,2000)

        # ----------------------------------------------------------------------------------------------------------------
        # main_ax_4 = fig.add_subplot(grid[2, 2])
        #
        # zdp_zhist = fig.add_subplot(grid[2, 2], sharex=main_ax_4)
        # zdp_zhist.axis('off')  # no axis for x histogram
        # zdp_dphist = fig.add_subplot(grid[2, 2], sharey=main_ax_4)
        # zdp_dphist.axis('off')  # no axis for y historgram
        #
        # main_ax_4.scatter(self.z, self.deltap, c=np.array(dens_zdp).reshape(len(self.z)), s=0.2, cmap=plt.cm.jet, alpha=0.4)
        # main_ax_4.grid(alpha=0.5)
        #
        # main_ax_4.set_xlim(zmin, zmax)
        # main_ax_4.set_ylim(dpmin, dpmax)
        # main_ax_4.set_xlabel("z(m)")
        # main_ax_4.set_ylabel("dp/p")
        #
        # zdp_zhist.hist(self.z, nbinsz, histtype='step', orientation='vertical', alpha=0.3, color='b', edgecolor='k');
        # zdp_dphist.hist(self.deltap, nbinsdp, histtype='step', orientation='horizontal', alpha=0.3, color='b', edgecolor='k');
        return fig

