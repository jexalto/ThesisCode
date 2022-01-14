""" Extension of Q_prop_wing_arb_1W
 2-way interaction
 inboard and tip mounted prop
 different BEM instance
  Copy of config3 (version 3) but rewritten as function"""

import copy
import numpy as np
import pandas as pd
import sys

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

if 'Configurations/' not in sys.path:
    sys.path.append('Configurations/')

from BEM import BEM
from bezier import bezier
from create_unsteady_rotor_model import create_unsteady_rotor_model
from vlm import PyVLM

#from TipMountedProp_Prop2 import PropWing
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time
#from config3 import InboardPropeller, WTMP, Wing, FlightConditions, Options


class WingSys(object):
    def __init__(self, fc, ib, wt, wing, options, propellers, plot=False):
        # ----------------- Flight conditions ---------------------
        Vinf = fc.Vinf  # freestream velocity
        a = fc.a  # speed of sound
        rho = fc.rho  # density
        mu = fc.mu  # viscosity

        # ------------- Wing --------------------
        c_r = wing.c_r  # root chord
        c_t = wing.c_t  # tip chrord
        b = wing.b  # span
        wing_ang = wing.alpha

        # ------------ Propellers-----------------
        for prop in propellers:
            prop.loc = np.array([prop.x, prop.y , prop.z])  # inboard prop location
            prop.bem.Vinf = Vinf
            prop.bem.omega = prop.omega
            prop.bem.rho = rho
            prop.bem.a = a
            prop.bem.mu = mu
            prop.alpha_prop = wing_ang + prop.prop_ang
            prop.bem.airfoil_folder = prop.bem.airfoil_folder
            prop.bem.Vinf = Vinf

        if ib not in propellers:
            ib.loc = False
        if wt not in propellers:
            wt.R = False



        # ------------------- VLM -------------------------
        n = wing.n  # chordwise panels per section
        m = wing.m  # spanwise panels per section
        opt = wing.opt  # spacing type per section

        le_1 = le_2 = le_3 = le_4 = le_5 = np.array([0, 0])
        if wing.le_t:
            le_6 = np.array([wing.le_t,b/2])
        else:
            le_6 = np.array([(c_r-c_t)/2,b/2])
        tap = (le_6[0]-le_1[0]) / b * 2  #change in le pos per m

        if ib in propellers:
            le_2 = np.array([0, ib.loc[1] - ib.R])
            le_3 = np.array([0, ib.loc[1] + ib.R])
            le_4 = np.array([0, le_3[1] + le_2[1]])

            le_2[0] = le_2[1] * tap
            le_3[0] = le_3[1] * tap
            le_4[0] = le_4[1] * tap

        if wt in propellers:
            le_5 = np.array([0, wt.loc[1] - wt.R])
            le_5[0] = le_5[1] *tap


        c_2 = c_r - (le_2[1] - le_1[1]) * (c_r - c_t) / b * 2
        c_3 = c_r - (le_3[1] - le_1[1]) * (c_r - c_t) / b * 2
        c_4 = c_r - (le_4[1] - le_1[1]) * (c_r - c_t) / b * 2
        c_5 = c_r - (le_5[1] - le_1[1]) * (c_r - c_t) / b * 2

        le = [le_1, le_2, le_3, le_4, le_5, le_6]  # wing leading edge coordinates
        le_idx = [0] + np.nonzero(le)[0].tolist()
        chord = [c_r, c_2, c_3, c_4, c_5, c_t]  # corresponding chord

        le = [le[i] for i in np.unique(le_idx)]
        chord = [chord[i] for i in np.unique(le_idx)]

        if ib in propellers:
            if not np.isclose(le_2[1], le_4[1] - le_3[1]) or wing.m[0] != wing.m[2] or wing.opt[0] != 'n' + wing.opt[2]:
                raise Exception(
                    'Please check that section 1 has the same characteristics as section 3 and an opossing spacing')
                # Check to see if panel symmetry around inboard propeller

            if m[1]%2==0:
                raise Exception('Use uneven number of panels for inboard propeller section')

        vlm = PyVLM()  # initiate VLM
        airfoil = np.loadtxt(wing.airfoil_dir + wing.airfoil + '.dat', skiprows=1)
        vlm.airfoil.input_selig(airfoil[:, 0], airfoil[:, 1], wing.airfoil)  # airfoil name!
        vlm.add_wing(le, chord, n, m, opt)
        vlm.polar_dir = wing.polar_dir

        vlm.Vinf = Vinf
        vlm.rho = rho
        vlm.a = a
        vlm.mu = mu
        vlm.wt_R = wt.R
        vlm.ib_loc = ib.loc
        vlm.ib_R = ib.R
        vlm.alpha = wing_ang

        # -----------------  Settings --------------------------------------------------
        settings = options.settings
        N_iter_max = options.N_iter_max
        visc = options.visc

        turb = options.turb
        skinf = options.skinf

        # pass settings
        for key in vlm.jet_corr_settings.keys():
            vlm.jet_corr_settings[key] = settings[key]

        self.propellers = propellers
        self.vlm = vlm
        self.settings = settings
        self.N_iter_max = N_iter_max
        self.vlm_clean = copy.deepcopy(vlm)  # copy of clean vlm = clean wing
        self.conv = None
        self.visc = visc
        self.turb = turb
        self.skinf = skinf
        self.le = le
        self.chord = chord
        self.tmp_dct = {}
        self.plot = plot

        if self.plot == True:
            self.plot_geom()


    def analyse(self):
        t0 = time.time()  # get time
        self.analyse_bem()                      # analyse BEM(s)
        if self.conv_bem is True:
            self.add_slipstream()               # add propeller slipstream(s)

            # data frame to store convergence data
            conv = pd.DataFrame(columns=['dCT', 'dCQ', 'dCL', 'dCD'])
            ct0, cq0, cl0, cd0 = 10, 10, 10, 10  # guess
            conv_iter = False

            for N in range(self.N_iter_max):
                self.analyse_prop() # anlayse propeller
                self.prop_induced_velocities()      # introduce propeller induced velocities to vlm
                self.run_vlm()
                self.wing_induced_velocities()
                self.slipstream_deflection()


                cl1 = self.vlm.res['CL']
                cd1 = self.vlm.res['CD']

                # convergence data
                if len(self.propellers) > 0:
                    ct1 = self.propellers[0].urot.prop_us['integral_CT']
                    cq1 = self.propellers[0].urot.prop_us['integral_CQ']

                    conv.loc[N] = [np.abs(ct0 - ct1), np.abs(cq0 - cq1), np.abs(cl0 - cl1), np.abs(cd0 - cd1)]
                else:
                    conv.loc[N] = [0, 0, 0, 0]  # no iter required as no propeller

                self.tmp_dct['urot_1_%d' % N] = copy.deepcopy(self.propellers[0].urot)  # dictionary of ib urot

                if len(self.propellers) >1 :
                    self.tmp_dct['urot_2_%d' % N] = copy.deepcopy(self.propellers[1].urot)  # dictionary of wt urot
                self.tmp_dct['vlm_%d' % N] = copy.deepcopy(self.vlm)

                if conv.loc[N, 'dCT'] < 1e-4 and conv.loc[N, 'dCQ'] < 1e-4 and conv.loc[N, 'dCL'] < 1e-3 and conv.loc[
                    N, 'dCD'] < 1e-4:
                    print(N)
                    conv_iter = True
                    print('iteration done')
                    break
                else:
                    ct0 = ct1
                    cq0 = cq1
                    cl0 = cl1
                    cd0 = cd1
                    print(conv)
                    print(self.vlm.res)
                    print('next iteration')

            if not conv_iter:
                N = conv['dCD'].idxmin()

            self.conv = conv.iloc[:N + 1]
            print('Done')
            self.converged = conv_iter

            print((time.time() - t0) / 60)
            #return vlm_clean, tmp_dct, conv



# --------------- Initiate --------------------------------------------------
    def analyse_bem(self):
        conv_bem = True
        for prop in self.propellers:
            prop.bem.analysis()
            print('bem analysed')
            if not  prop.bem.res_sum['converged']:
                conv_bem = False
      #  self.propellers = propellers
        self.conv_bem = conv_bem


    def add_slipstream(self):
        # create a propeller map
        for prop in self.propellers:
            if prop.urot is None:
                dJ = self.settings['dJ']
                N_prop_map = self.settings['N_prop_map']
                prop.urot = create_unsteady_rotor_model(prop.bem, dJ, N_prop_map)
                print('urot created')
                prop.urot.N_time_step = self.settings['N_time_step']
                prop.urot.sign_rotation = prop.sign_rotation

            else:
                prop.urot = prop.urot

            # slipstream geometry
            sl = self.settings['slipstream_length']  # slipstream length
            dsl = prop.bem.Vinf / (prop.bem.omega / (2 * np.pi)) / self.settings['N_slipstream_x']
            # prevent slipstream with only one x location
            if dsl >= sl:
                dsl = sl / 2.001
            prop.x_st = np.arange(0, sl, dsl)
            # prop.x_st = (1 - np.cos(np.linspace(0, 0.5 * np.pi, 100))) * 2
            # prop.z_st = np.zeros(prop.x_st.shape)
            prop.z_st = np.zeros(prop.x_st.shape)

            # slipstream settings
            prop.urot.slipstream.cnw.ds = self.settings['conway_ds']
            prop.urot.slipstream.cnw.s_max = self.settings['conway_s_max']

            # set propeller inflow field
            prop.u_aoa = np.cos(np.deg2rad(prop.alpha_prop))
            prop.w_aoa = np.sin(np.deg2rad(prop.alpha_prop))
            prop.grid_y = np.linspace(-prop.R, prop.R, prop.grid_N)
            prop.grid_z = np.linspace(-prop.R, prop.R, prop.grid_N)
            prop.u_in = np.ones((prop.grid_y.size, prop.grid_z.size)) * prop.u_aoa
            prop.v_in = np.zeros((prop.grid_y.size, prop.grid_z.size))
            prop.w_in = np.ones((prop.grid_y.size, prop.grid_z.size)) * prop.w_aoa



    def analyse_prop(self):
    # analyse prop
        for prop in self.propellers:  # initialize slipstream
            print('analysing prop started')
            prop.urot.inflow.grid_input(prop.grid_y / prop.R, prop.grid_z / prop.R, prop.u_in, prop.v_in, prop.w_in)
            prop.urot.analysis()

            prop.urot.calculate_circulation(prop.bem.res['Va_i'].to_numpy(),
                                      prop.bem.res['Vt_i'].to_numpy(),
                                      prop.bem.res['rR'].to_numpy())

            prop.urot.init_slipstream(prop.x_st, prop.z_st)

            print('analysing prop done')



    def prop_induced_velocities(self):
        # ------ propeller induced velocities on wing ---------------------------------
        Vpert = np.zeros((self.vlm.Vpert.shape[0], 3))
        Vperti = np.zeros((self.vlm.Vpert.shape[0], 3))
        P_lst = []
        P1_lst = []
        P2_lst = []
        V1x_lst = []
        V2x_lst = []
        V1z_lst=[]
        V2z_lst=[]
        V1x_lst_wt = []
        V2x_lst_wt = []
        V1z_lst_wt=[]
        V2z_lst_wt=[]

        for i in range(self.vlm.Vpert.shape[0]):
            P_ = np.array([self.vlm.Vpert.iloc[i]['x'],
                           self.vlm.Vpert.iloc[i]['y'],
                           0])
            P_lst.append(P_[1])
            i = 0
            for prop in self.propellers:
                i = i+1
                urot = prop.urot
                P1, P2 = wing_to_prop(prop, P_)
                P1_lst.append(P1[1])
                P2_lst.append((P2[1]))

                Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
                Vt2, Vw2 = urot.slipstream.induced_velocity(P2)

                V1 = np.matmul(prop.prop_rot, Vt1)
                V2 = np.matmul(prop.prop_rot, Vt2)

                V2 = V2 * np.array([1, -1, 1])  # left prop

                V1 = V1 *np.array([1,-1*prop.sign_rotation, -1*prop.sign_rotation])
                V2 = V2 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])

                if i ==1:
                    V1x_lst.append(V1[0])
                    V2x_lst.append(V2[0])
                    V1z_lst.append(V1[2])
                    V2z_lst.append(V2[2])
                else:
                    V1x_lst_wt.append(V1[0])
                    V2x_lst_wt.append(V2[0])
                    V1z_lst_wt.append(V1[2])
                    V2z_lst_wt.append(V2[2])

                Vpert[i] = Vpert[i] + V1 + V2  # add induced velocities by both props

                # velocities for trefftz plane
                xi = np.max(urot.slipstream.x)/2
                x0 = P1[0]
                if urot.slipstream.r_r0 is not None:

                    sl_r_r0 = urot.slipstream.r_r0
                    sl_x = urot.slipstream.x.flatten()
                    sl_z = urot.slipstream.z.flatten()
                    i_xi = np.argmin(np.abs(sl_x-xi))
                    i_x0 = np.argmin(np.abs(sl_x-x0))

                    z0 = np.interp(x0, sl_x, sl_z)
                    zi = np.interp(xi, sl_x, sl_z)
                    r_r0 = sl_r_r0[i_x0].flatten()

                    r0 = (P1[1]**2+P1[2]**2)**0.5

                    if r0>np.max(r_r0)*prop.R:
                        r_r = 1
                    else:
                        i_r = np.argmin(np.abs(r0-r_r0*prop.R))
                        r_r = sl_r_r0[i_xi, 0, i_r]/r_r0[i_r]

                    P1[1] = P1[1]*r_r
                    P1[2] = (P1[2]-z0)*r_r+zi
                    P2[1] = P2[1]*r_r
                    P2[2] = (P2[2]-z0)*r_r+zi

                P1[0] = xi
                P2[0] = xi

                Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
                Vt2, Vw2 = urot.slipstream.induced_velocity(P2)

                V1 = np.matmul(prop.prop_rot, Vt1)
                V2 = np.matmul(prop.prop_rot, Vt2)

                V2 = V2*np.array([1, -1, 1]) #left prop
                V1 = V1 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])
                V2 = V2 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])

                Vperti[i] =  Vperti[i] + V1+V2

        self.vlm.Vpert['Vx']   = Vpert[:, 0]
        self.vlm.Vpert['Vz']   = Vpert[:, 2]
        self.vlm.Vpert['Vzi']   = Vperti[:, 2]
        self.V1x_lst_wt = V1x_lst_wt
        self.V2x_lst_wt = V2x_lst_wt
        self.V1z_lst_wt = V1z_lst_wt
        self.V2z_lst_wt = V2z_lst_wt

        self.V1x_lst_ib = V1x_lst
        self.V2x_lst_ib = V2x_lst
        self.V1z_lst_ib = V1z_lst
        self.V2z_lst_ib = V2z_lst



    def run_vlm(self):
        self.vlm.G = 0
        if self.turb:
            self.vlm.Vpert['turb'] = True

        if self.visc:
            self.vlm.vlm_visc()
        else:
            self.vlm.vlm()
            if self.skinf:
                self.vlm.strip_properties(True)
            else:
                self.vlm.strip_properties(False)



    def wing_induced_velocities(self):
        for prop in self.propellers:
            for iy, y_ in enumerate(prop.grid_y):
                for iz, z_ in enumerate(prop.grid_z):
                    P_ = np.array([0, y_, z_])
                    P = prop_to_wing(prop, P_)

                    V = self.vlm.induced_velocity(P)                 # this takes time ~0.15 sec
                    V = np.matmul(prop.prop_rot.T, V)

                    prop.u_in[iy, iz] = V[0]/self.vlm.Vinf+prop.u_aoa
                    prop.v_in[iy, iz] = V[1]/self.vlm.Vinf
                    prop.w_in[iy, iz] = V[2]/self.vlm.Vinf+prop.w_aoa



    def slipstream_deflection(self):
        # slipstream deflection
        for prop in self.propellers:
            for i in range(prop.z_st.size-1):
                P_ = np.array([prop.x_st[i], 0, 0])
                P = prop_to_wing(prop, P_)

                w = self.vlm.induced_velocity(P)
                w = np.matmul(prop.prop_rot.T, w)

                dx = prop.x_st[i+1]-prop.x_st[i]
                prop.z_st[i+1] = prop.z_st[i]+dx*w[2]/self.vlm.Vinf


    # self.urot = urot
    # self.vlm = vlm
    # self.conv = conv

    def plot_geom(self):

        # Plot wing plan form
        plt.figure()
        y_wing = []
        x_wing = []
        xc_wing =[]
        for i in np.arange(len(self.le)):
            le = self.le[i]
            y_wing.append(le[1])
            x_wing.append(le[0])
            xc_wing.append(le[0]+self.chord[i])
        chord = self.chord
        y_wing.extend( y_wing[::-1])
        x_wing.extend(xc_wing[::-1])

        plt.plot(y_wing,x_wing)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.gca().invert_yaxis()
        plt.xlabel('y [m]')
        plt.ylabel('x [m]')

        # Plot front view, normalised
        plt.figure()
        for prop in self.propellers:
            span = self.le[-1][1] *2
            radius = prop.bem.R
            prop_loc_y = prop.loc[1]
            prop_loc_z = prop.loc[2]

            plt.plot([0, span / 2], [0, 0], color='black')
            pr = plt.Circle((prop_loc_y, prop_loc_z), radius, color='r', fill=False)
            plt.gcf().gca().add_artist(pr)

        plt.xlim(0, self.le[-1][1])
        plt.ylim(-.2, 0.2)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('y/b [-]')
        plt.ylabel('z/b [-]')
        plt.show()

def prop_to_wing(propeller, P):
    """Coordinate transform from propeller to wing coordinate system

    Args:
        propeller(): propeller under consideration (ib or wt)
        P (numpy.array): point in propeller coordinate system

    Returns:
        p (numpy.array): point in wing coordinate system
    """

    p = P + propeller.loc
    p = np.matmul(propeller.prop_rot, p)

    return p


def wing_to_prop(propeller, P):
    """Coordinate transform from wing to propellercoordinate system

    Args:
        propeller: propeller under consideration (ib or wt)
        P (numpy.array): point in wing coordinate system

    Returns:
        p1,p2 (numpy.array): points in propeller coordinate system for the
            left and right propeller
    """

    p1 = P - propeller.loc
    p2 = P - propeller.loc * np.array([1, -1, 1])
    p1 = np.matmul(propeller.prop_rot.T, p1)  # right prop
    p2 = np.matmul(propeller.prop_rot.T, p2)  # left prop
    p2 = p2 * np.array([1, -1, 1])

    return p1, p2

if __name__ == '__main__':
    import test_config as tc

    ib = tc.InboardPropeller
    wt = tc.WTMP
    wing = tc.Wing
    options = tc.Options
    fc = tc.FlightConditions
    propellers = [ib]
    wingsys = WingSys( fc, ib, wt, wing, options, propellers, plot=True)


#plt.plot(P_lst,P1_lst)
#plt.plot(P_lst,P2_lst)


# import pickle
# with open('Vpert', 'wb') as Vpert_file:
#     pickle.dump(vlm.Vpert, Vpert_file)

# with open('Vpert', 'rb') as Vpert_file:
#    vlm_load = pickle.load(Vpert_file)
#     vlm.Vpert['Vx'] = vlm_load['Vx']
#     vlm.Vpert['Vz'] = vlm_load['Vz']
#     vlm.Vpert['Vzi'] = vlm_load['Vzi']

# import os
# airfoil_folder = 'data_aero_wing/laminar/'
# dir_lst = os.listdir(airfoil_folder)
# dir_lst = [i for i in dir_lst]
# d1 = dir_lst[2]
# d2 = dir_lst[3]
# df1 = pd.read_csv(airfoil_folder+'/'+d1, header=6)
# df2 = pd.read_csv(airfoil_folder+'/'+d2, header=6)
# df2 = df2.set_index('alpha')
# df1 = df1.set_index('alpha')
# df2 = df2.astype('float64')
# df1 = df1.astype('float64')