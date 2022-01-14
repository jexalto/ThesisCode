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


def main(fc, ib, wt, wing, options, propellers):
    t0 = time.time()  # get time
    ib, wt, vlm, settings = get_instances(fc, ib, wt, wing, options, propellers)    # initialise
    propellers, conv_bem = analyse_bem(propellers)                      # analyse BEM(s)
    if conv_bem is True:
        propellers = add_slipstream(propellers, settings)               # add propeller slipstream(s)

        # data frame to store convergence data
        conv = pd.DataFrame(columns=['dCT', 'dCQ', 'dCL', 'dCD'])
        ct0, cq0, cl0, cd0 = 10, 10, 10, 10  # guess
        tmp_dct = {}
        conv_iter = False

        vlm_clean = copy.deepcopy(vlm)  # copy of clean vlm = clean wing
        for N in range(options.N_iter_max):
            propellers = analyse_prop(propellers)
            vlm = prop_induced_velocities(propellers, vlm, ib)      # introduce propeller induced velocities to vlm
            vlm = run_vlm(vlm, options.turb)
            propellers = wing_induced_velocities(vlm, propellers)
            propellers = slipstream_deflection(vlm, propellers)


            # convergence data

            ct1 = propellers[0].urot.prop_us['integral_CT']
            cq1 = propellers[0].urot.prop_us['integral_CQ']
            cl1 = vlm.res['CL']
            cd1 = vlm.res['CD']
            conv.loc[N] = [np.abs(ct0 - ct1), np.abs(cq0 - cq1), np.abs(cl0 - cl1), np.abs(cd0 - cd1)]

            tmp_dct['urot_1_%d' % N] = copy.deepcopy(propellers[0].urot)  # dictionary of ib urot
            if len(propellers) >1 :
                tmp_dct['urot_2_%d' % N] = copy.deepcopy(propellers[1].urot)  # dictionary of wt urot
            tmp_dct['vlm_%d' % N] = copy.deepcopy(vlm)

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
                print('next iteration')

        if not conv_iter:
            N = conv['dCD'].idxmin()

        conv = conv.iloc[:N + 1]
        print('Done')

        print((time.time() - t0) / 60)
        return vlm_clean, tmp_dct, conv
    else:
        return None


def get_instances(fc, ib, wt, wing, options, propellers):
    # ----------------- Flight conditions ---------------------
    Vinf = fc.Vinf  # freestream velocity
    a = fc.a  # speed of sound
    alt = fc.alt  # altitude
    rho = fc.rho  # density
    mu = fc.mu  # viscosity
    Tinf = fc.Tinf

    # ------------- Wing --------------------
    c_r = wing.c_r  # root chord
    c_t = wing.c_t  # tip chrord
    b = wing.b  # span
    wing_ang = wing.alpha

    # ------------ Propellers-----------------
    for prop in propellers:
        prop.loc = np.array([-prop.x, prop.y_b * b, prop.z])  # inboard prop location
        prop.bem.Vinf = Vinf
        prop.bem.omega = prop.omega
        prop.bem.rho = rho
        prop.bem.a = a
        prop.bem.mu = mu
        prop.alpha_prop = wing_ang + prop.prop_ang
        prop.bem.airfoil_folder = prop.bem.airfoil_folder
        prop.bem.Vinf = Vinf

    # ------------------- VLM -------------------------
    n = wing.n  # chordwise panels per section
    m = wing.m  # spanwise panels per section
    opt = wing.opt  # spacing type per section

    le_1 = le_2 = le_3 = le_4 = le_5 = np.array([0, 0])
    if ib in propellers:
        le_2 = np.array([0, ib.loc[1] - ib.R])
        le_3 = np.array([0, ib.loc[1] + ib.R])
        le_4 = np.array([0, le_3[1] + le_2[1]])
    if wt in propellers:
        le_5 = np.array([0, wt.loc[1] - wt.R])
    le_6 = np.array([0, b / 2])

    c_2 = c_r - (le_2[1] - le_1[1]) * (c_r - c_t) / b * 2
    c_3 = c_r - (le_3[1] - le_1[1]) * (c_r - c_t) / b * 2
    c_4 = c_r - (le_4[1] - le_1[1]) * (c_r - c_t) / b * 2
    c_5 = c_r - (le_5[1] - le_1[1]) * (c_r - c_t) / b * 2

    le = [le_1, le_2, le_3, le_4, le_5, le_6]  # wing leading edge coordinates
    le_idx = [0]+ np.nonzero(le)[0].tolist()
    chord = [c_r, c_r, c_r, c_r, c_r, c_r]  # corresponding chord               #TODO: change to local chord

    le = [le[i] for i in le_idx]
    chord = [chord[i] for i in le_idx]

    if ib in propellers:
        if le_2[1] != le_4[1] - le_3[1] or wing.m[0] != wing.m[2] or wing.opt[0] != 'n' + wing.opt[2]:
            raise Exception('Please check that section 1 has the same characteristics as section 3 and an opossing spacing')
            # Check to see if panel symmetry around inboard propeller

    # if m[1]%2==0:
    #    raise Exception('Use uneven number of panels for inboard propeller section')

    vlm = PyVLM()  # initiate VLM
    airfoil = np.loadtxt(wing.airfoil_dir + wing.airfoil + '.dat', skiprows=1)
    vlm.airfoil.input_selig(airfoil[:, 0], airfoil[:, 1], wing.airfoil)  # airfoil name!
    # vlm.NACA = wing.airfoil #'2412'  # set airfoil #TODO
    vlm.add_wing(le, chord, n, m, opt)
    vlm.polar_dir = wing.polar_dir

    vlm.Vinf = Vinf
    vlm.rho = rho
    vlm.a = a
    vlm.mu = mu
    vlm.wt_R = wt.R
    vlm.ib_loc = ib.loc
    vlm.alpha = wing_ang

    # -----------------  Settings --------------------------------------------------
    settings = options.settings
    N_iter_max = options.N_iter_max

    # pass settings
    for key in vlm.jet_corr_settings.keys():
        vlm.jet_corr_settings[key] = settings[key]

    return ib, wt, vlm, settings


# --------------- Initiate --------------------------------------------------
def analyse_bem(propellers):
    conv_bem = True
    for prop in propellers:
        prop.bem.analysis()
        print('bem analysed')
        if not  prop.bem.res_sum['converged']:
            conv_bem = False
    return propellers, conv_bem


def add_slipstream(propellers, settings):
    # create a propeller map
    for prop in propellers:
        if prop.urot is None:
            dJ = settings['dJ']
            N_prop_map = settings['N_prop_map']
            prop.urot = create_unsteady_rotor_model(prop.bem, dJ, N_prop_map)
            print('urot created')
            prop.urot.N_time_step = settings['N_time_step']
            prop.urot.sign_rotation = prop.sign_rotation

        else:
            prop.urot = prop.urot

        # slipstream geometry
        sl = settings['slipstream_length']  # slipstream length
        dsl = prop.bem.Vinf / (prop.bem.omega / (2 * np.pi)) / settings['N_slipstream_x']
        # prevent slipstream with only one x location
        if dsl >= sl:
            dsl = sl / 2.001
        prop.x_st = np.arange(0, sl, dsl)      #TODO
        # prop.x_st = (1 - np.cos(np.linspace(0, 0.5 * np.pi, 100))) * 2
        # prop.z_st = np.zeros(prop.x_st.shape)
        prop.z_st = np.zeros(prop.x_st.shape)

        # slipstream settings
        prop.urot.slipstream.cnw.ds = settings['conway_ds']
        prop.urot.slipstream.cnw.s_max = settings['conway_s_max']

        # set propeller inflow field
        prop.u_aoa = np.cos(np.deg2rad(prop.alpha_prop))
        prop.w_aoa = np.sin(np.deg2rad(prop.alpha_prop))
        prop.grid_y = np.linspace(-prop.R, prop.R, prop.grid_N)
        prop.grid_z = np.linspace(-prop.R, prop.R, prop.grid_N)
        prop.u_in = np.ones((prop.grid_y.size, prop.grid_z.size)) * prop.u_aoa
        prop.v_in = np.zeros((prop.grid_y.size, prop.grid_z.size))
        prop.w_in = np.ones((prop.grid_y.size, prop.grid_z.size)) * prop.w_aoa
    return propellers


def analyse_prop(propellers):
# analyse prop
    for prop in propellers:  # initialize slipstream
        print('analysing prop started')
        prop.urot.inflow.grid_input(prop.grid_y / prop.R, prop.grid_z / prop.R, prop.u_in, prop.v_in, prop.w_in)
        prop.urot.analysis()

        prop.urot.calculate_circulation(prop.bem.res['Va_i'].to_numpy(),
                                  prop.bem.res['Vt_i'].to_numpy(),
                                  prop.bem.res['rR'].to_numpy())

        prop.urot.init_slipstream(prop.x_st, prop.z_st)

        print('analysing prop done')
    return propellers


def prop_induced_velocities(propellers, vlm, ib):
    # ------ propeller induced velocities on wing ---------------------------------
    Vpert = np.zeros((vlm.Vpert.shape[0], 3))
    Vperti = np.zeros((vlm.Vpert.shape[0], 3))
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

    for i in range(vlm.Vpert.shape[0]):
        P_ = np.array([vlm.Vpert.iloc[i]['x'],
                       vlm.Vpert.iloc[i]['y'],
                       0])
        P_lst.append(P_[1])
        for prop in propellers:
            urot = prop.urot
            P1, P2 = wing_to_prop(prop, P_)
            P1_lst.append(P1[1])
            P2_lst.append((P2[1]))

            Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
            Vt2, Vw2 = urot.slipstream.induced_velocity(P2)

            V1 = np.matmul(prop.prop_rot, Vt1)
            V2 = np.matmul(prop.prop_rot, Vt2)

            V2 = V2 * np.array([1, -1, 1])  # left prop

            if urot == ib.urot:
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

            Vperti[i] =  Vperti[i] + V1+V2

    vlm.Vpert['Vx']   = Vpert[:, 0]
    vlm.Vpert['Vz']   = Vpert[:, 2]
    vlm.Vpert['Vzi']   = Vperti[:, 2]

    return vlm


def run_vlm(vlm, turb):
    vlm.G = 0
    if turb:
        vlm.Vpert['turb'] = True       #TODO: turbulent or not?
    vlm.vlm_visc()
    # vlm.vlm()
    # vlm.strip_properties(False)
    return vlm


def wing_induced_velocities(vlm, propellers):
    for prop in propellers:
        for iy, y_ in enumerate(prop.grid_y):
            for iz, z_ in enumerate(prop.grid_z):
                P_ = np.array([0, y_, z_])
                P = prop_to_wing(prop, P_)

                V = vlm.induced_velocity(P)                 # this takes time ~0.15 sec
                V = np.matmul(prop.prop_rot.T, V)

                prop.u_in[iy, iz] = V[0]/vlm.Vinf+prop.u_aoa
                prop.v_in[iy, iz] = V[1]/vlm.Vinf
                prop.w_in[iy, iz] = V[2]/vlm.Vinf+prop.w_aoa
    return propellers


def slipstream_deflection(vlm, propellers):
    # slipstream deflection
    for prop in propellers:
        for i in range(prop.z_st.size-1):
            P_ = np.array([prop.x_st[i], 0, 0])
            P = prop_to_wing(prop, P_)

            w = vlm.induced_velocity(P)
            w = np.matmul(prop.prop_rot.T, w)

            dx = prop.x_st[i+1]-prop.x_st[i]
            prop.z_st[i+1] = prop.z_st[i]+dx*w[2]/vlm.Vinf
    return propellers


    # self.urot = urot
    # self.vlm = vlm
    # self.conv = conv


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

#plt.plot(P_lst,P1_lst)
#plt.plot(P_lst,P2_lst)
plot = False
if plot:
    plt.figure()
    plt.scatter(ib.loc[1], 0, marker='x', color='red')
    plt.plot(P_lst, V1x_lst)
    plt.xlabel('y')
    plt.ylabel('Induced velocity x [m/s]')
    plt.grid()

    plt.figure()
    plt.scatter(ib.loc[1], 0, marker='x', color='red')
    plt.plot(P_lst, V1z_lst)
    plt.xlabel('y')
    plt.ylabel('Induced velocity z [m/s]')
    plt.grid()


    # Plot wing plan form
    plt.figure()
    y_wing = np.array([le_1[1], le_2[1], le_3[1], le_4[1], le_5[1], le_5[1], le_4[1], le_3[1], le_2[1], le_1[1], le_1[1]])
    x_wing = np.array([le_1[0], le_2[0], le_3[0], le_4[0], le_5[0], le_5[0] + c_t, le_4[0] + c_4, le_3[0] + c_3, le_2[0] + c_2, le_1[0] + c_r, le_1[0]])
    plt.plot(y_wing, x_wing)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('y [m]')
    plt.ylabel('x [m]')

    # Plot front view, normalised
    plt.figure()
    for prop in propellers:
        span = wing.b / wing.b
        radius = prop.bem.R / wing.b
        prop_loc_y = prop.loc[1] / wing.b
        prop_loc_z = prop.loc[2] / wing.b

        plt.plot([0, span / 2], [0, 0], color='black')
        pr = plt.Circle((prop_loc_y, prop_loc_z), radius, color='r', fill=False)
        plt.gcf().gca().add_artist(pr)

    plt.xlim(0, 1.1)
    plt.ylim(-wt.R * 1.1, wt.R * 1.1)
    plt.xlabel('y/b [-]')
    plt.ylabel('z/b [-]')





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