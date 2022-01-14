""" Tests Q_prop+ (updated vlm) against TipMountedProp + original vlm
The configuration is based on the tip mounted propeller configuration from 'Wingtip-Mounted Propellers:
Aerodynamic Analysis of Interaction Effects and Comparison with Conventional Layout' """
import numpy as np
import sys
import pandas as pd

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from vlm_or import PyVLM

from BEM import BEM
from bezier import bezier
from Q_prop2 import WingSys


dir_prop_wing = 'Validation/SecIIIC_ModelII/'

fig24 = pd.read_csv(dir_prop_wing+'SecIIIC_Fig24_CLCD_ModelII_tipMounted_ReD640k.txt',
                   header=20)
fig24 = fig24.drop(0)
fig24 = fig24.drop(columns=['config'])
fig24 = fig24.astype('float')

p = 4  # J = 0.8, change when desired


# --------------- Inboard Propeller ---------------- #ignored anyway
class InboardPropeller:
    loc = False
    # # Vinf = 30
    # # J = 1.2
    # B = 4  # number of blades
    # R = False # propeller radius
    # D = 2 * R  # propeller diameter
    # pitch = 25  # pitch angle
    # rR_hub = 0.15  # hub size [r/R]
    # alpha_prop = 0  # angle of attack propeller
    # RPM = 964.1952  # rotations per minute
    # omega = RPM/60 * (2*np.pi)#Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    # sign_rotation = 1  # 1 = outboard up, -1 = inboard up
    #
    # rR_beta = 0.75
    #
    # rR = np.array([0.150000000000000,0.405000000000000,0.510000000000000,0.615000000000000,0.720000000000000,0.790000000000000,0.860000000000000,0.930000000000000,0.965000000000000,0.999000000000000]) # blade section locations [r/R]
    # #t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    # cR = np.array([0.165334132226411,0.145600000000000,0.146800000000000,0.151500000000000,0.153400000000000,0.151500000000000,0.138000000000000,0.109400000000000,0.091200000000000,0.068760000000000])  # chord distribution
    # theta =  np.array([52.256386365251190,38.870000000000000,33.730000000000000,28.750000000000000,24.270000000000000,22.090000000000000,20.200000000000000,18.470000000000000,17.630000000000000,17.017999999999997]) # twist distribution
    # #theta += pitch  # blade twist
    #
    # theta0 = 25  # np.interp(rR_beta, rR, theta)
    # # theta = theta - theta0
    # x = False # spanwise position prop wrt LE [m]
    # z = False  # vertical position prop
    # y_b = False  # spanwise position prop y/b/2
    #
    # prop_ang = 0
    # prop_angle_ = -np.deg2rad(prop_ang)
    # prop_rot = np.array([[np.cos(prop_angle_), 0, np.sin(prop_angle_)],
    #                      [0, 1, 0],
    #                      [-np.sin(prop_angle_), 0, np.cos(prop_angle_)]])
    #
    #
    # bem = BEM()
    # bem.R = R
    # bem.B = B
    # bem.beta = theta0
    #
    # bem.theta_b = theta
    # bem.cR_b = cR
    # bem.rR_b = rR
    # bem.rR_hub = rR_hub
    # # bem.airfoil_folder = 'data_aero_prop/ATR72/'
    # # Airfoil not known per section so polar of complete blade is used
    # polar_mat_file = 'PolarsData.mat'
    # get_f = Get_f()
    # get_f.analyse()
    # bem.f_cl = get_f.f_cl
    # bem.f_cd = get_f.f_cd
    #
    # urot = None
    # loc = False
    #
    # grid_N = 20  # propeller split in number of panels (N x N)


# --------------- WTM Propeller ----------------
class WTMP:
    # Vinf = 30
    # J = 1.2
    B = 4  # number of blades
    R = 0.1184 # propeller radius
    D = 2 * R  # propeller diameter
    pitch = 23.9 - 1.2  # pitch angle
    rR_hub = 0.148  # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    #RPM = 2519
    omega = fig24[fig24['polar'] == p + 1]['n'].mean() * 2 * np.pi #Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = 1  # 1 = outboard up, -1 = inboard up

    rR_beta = 0.75

    prop_geom = pd.read_csv('Validation/SecIIIC_ModelII/prop_geom.txt')
    theta = prop_geom['theta'].to_numpy()
    rR = prop_geom['rR'].to_numpy()
    cR = prop_geom['cR'].to_numpy()

    theta += pitch  # blade twist

    theta0 = np.interp(rR_beta, rR, theta)
    theta = theta - theta0

    x = 0.853 * R  # spanwise position prop wrt LE [m]
    z = 0.042 * R  # vertical position prop
    y_b = 1/2  # spanwise position prop y/b

    prop_ang = 0

    prop_angle_ = -np.deg2rad(prop_ang)
    prop_rot = np.array([[np.cos(prop_angle_), 0, np.sin(prop_angle_)],
                         [0, 1, 0],
                         [-np.sin(prop_angle_), 0, np.cos(prop_angle_)]])

    airfoil_folder = 'data_aero_prop/val/'

    bem = BEM()
    bem.R = R
    bem.B = B
    bem.beta = pitch

    bem.theta_b = theta
    bem.cR_b = cR
    bem.rR_b = rR
    bem.rR_hub = rR_hub
    bem.airfoil_folder = airfoil_folder  # storage of airfoil polars

    urot = None
    loc = None
    grid_N = 20  # propeller split in number of panels (N x N)


class Wing:
    c_r = 0.24  # root chord
    c_t = 0.24 #1.624  # tip chrord
    b = 0.730*2 *0.952 # span

    alpha = -2  # AoA, change when desired

    airfoil_dir = 'prop_airfoils/'
    polar_dir = 'data_aero_wing/val/'
    airfoil = 'NACA64(2)-A015'            #
    n = [24,24]
    #n = [8, 8]  # chordwise panels per section             only sec 1 and 6
    m = [10, 10]  # spanwise panels per section
    opt = ['nsin', 'uni']  # spacing type per section


class FlightConditions:
    # alt = 0  # altitude
    rho = fig24[fig24['polar'] == p + 1]['rhoInf'].mean() # density
    mu = 1.8121e-05  # viscosity
    Tinf = None
    a = (1.4 * 287 * fig24[fig24['polar'] == p + 1]['Tinf'].mean()) ** 0.5  # speed of sound
    Vinf = fig24[fig24['polar']==p+1]['Vinf'].mean()


class Options:
    settings = {'dJ': 0.5,                      # 0.5
                     'N_prop_map': 3,           # 3
                     'slipstream_length': 2,    # 2 # at least (c4 + x)*2
                     'dlbda': 0.1,              # 0.1
                     'dlb': 0.005,              # 0.005
                     'p_max': 10,               # 10
                     'lbda_max': 5,             # 5
                     'N_time_step': 30,         # 30
                     'N_slipstream_x': 12,      # 12 # 100 !!!!!!!!!!!!!!!!
                     'conway_s_max': 500,       # 500
                     'conway_ds': 0.1}          # 0.1
    N_iter_max = 4
    turb = True
    visc=True


if __name__ == '__main__':
    fc = FlightConditions()
    ib = InboardPropeller()
    wt = WTMP()
    wing = Wing()
    options = Options()
    propellers = [wt]  # propellers to be taken into account
    import copy

    # from Q_prop import main, get_instances

    wingsys = WingSys(fc, ib, wt, wing, options, propellers)
    wingsys2 = copy.deepcopy(wingsys)
    wingsys.analyse()

    from TipMountedProp import PropWing

    bem2 = wingsys2.propellers[0].bem

    vlm2 = PyVLM()
    airfoil = np.loadtxt(wing.airfoil_dir + wing.airfoil + '.dat', skiprows=1)
    vlm2.airfoil.input_selig(airfoil[:, 0], airfoil[:, 1], wing.airfoil)
    le = wingsys2.le
    chord = wingsys2.chord
    n = wing.n
    m = wing.m
    opt = wing.opt
    vlm2.Vinf = fc.Vinf
    vlm2.rho = fc.rho
    vlm2.a = fc.a
    vlm2.mu = fc.mu

    vlm2.add_wing(le, chord, n, m, opt)
    vlm2.polar_dir = wing.polar_dir


    propwing = PropWing(vlm2, bem2, wingsys2.propellers[0].loc)
    propwing.N_iter_max = options.N_iter_max
    propwing.settings = options.settings
    alpha = wing.alpha
    propwing.analysis(alpha)
    prop_us = propwing.urot.prop_us
    vlm2 = propwing.vlm

    vlm1 = wingsys.vlm
    print(vlm1.res)
    print(vlm2.res)
    print(wingsys.propellers[0].urot.prop_us['integral_T']*np.sin(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.cos(np.deg2rad(alpha)))
    print(propwing.urot.prop_us['integral_T']*np.sin(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.cos(np.deg2rad(alpha)))
    print(-wingsys.propellers[0].urot.prop_us['integral_T']*np.cos(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.sin(np.deg2rad(alpha)))
    print(-propwing.urot.prop_us['integral_T']*np.cos(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.sin(np.deg2rad(alpha)))
