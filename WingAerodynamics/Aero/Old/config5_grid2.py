""" Based on  TP_tip_med_02.xml run
Copy of config4 but data read from mat file"""
import numpy as np
import sys
import scipy.io
import numpy as np
from scipy import interpolate as si
import pandas as pd
import mat4py

data = mat4py.loadmat('ATR72.mat')

fc = data['data']['fc']
wing = data['data']['wing']
ib = data['data']['ib']
wt = data['data']['wt']


if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

from BEM import BEM

from Q.get_f import Get_f


# --------------- Inboard Propeller ----------------
class InboardPropeller:
    # Vinf = 30
    # J = 1.2
    B = ib['B']  # number of blades
    R = ib['R'] # propeller radius
    D = 2 * R  # propeller diameter
    pitch = ib['pitch']  # pitch angle
    rR_hub = ib['rR_hub']  # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    RPM = ib['RPM']  # rotations per minute
    omega = RPM/60 * (2*np.pi)#Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = 1  # 1 = outboard up, -1 = inboard up

    #rR_beta = 0.75

    rR = np.array(ib['rR']) # blade section locations [r/R]
    #t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    cR = np.array(ib['cR'])  # chord distribution
    theta =  np.array(ib['theta']) # twist distribution
    #theta += pitch  # blade twist

    #theta0 = np.interp(rR_beta, rR, theta)
    #theta = theta - theta0
    x = ib['x']    # spanwise position prop wrt LE [m]
    z = ib['z'] # vertical position prop
    y = ib['y']
    prop_ang = 0
    prop_angle_ = -np.deg2rad(prop_ang)
    prop_rot = np.array([[np.cos(prop_angle_), 0, np.sin(prop_angle_)],
                         [0, 1, 0],
                         [-np.sin(prop_angle_), 0, np.cos(prop_angle_)]])


    bem = BEM()
    bem.R = R
    bem.B = B
    bem.beta = pitch

    bem.theta_b = theta
    bem.cR_b = cR
    bem.rR_b = rR
    bem.rR_hub = rR_hub
    # bem.airfoil_folder = 'data_aero_prop/ATR72/'
    # Airfoil not known per section so polar of complete blade is used
    polar_mat_file = 'PolarsData.mat'
    get_f = Get_f()
    get_f.analyse()
    bem.f_cl = get_f.f_cl
    bem.f_cd = get_f.f_cd

    urot = None
    loc = None

    grid_N = 20  # propeller split in number of panels (N x N)


# --------------- WTM Propeller ----------------
class WTMP:
    # Vinf = 30
    # J = 1.2
    B = wt['B']  # number of blades
    R = wt['R'] # propeller radius
    D = 2 * R  # propeller diameter
    pitch = wt['pitch']  # pitch angle
    rR_hub = wt['rR_hub'] # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    RPM = wt['RPM']
    omega = RPM/60 * 2 *np.pi #Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = 1  # 1 = outboard up, -1 = inboard up

    #rR_beta = 0.75

    rR = np.array(wt['rR']) # blade section locations [r/R]
    # t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    cR = np.array(wt['cR'])  # chord distribution
    theta = np.array(wt['theta'])  # twist distribution
    #theta += pitch  # blade twist

    #theta0 = np.interp(rR_beta, rR, theta )
    #theta = theta - theta0
    x = wt['x']  # spanwise position prop wrt LE [m]
    z = wt['z']  # vertical position prop
    y = wt['y']  # half spanwise position prop y/b/2

    prop_ang = 0

    prop_angle_ = -np.deg2rad(prop_ang)
    prop_rot = np.array([[np.cos(prop_angle_), 0, np.sin(prop_angle_)],
                         [0, 1, 0],
                         [-np.sin(prop_angle_), 0, np.cos(prop_angle_)]])

    bem = BEM()
    bem.R = R
    bem.B = B
    bem.beta = pitch

    bem.theta_b = theta
    bem.cR_b = cR
    bem.rR_b = rR
    bem.rR_hub = rR_hub
    bem.airfoil_folder = 'data_aero_prop_wt/ATR72/'  # storage of airfoil polars

    urot = None
    loc = None
    grid_N = 20  # propeller split in number of panels (N x N)


class Wing:
    c_r = wing['c_r']  # root chord
    c_t = wing['c_t']  # tip chrord
    b = wing['b']    # span
    le_t = wing['le_t'] # value or false, x location of le tip

    alpha = wing['alpha']  #

    airfoil_dir = 'prop_airfoils/'
    polar_dir = 'new_polars' #'data_aero_wing/ATR72/'
    airfoil = 'NACA663-418'

    n = [8, 8, 8, 8, 8]  # chordwise panels per section
    m = [8, 29, 8, 30, 15]  # spanwise panels per section
    opt = ['nsin', 'uni', 'sin', 'nsin', 'uni']  # spacing type per section


class FlightConditions:
    alt = fc['alt']  # altitude
    rho = fc['rho']  # density
    mu = fc['mu']  # viscosity
    #Tinf = 238.62
    a = fc['a']  # speed of sound
    M = fc['M']
    Vinf = fc['Vinf']


class Options:
    settings = {'dJ': 0.8,                      # 0.5
                     'N_prop_map': 3,           # 3
                     'slipstream_length': 8,    # 2 # at least (c4 + x)*2
                     'dlbda': 0.1,              # 0.1
                     'dlb': 0.005,              # 0.005
                     'p_max': 10,               # 10
                     'lbda_max': 5,             # 5
                     'N_time_step': 30,         # 30
                     'N_slipstream_x': 100,      # 12 # 100 !!!!!!!!!!!!!!!!
                     'conway_s_max': 500,       # 500
                     'conway_ds': 0.1}          # 0.1
    N_iter_max = 3
    turb = True
    visc = True



if __name__ == '__main__':
    vlm_dct = {}
    import time
    import pickle
    for i in np.array([20, 16, 14, 12 , 10, 8]):
        fc = FlightConditions()
        ib = InboardPropeller()
        wt = WTMP()
        wing = Wing()
        options = Options()
        ib.grid_N = i
        wt.grid_N = i
        propellers = [ib, wt]  # propellers to be taken into account
        from Aero.VLM.Q_prop2 import WingSys
        wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
        t0 = time.time()
        res = wingsys.analyse()
        vlm = wingsys.vlm
        vlm_res = vlm.res
        vlm.time = (time.time() - t0)/60
        vlm_dct[i] = vlm

        print(vlm.res)
        print(vlm.time)
        with open('vlm_convstudy_grid', 'wb') as conv_file:
            pickle.dump(vlm_dct, conv_file)

        print('Next')
        del fc, ib, wt, wing, wingsys
    print('Ready')

