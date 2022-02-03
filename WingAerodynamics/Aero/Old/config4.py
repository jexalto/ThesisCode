""" Based on  TP_tip_med_02.xml run
Copy of config2 but with extra wing section to enable symmetric hs around propeller centre line"""
import numpy as np
import sys

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

from BEM import BEM

from Q.get_f import Get_f


# --------------- Inboard Propeller ----------------
class InboardPropeller:
    # Vinf = 30
    # J = 1.2
    B = 6  # number of blades
    R = 4.24/2 # propeller radius
    D = 2 * R  # propeller diameter
    pitch = 25  # pitch angle
    rR_hub = 0.15  # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    RPM = 964.1952  # rotations per minute
    omega = RPM/60 * (2*np.pi)#Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = 1  # 1 = outboard up, -1 = inboard up

    #rR_beta = 0.75

    rR = np.array([0.150000000000000,0.405000000000000,0.510000000000000,0.615000000000000,0.720000000000000,0.790000000000000,0.860000000000000,0.930000000000000,0.965000000000000,0.999000000000000]) # blade section locations [r/R]
    #t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    cR = np.array([0.165334132226411,0.145600000000000,0.146800000000000,0.151500000000000,0.153400000000000,0.151500000000000,0.138000000000000,0.109400000000000,0.091200000000000,0.068760000000000])  # chord distribution
    theta =  np.array([52.256386365251190,38.870000000000000,33.730000000000000,28.750000000000000,24.270000000000000,22.090000000000000,20.200000000000000,18.470000000000000,17.630000000000000,17.017999999999997]) # twist distribution
    #theta += pitch  # blade twist

    #theta0 = np.interp(rR_beta, rR, theta)
    #theta = theta - theta0
    x = -2.9488    # spanwise position prop wrt LE [m]    #TODO: update
    z = -0.3311  # vertical position prop                 #TODO: update
    y = 4.2719
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
    B = 6  # number of blades
    R = 1.6238/2 # propeller radius
    D = 2 * R  # propeller diameter
    pitch = 25  # pitch angle
    rR_hub = 0.15  # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    RPM = 2519
    omega = RPM/60 * 2 *np.pi #Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = 1  # 1 = outboard up, -1 = inboard up

    #rR_beta = 0.75

    rR = np.array([0.150000000000000,0.405000000000000,0.510000000000000,0.615000000000000,0.720000000000000,0.790000000000000,0.860000000000000,0.930000000000000,0.965000000000000,0.999000000000000]) # blade section locations [r/R]
    # t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    cR = np.array([0.165334132226411,0.145600000000000,0.146800000000000,0.151500000000000,0.153400000000000,0.151500000000000,0.138000000000000,0.109400000000000,0.091200000000000,0.068760000000000])  # chord distribution
    theta = np.array([52.256386365251190,38.870000000000000,33.730000000000000,28.750000000000000,24.270000000000000,22.090000000000000,20.200000000000000,18.470000000000000,17.630000000000000,17.017999999999997])  # twist distribution
    #theta += pitch  # blade twist

    #theta0 = np.interp(rR_beta, rR, theta )
    #theta = theta - theta0
    x = 0.1986  # spanwise position prop wrt LE [m] #TODO: update
    z = 0  # vertical position prop
    y = 17.05  # half spanwise position prop y/b/2

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
    bem.airfoil_folder = 'data_aero_prop_wt/ATR72/'  # storage of airfoil polars #TODO: update polars to 'N250'

    urot = None
    loc = None
    grid_N = 20  # propeller split in number of panels (N x N)


class Wing:
    c_r = 3.615  # root chord
    c_t = 1.624  # tip chrord
    b = 34.1    # span
    le_t = 0.6088  # value or false, x location of le tip

    alpha = 1.2  #

    airfoil_dir = 'prop_airfoils/'
    polar_dir = 'new_polars' #'data_aero_wing/ATR72/'
    airfoil = 'NACA663-418'            # TODO: change to 'N663418'

    n = [8, 8, 8, 8, 8]  # chordwise panels per section
    m = [8, 29, 8, 30, 15]  # spanwise panels per section
    opt = ['nsin', 'uni', 'sin', 'nsin', 'uni']  # spacing type per section


class FlightConditions:
    alt = 7620  # altitude
    rho = 0.5489  # density
    mu = 1.55e-05  # viscosity
    Tinf = 238.62
    a = (1.4 * 287 * Tinf) ** 0.5  # speed of sound
    M = 0.6
    Vinf = M * a


class Options:
    settings = {'dJ': 0.5,                      # 0.5
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
    fc = FlightConditions()
    ib = InboardPropeller()
    wt = WTMP()
    wing = Wing()
    options = Options()
    propellers = [ib, wt]  # propellers to be taken into account
    from Aero.VLM.Q_prop2 import WingSys
    wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=True)
    res = wingsys.analyse()
    print('Done')

