""" Reads configuration data from mat file
Make sure the wing and propeller polars are stored in the right place
Settings are base settings"""
import numpy as np
import sys
import scipy.io
import numpy as np
from scipy import interpolate as si
import pandas as pd
import mat4py
from Q_prop2 import WingSys


def load_config(matfile, viscous=False, iter_max=1):

    #matfile = 'ATR72.mat'
    dir = '../Data/Design/'
    file = dir + matfile
    data = mat4py.loadmat(file)  # input file which is used to create object

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
        omega = 1*RPM/60 * (2*np.pi)#Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
        sign_rotation = -1  # 1 = outboard up, -1 = inboard up

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

        # If airfoil not known per section, polar of complete blade is used
        blade_name = ib['blade_name']
        if not ib['airfoils']:
            polar_mat_file = '../Data/Polars/Propeller/' + blade_name + '/PolarsData.mat'
            get_f = Get_f(polar_mat_file)
            get_f.analyse()
            bem.f_cl = get_f.f_cl
            bem.f_cd = get_f.f_cd
        else:
            bem.airfoil_folder = '../Data/Polars/Propeller/' + blade_name  # storage of airfoil polars

        urot = None
        loc = None

        grid_N = 10  # propeller split in number of panels (N x N)


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
        sign_rotation = -1  # 1 = outboard up, -1 = inboard up

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
        blade_name = wt['blade_name']
        if not wt['airfoils']:
            polar_mat_file = '../Data/Polars/Propeller/' + blade_name + 'PolarsData.mat'
            get_f = Get_f()
            get_f.analyse()
            bem.f_cl = get_f.f_cl
            bem.f_cd = get_f.f_cd
        else:
            bem.airfoil_folder = '../Data/Polars/Propeller/' + blade_name  # storage of airfoil polars

        urot = None
        loc = None
        grid_N = 20  # propeller split in number of panels (N x N)


    class Wing:
        c_r = wing['c_r']  # root chord
        c_t = wing['c_t']  # tip chord
        b = wing['b']    # span
        le_t = wing['le_t'] # value or false, x location of le tip

        alpha = wing['alpha']  # + wing['incidence']

        airfoil_dir = '../Data/Airfoil/Wing/'
        airfoil = wing['airfoil']
        polar_dir = '../Data/Polars/Wing/' + airfoil


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
        turb = True
        visc = viscous
        N_iter_max = iter_max

    fc = FlightConditions()
    ib = InboardPropeller()
    wt = WTMP()
    wing = Wing()
    options = Options()

    return fc, ib, wt, wing, options


def main(matfile, alpha, viscous=False, iter_max=1, props=['ib', 'wt']):
    fc, ib, wt, wing, options = load_config(matfile, viscous, iter_max)
    wing.alpha = alpha
    if 'ib' in props and 'wt' in props:
        propellers = [ib, wt]  # propellers to be taken into account
    elif props ==[]:
        wing.n = np.array([12])
        wing.m = np.array([30])
        wing.opt = ['cos']
        propellers = []

    wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
    if propellers==[]:

        wingsys.run_vlm()
    else:
        wingsys.analyse()

    print('Done')
    print(wingsys.vlm.res)
    return wingsys


if __name__ == '__main__':
    res = main('ATR72.mat', 0, props=[])
