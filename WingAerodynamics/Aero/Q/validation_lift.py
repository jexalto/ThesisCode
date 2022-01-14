import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import copy
import time
import pickle

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from vlm import PyVLM
from BEM import BEM
from bezier import bezier
from Q_prop2 import WingSys

bem_dct = {}
vlm_dct = {}
vlm_dct2 = {}
conv_dct = {}

bem_dct_ = {}
vlm_dct_ = {}
vlm_dct2_ = {}
conv_dct_ ={}


### ISOLATED PROP
dir_prop_wing = 'Validation/SecIIIC_ModelII/'

fig24 = pd.read_csv(dir_prop_wing + 'SecIIIC_Fig24_CLCD_ModelII_conventional_ReD640k.txt',
                    header=20)
fig24 = fig24.drop(0)
fig24 = fig24.drop(columns=['config'])
fig24 = fig24.astype('float')

fig24_val = pd.DataFrame(columns=['Polar', 'AoA', 'CL', 'CD', 'conv'], index=np.arange(0, 49, 1))
aoa = [4] #np.linspace(-8, 10, 10)

import numpy as np
import sys

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

from BEM import BEM

#from Q.get_f import Get_f


# --------------- Inboard Propeller ----------------
class InboardPropeller(object):
    # Vinf = 30
    # J = 1.2
    B = 4  # number of blades
    R = 0.1184  # propeller radius
    D = 2 * R  # propeller diameter
    pitch = 23.9 - 1.2  # pitch angle
    rR_hub = 0.148  # hub size [r/R]
    alpha_prop = 0  # angle of attack propeller
    # RPM = required rotations per minute
    omega = None  # required RPM/60 * (2*np.pi)#Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]
    sign_rotation = -1  # 1 = outboard up, -1 = inboard up

    rR_beta = 0.75

    prop_geom = pd.read_csv('Validation/SecIIIC_ModelII/prop_geom.txt')
    theta = prop_geom['theta'].to_numpy()
    rR = prop_geom['rR'].to_numpy()
    cR = prop_geom['cR'].to_numpy()

    # rR = np.array([0.150000000000000,0.405000000000000,0.510000000000000,0.615000000000000,0.720000000000000,0.790000000000000,0.860000000000000,0.930000000000000,0.965000000000000,0.999000000000000]) # blade section locations [r/R]
    # #t = 1 / (max(rR) - min(rR)) * (rR - min(rR))  # distance to center from blade section [r/R]
    # cR = np.array([0.165334132226411,0.145600000000000,0.146800000000000,0.151500000000000,0.153400000000000,0.151500000000000,0.138000000000000,0.109400000000000,0.091200000000000,0.068760000000000])  # chord distribution
    # theta =  np.array([52.256386365251190,38.870000000000000,33.730000000000000,28.750000000000000,24.270000000000000,22.090000000000000,20.200000000000000,18.470000000000000,17.630000000000000,17.017999999999997]) # twist distribution
    # #theta += pitch  # blade twist
    theta += pitch  # blade twist

    theta0 = np.interp(rR_beta, rR, theta)
    theta = theta - theta0

    x = -0.853 * R  # spanwise position prop wrt LE [m]
    z = 0.042 * R  # vertical position prop
    y = 0.444 * 0.748  # spanwise position prop y/b/2

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
    bem.airfoil_folder = airfoil_folder  # storage folder of polars

    # # Airfoil not known per section so polar of complete blade is used
    # polar_mat_file = 'PolarsData.mat'
    # get_f = Get_f()
    # get_f.analyse()
    # bem.f_cl = get_f.f_cl
    # bem.f_cd = get_f.f_cd

    urot = None
    loc = None

    grid_N = 20  # propeller split in number of panels (N x N)


# --------------- WTM Propeller ----------------
class WTMP:  # disregarded
    R = False  # propeller radius
    x = False
    y = False
    z = False


class Wing:
    c_r = 0.24  # root chord
    c_t = 0.24  # tip chrord
    b = 0.730 * 2   # span
    le_t = False

    alpha = None  # variable

    airfoil_dir = 'prop_airfoils/'
    polar_dir = 'data_aero_wing/val/'
    airfoil = 'NACA64(2)-A015'

    n = [24, 24, 24, 24]  # chordwise panels per section
    m = [10, 11, 10, 5]  # spanwise panels per section
    opt = ['nsin', 'uni', 'sin', 'nsin']  # spacing type per section


class FlightConditions:
    # alt = 7620  # altitude
    rho = None  # density
    mu = 1.8121e-05  # viscosity
    Tinf = None
    a = None  # (1.4 * 287 * Tinf) ** 0.5  # speed of sound
    Vinf = None


class Options:
    settings = {'dJ': 0.65,  # 0.5
                'N_prop_map': 3,
                'slipstream_length': 2,
                'dlbda': 0.1 / 0.8,
                'dlb': 0.005,
                'p_max': int(10 * 0.8),
                'lbda_max': 5 * 0.8,
                'N_time_step': int(30 * 0.8),
                'N_slipstream_x': int(12 * 0.8),
                'conway_s_max': 500 * 0.8,
                'conway_ds': 0.1 / 0.8,
                'conway_N_poly': 8}
    N_iter_max = 6
    turb = True
    visc = True


ind = 0

fc = FlightConditions()
ib = InboardPropeller()
wt = WTMP()
wing = Wing()
options = Options()

# for p in [0, 2, 4]:
# aoa = np.linspace(-12, 25, 20)
# aoa = np.linspace(-10, 34, 20)
# aoa = [5]

propellers = [ib]  # propellers to be taken into account

ps = [0,1,2,3,4]
#aoa = np.linspace(-8,10,10)
aoa = np.linspace(10,-8,10)

def val_ib(ps, aoa, name =None):
    ind = 0
    for p in ps:
        # for p in [0, 2, 4]:
        if name is None:
            name = str(p+1)
        if p == 0:
            for alpha in aoa:
                fc = FlightConditions()
                ib = InboardPropeller()
                wt = WTMP()
                wing = Wing()
                options = Options()
                fc.Vinf = fig24[fig24['polar'] == p + 1]['Vinf'].mean()
                fc.rho = fig24[fig24['polar'] == p + 1]['rhoInf'].mean()
                fc.a = (1.4 * 287 * fig24[fig24['polar'] == p + 1]['Tinf'].mean()) ** 0.5
                wing.alpha = alpha
                propellers = []
                wingsys = WingSys(fc, ib, wt, wing, options, propellers)  # initialise

                wingsys.run_vlm()
                vlm = wingsys.vlm

                fig24_val.loc[ind]['CL'] = vlm.res['CL']
                fig24_val.loc[ind]['CD'] = vlm.res['CD']
                fig24_val.loc[ind]['Polar'] = p + 1
                fig24_val.loc[ind]['AoA'] = alpha

                vlm_dct2[alpha] = copy.deepcopy(vlm)
                with open('fig24_data_'+ name, 'wb') as fig_data_file:
                    pickle.dump(fig24_val, fig_data_file)
                ind += 1
                del wingsys, ib, wt, wing, options, propellers  # delete
        else:

            for alpha in aoa:
                fc = FlightConditions()
                ib = InboardPropeller()
                wt = WTMP()
                wing = Wing()
                options = Options()
                propellers = [ib]

                fc.Vinf = fig24[fig24['polar'] == p + 1]['Vinf'].mean()
                fc.rho = fig24[fig24['polar'] == p + 1]['rhoInf'].mean()
                fc.a = (1.4 * 287 * fig24[fig24['polar'] == p + 1]['Tinf'].mean()) ** 0.5

                ib.Vinf = fig24[fig24['polar'] == p + 1]['Vinf'].mean()
                ib.omega = fig24[fig24['polar'] == p + 1]['n'].mean() * 2 * np.pi

                wing.alpha = alpha
                t0 = time.time()
                wingsys = WingSys(fc, ib, wt, wing, options, propellers)  # initialise
                wingsys.analyse()

                print(t0 - time.time())

                ib_prop = wingsys.propellers[0].urot.prop_us
                bem_ib = wingsys.propellers[0].bem
                vlm = wingsys.vlm

                L_t = ib_prop['integral_T'] * np.sin(np.deg2rad(alpha)) + ib_prop['integral_delta_Fz'] * np.cos(
                    np.deg2rad(alpha))
                D_t = -ib_prop['integral_T'] * np.cos(np.deg2rad(alpha)) + ib_prop['integral_delta_Fz'] * np.sin(
                    np.deg2rad(alpha))

                qS = 0.5 * vlm.res['rho'] * vlm.res['Vinf'] ** 2 * vlm.res['S']

                fig24_val.loc[ind]['CL'] = vlm.res['CL'] + 2 * L_t / qS
                fig24_val.loc[ind]['CD'] = vlm.res['CD'] + 2 * D_t / qS
                fig24_val.loc[ind]['Polar'] = p + 1
                fig24_val.loc[ind]['AoA'] = alpha
                fig24_val.loc[ind]['conv'] = wingsys.converged

                vlm_dct[alpha] = copy.deepcopy(vlm)
                bem_dct[alpha] = copy.deepcopy(bem_ib)
                conv_dct[alpha] = copy.deepcopy(wingsys.conv)

                with open('fig24_data_'+ name, 'wb') as fig_data_file:
                    pickle.dump(fig24_val, fig_data_file)

                ind += 1
                del wingsys, ib, wt, wing, options, propellers  # delete

            vlm_dct_[p] = vlm_dct
            bem_dct_[p] = bem_dct
            conv_dct_[p] = conv_dct

            with open('vlm_dct_'+ name, 'wb') as vlm_file:
                pickle.dump(vlm_dct_, vlm_file)

            with open('bem_dct_'+name, 'wb') as bem_file:
                pickle.dump(bem_dct_, bem_file)

            with open('conv_dct_'+name, 'wb') as conv_file:
                pickle.dump(conv_dct_, conv_file)

#val_ib([3], [4])
