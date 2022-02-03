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
from Aero.VLM.Q_prop2 import wing_to_prop


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
    bem.airfoil_folder = 'data_aero_prop_wt/ATR72/'  # storage of airfoil polars

    urot = None
    loc = None
    grid_N = 20  # propeller split in number of panels (N x N)


class Wing:
    c_r = wing['c_r']  # root chord
    c_t = wing['c_t']  # tip chrord
    b = wing['b']    # span
    le_t = wing['le_t'] # value or false, x location of le tip

    alpha = wing['alpha'] +1.5 #

    airfoil_dir = 'prop_airfoils/'
    polar_dir = 'data_aero_wing/ATR72/'
    airfoil = 'NACA663-418'

    n = [8, 8, 8, 8, 8]  # chordwise panels per section
    m = [8, 19, 8, 30, 15]  # spanwise panels per section
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
    N_iter_max = 1
    turb = True
    visc = True




fc = FlightConditions()
ib = InboardPropeller()
wt = WTMP()
wing = Wing()
options = Options()
propellers = [ib, wt]  # propellers to be taken into account
from Aero.VLM.Q_prop2 import WingSys
import matplotlib.pyplot as plt
wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
self = wingsys
self.analyse_bem()
self.add_slipstream()
self.analyse_prop()


# ------ propeller induced velocities on wing ---------------------------------
Vpert = np.zeros((self.vlm.Vpert.shape[0], 3))
Vperti = np.zeros((self.vlm.Vpert.shape[0], 3))
P_lst = []
P1_lst = []
P2_lst = []
V1x_lst = []
V2x_lst = []
V1z_lst = []
V2z_lst = []
V1x_lst_wt = []
V2x_lst_wt = []
V1z_lst_wt = []
V2z_lst_wt = []
Vx_lst =[]
Vz_lst =[]
Vx_lst_wt =[]
Vz_lst_wt =[]
for i in range(self.vlm.Vpert.shape[0]):
    P_ = np.array([self.vlm.Vpert.iloc[i]['x'],
                   self.vlm.Vpert.iloc[i]['y'],
                   0])
    P_lst.append(P_[1])
    for prop in self.propellers:
        urot = prop.urot
        P1, P2 = wing_to_prop(prop, P_)
        P1_lst.append(P1[1])
        P2_lst.append((P2[1]))

        Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
        Vt2, Vw2 = urot.slipstream.induced_velocity(P2)

        V1 = np.matmul(prop.prop_rot, Vt1)
        V2 = np.matmul(prop.prop_rot, Vt2)

        V2 = V2 * np.array([1, -1, 1])  # left prop

        V1 = V1 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])
        V2 = V2 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])

        if urot == ib.urot:
            V1x_lst.append(V1[0])
            V2x_lst.append(V2[0])
            V1z_lst.append(V1[2])
            V2z_lst.append(V2[2])
            Vx_lst.append(V1[0]+V2[0])
            Vz_lst.append(V1[2] + V2[2])
        else:
            V1x_lst_wt.append(V1[0])
            V2x_lst_wt.append(V2[0])
            V1z_lst_wt.append(V1[2])
            V2z_lst_wt.append(V2[2])
            Vx_lst_wt.append(V1[0] + V2[0])
            Vz_lst_wt.append(V1[2] + V2[2])

        Vpert[i] = Vpert[i] + V1 + V2  # add induced velocities by both props

        # velocities for trefftz plane
        xi = np.max(urot.slipstream.x) / 2
        x0 = P1[0]
        if urot.slipstream.r_r0 is not None:

            sl_r_r0 = urot.slipstream.r_r0
            sl_x = urot.slipstream.x.flatten()
            sl_z = urot.slipstream.z.flatten()
            i_xi = np.argmin(np.abs(sl_x - xi))
            i_x0 = np.argmin(np.abs(sl_x - x0))

            z0 = np.interp(x0, sl_x, sl_z)
            zi = np.interp(xi, sl_x, sl_z)
            r_r0 = sl_r_r0[i_x0].flatten()

            r0 = (P1[1] ** 2 + P1[2] ** 2) ** 0.5

            if r0 > np.max(r_r0) * prop.R:
                r_r = 1
            else:
                i_r = np.argmin(np.abs(r0 - r_r0 * prop.R))
                r_r = sl_r_r0[i_xi, 0, i_r] / r_r0[i_r]

            P1[1] = P1[1] * r_r
            P1[2] = (P1[2] - z0) * r_r + zi
            P2[1] = P2[1] * r_r
            P2[2] = (P2[2] - z0) * r_r + zi

        P1[0] = xi
        P2[0] = xi

        Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
        Vt2, Vw2 = urot.slipstream.induced_velocity(P2)

        V1 = np.matmul(prop.prop_rot, Vt1)
        V2 = np.matmul(prop.prop_rot, Vt2)

        V2 = V2 * np.array([1, -1, 1])  # left prop
        V1 = V1 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])
        V2 = V2 * np.array([1, -1 * prop.sign_rotation, -1 * prop.sign_rotation])

        Vperti[i] = Vperti[i] + V1 + V2



self.vlm.Vpert['Vx'] = Vpert[:, 0]
self.vlm.Vpert['Vz'] = Vpert[:, 2]
self.vlm.Vpert['Vzi'] = Vperti[:,2]
print('Done')

plt.figure()
plt.plot(P_lst,Vx_lst_wt )
plt.plot(P_lst,Vx_lst )
plt.xlabel('y [m]')
plt.ylabel('Axial induced velocity [m/s]')
plt.legend(['Tip-mounted', 'Inboard'])


plt.figure()
plt.plot(P_lst,Vz_lst_wt )
plt.plot(P_lst,Vz_lst )
plt.xlabel('y [m]')
plt.ylabel('Vertical induced velocity [m/s]')
plt.legend(['Tip-mounted', 'Inboard'])
