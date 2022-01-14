import numpy as np
import copy
import sys
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as si

from Q.get_f import Get_f
import mat4py

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

from BEM import BEM


def load_propeller(matfile):
    # dir = '../Data/Design/Propellers/'
    # file = dir + matfile
    # print('loading propeller...')
    # data = mat4py.loadmat(file)  # input file which is used to create object
    print('propeller loaded')

    # load mat data
    # prop = data['data']['propeller']
    # fc = data['data']['fc']

    # load flight conditions
    Vinf = 60           # fc['Vinf']
    rho = 1.225         # fc['rho']
    a = fc['a']
    mu = fc['mu']

    B = prop['B']  # number of blades
    R = prop['R']  # propeller radius
    D = 2 * R  # propeller diameter
    pitch = prop['pitch']  # pitch angle
    rR_hub = prop['rR_hub']  # hub size [r/R]
    RPM = prop['RPM']  # rotations per minute
    omega = 1 * RPM / 60 * (2 * np.pi)  # Vinf / (J * D) * 2 * np.pi  # rotational speed [rad/s]

    # load blade
    blade_name = prop['blade_name']
    rR_beta = 0.7
    rR = np.array(prop['rR'])  # blade section locations [r/R]
    cR = np.array(prop['cR'])  # chord distribution
    theta = np.array(prop['theta'])  # twist distribution
    theta0 = np.interp(rR_beta, rR, theta) # twist is zero @ rR_beta, here pitch is given
    theta = theta - theta0

    # initiate BEM
    bem = BEM()
    bem.Vinf = Vinf
    bem.omega = omega
    bem.rho = rho
    bem.a = a
    bem.mu = mu
    bem.airfoil_folder = '../Data/Polars/Propeller/' + blade_name  # storage of airfoil polars
    bem.R = R
    bem.B = B
    bem.beta = pitch
    bem.theta_b = theta #+pitch
    bem.cR_b = cR
    bem.rR_b = rR
    bem.rR_hub = rR_hub

    return bem


def sweep(bem, type, range, plot=True):
    n_lst = []
    J_lst = []
    T_lst = []
    Ct_lst = []
    P_lst = []
    eta_lst =[]
    Cp_lst = []
    beta_lst =[]
    beta = bem.beta

    cols = ['n', 'J', 'T', 'Ct','P', 'Cp', 'eta']
    res = pd.DataFrame(columns=cols)

    for J in range:

        if type=='Vinf':
            RPM = bem.omega*60/(2 * np.pi)
            bem.Vinf = J*RPM/60*bem.R*2
        elif type =='n' or 'omega':
            bem.omega = bem.Vinf/(bem.R*2*J)*2*np.pi

        bem.analysis()
        if bem.res_sum['converged']:
            n = bem.omega/(2*np.pi)
            T = bem.res_sum['T']
            P = bem.res_sum['P']
            Ct = bem.res_sum['Ct']
            Cp = bem.res_sum['Cp']
            J = bem.res_sum['J']
            eta = bem.res_sum['eta']

            n_lst.append(n)
            T_lst.append(T)
            J_lst.append(J)
            Ct_lst.append(Ct)
            P_lst.append(P)
            eta_lst.append(eta)
            Cp_lst.append(Cp)
            beta_lst.append(beta)

    res['n'] = n_lst
    res['J'] = J_lst
    res['T'] = T_lst
    res['Ct'] = Ct_lst
    res['P'] = P_lst
    res['Cp'] = Cp_lst
    res['eta'] = eta_lst
    res['beta'] = beta_lst

    if plot:
        plt.title(bem.R*2)
        plt.plot(n_lst, T_lst)
        plt.xlabel('n')
        plt.ylabel('T')

        plt.figure()
        plt.title(bem.R * 2)
        plt.plot(J_lst, T_lst)
        plt.xlabel('J')
        plt.ylabel('T')

        plt.figure()
        plt.title(bem.R * 2)
        plt.plot(J_lst, P_lst)
        plt.xlabel('J')
        plt.ylabel('P')

        plt.figure()
        plt.title(bem.R * 2)
        plt.plot(J_lst, Ct_lst)
        plt.xlabel('J')
        plt.ylabel('Ct')

        plt.figure()
        plt.title(bem.R * 2)
        plt.plot(J_lst, eta_lst)
        plt.xlabel('J')
        plt.ylabel('eta')

        plt.figure()
        plt.title(bem.R * 2)
        plt.plot(J_lst, Cp_lst)
        plt.xlabel('J')
        plt.ylabel('Cp')
    return res

def validate():
    # validate the N250 propeller with experimental data
    bem = load_propeller('N250.mat')
    data = {}
    J_range = np.linspace(0.4, 1.8, 15)
    for beta in [25, 27.5, 30, 32.5, 35]:
        bem.beta = beta
        res = sweep(bem, 'Vinf', J_range, plot=False)
        data[beta] = res

    for beta in [25, 27.5, 30, 32.5, 35]:
        plt.figure(1)
        plt.plot(data[beta]['J'], data[beta]['Ct'])
        plt.figure(2)
        plt.plot(data[beta]['J'], data[beta]['Cp'])

    plt.figure(1)
    plt.ylim(0,0.4)
    plt.xlim(0.2, 1.8)
    plt.xlabel('J')
    plt.ylabel('Ct')
    #plt.legend([25, 27.5, 30, 32.5, 35])

    plt.figure(2)
    plt.ylim(0,0.5)
    plt.xlim(0.2, 1.8)
    plt.xlabel('J')
    plt.ylabel('Cp')
    #plt.legend([25, 27.5, 30, 32.5, 35])

    df = pd.concat([data[25],data[27.5], data[30],data[32.5], data[35]])
    return df


def get_rpm(bem, data, thrust, plot=False):
    n_lim = (0.91 * bem.a) / bem.R / (2 * np.pi)
    f = si.interp1d(data['T'], data['n'], kind='cubic')
    f2 = si.interp1d(data['n'], data['P'], kind = 'cubic')
    rps = f(thrust)
    P = f2(rps)

    if plot:
        plt.figure(1)
        i = np.linspace(min(f.x),max(f.x),200)
        r = f(i)
        plt.plot(r,i)
    if rps <=n_lim:
        return rps*60, P.tolist()
    else:
        raise ValueError('RPM exceeds allowed RPM for blade tip mach number ')


def size_prop(file,D, T, beta=55):
    # load propeller and set R, beta
    bem = load_propeller(file)
    bem.R = D/2
    bem.beta = beta

    # sweep propeller rotational speed to get performance map
    J_range = np.linspace(0.4, 5, 15)
    data = sweep(bem, 'n', J_range, plot=False)

    # get propeller rpm to get required thrust
    rpm, P = get_rpm(bem, data, T)

    return rpm, P



if __name__ == '__main__':

    df = validate()
    df.to_excel('Validation/N250.xlsx',index=False)
    #file = 'propeller.mat'
    #bem = load_propeller(file)
    #J_range = np.linspace(0.4, 4, 14)
    #bem.R = 4.3158/2 #inboard prop
    #bem.beta = 55 #bem.beta
    #bem.Vinf = bem.Vinf * np.cos(6/180*2*np.pi)
    #res = sweep(bem,'n', J_range, plot=True)


    # V in propeller plane depends on AoA

    #CD_cruise = 0.03
    #D_cruise = CD_cruise*0.5*0.54*185**2*99
    #phi = 0.2

    #T = D_cruise * (1-phi)/2

    #rpm = get_rpm(bem, res, T)

    #rpm = 2290
    #n = rpm/60
    #r = 2.158
    #v = n *2*np.pi*r

    #J range =

    #propeller.mat holds propeller design used in initiator
    # N250 holds propeller design for validation and ISA sealevel conditions