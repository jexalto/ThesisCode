# Perform convergence study on propeller wing system
# by varying number of panels
from runme import load_config
from Aero.VLM.Q_prop2 import WingSys
import numpy as np
import copy
import time
import pickle
import matplotlib.pyplot as plt
import pandas as pd

props = ['ib', 'wt']
case = 'n'
factors = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2] # 2.25, 2.5, 2.75, 3]
alpha = 12  # representable for 2.5G pull up

def convergence(props, case, factors, alpha):
    # load configuration and standard settings
    matfile = 'TP_tip_high_02_kink0_ATRairfoil__0.200_0.8735_0.030_.mat'
    fc_, ib_, wt_, wing_, options_ = load_config(matfile, viscous=False, skinfriction=False, iter_max=4)
    res_dct = {}

    wing_.alpha = alpha
    # set default settings to be studied
    if 'ib' in props and 'wt' in props:
        wing_.m = np.array([12, 33, 12, 22, 16])
        wing_.n = np.array([36, 36, 36, 36, 36])
        wing_.opt = ['nsin', 'uni', 'sin', 'nsin', 'uni']
    elif 'wt' in props:
        wing_.n = np.array([12, 12])
        wing_.m = np.array([30, 10])
        wing_.opt = ['cos', 'uni']
    elif 'ib' in props:
        wing_.n = np.array([32, 32, 32, 32])
        wing_.m = np.array([10, 11, 10, 15])
        wing_.opt = ['nsin', 'uni', 'sin', 'nsin']
    elif props == []:
        wing_.n = np.array([16])
        wing_.m = np.array([30])
        wing_.opt = ['cos']

    m = np.array(wing_.m)
    n = np.array(wing_.n)
    opt = wing_.opt

    print('Starting convergence study on ' + case)
    for i in factors:
        # re-load system
        fc = copy.deepcopy(fc_)
        ib = copy.deepcopy(ib_)
        wt = copy.deepcopy(wt_)
        wing = copy.deepcopy(wing_)
        options = copy.deepcopy(options_)
        if 'ib' in props and 'wt' in props:
            propellers = [ib, wt]
        elif 'ib' in props :
            propellers = [ib]
        elif 'wt' in props:
            propellers = [ wt]
        else:
            propellers = []
            m = np.array([30]) # if no propellers change number of strips and spacing to other default
            wing.opt = ['cos']
        res = {}

        # adapt settings for convergence study
        if case == 'n':
            n_ = (n * i).round()    # new n
            wing.n = n_.astype(int).tolist()
            print(wing.n)

        elif case == 'm':
            m_ = (m * i).round()    # new m
            if propellers != []:
                if m_[1] % 2 == 0:
                    m_[1] += 1          # uneven number of strips in ib propeller region
            wing.m = m_.astype(int).tolist()
            print(wing.m)

        # load wing system and analyse
        wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
        t0 = time.time()
        if propellers ==[]:
            wingsys.run_vlm()
        else:
            wingsys.analyse()
        toc = (time.time() - t0) / 60

        # collect data and save
        vlm = wingsys.vlm
        res['n'] = wing.n
        res['m'] = wing.m
        res['time'] = toc
        res['vlm_res'] = vlm.res
        res['strip'] = vlm.Strip
        res['vpert'] = vlm.Vpert
        if len(wingsys.propellers) > 0:
            res['ct1'] = wingsys.propellers[0].urot.prop_us['integral_CT']
        if len(wingsys.propellers) > 1:
            res['ct2'] = wingsys.propellers[1].urot.prop_us['integral_CT']

        res_dct[i] = res
        file_name = 'Convergence_study/convergence_' + case + 'cruise'
        if propellers == []:
            file_name = file_name + '_clean'

        with open(file_name, 'wb') as conv_file:
            pickle.dump(res_dct, conv_file)

        # delete wingsystem
        del fc, ib, wt, wing, options, propellers, wingsys, vlm
        print(case +' x' + str(i) + ' done')
        print('Next')
    return res_dct

# for r in res_dct.keys():
#     res = res_dct[r]
#     plt.scatter(r,res['vlm_res']['CL'])

#case: clean
# n = 32
# m = 60

# case full
#Run 1 :
    # wing_.m = [8, 21, 8, 15, 11]
    # wing_.n = np.array([24, 24, 24, 24, 24])

# Run 2 :
# wing_.m = [8, 21, 8, 15, 11] * 1.5
# wing_.n = np.array([24, 24, 24, 24, 24]) *1.5 as base


# cruise 1:
# n = 8 >> *3
# m = np.array([9,29,8,30,15])


res = convergence(props, case, factors, alpha)