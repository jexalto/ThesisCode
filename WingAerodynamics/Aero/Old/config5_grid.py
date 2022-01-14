""" Loads configuration and explores effect of grid size on results"""
import numpy as np
import sys
import scipy.io
import numpy as np
from scipy import interpolate as si
import pandas as pd
import mat4py
from config5 import config
import time
import pickle


if __name__ == '__main__':
    vlm_dct = {}
    for i in np.array([10, 15, 20, 25, 30, 32,36,38,40,45]):
        fc, ib, wt, wing, options = config()

        wing.n = np.array([24,24,24,24,24])
        wing.m = np.array([6,21,6,21,11])
        ib.grid_N = 10
        wt.grid = 10
        options.settings['N_slipstream_x'] = 50
        m = wing.m
        n = wing.n

        ib.grid_N = i
        wt.grid_N = i
        propellers = [ib, wt]  # propellers to be taken into account
        from Q_prop2 import WingSys
        wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
        t0 = time.time()
        res = wingsys.analyse()
        vlm = wingsys.vlm
        vlm_res = vlm.res
        vlm.time = (time.time() - t0)/60
        vlm_dct[i] = vlm

        print(vlm.res)
        print(vlm.time)
        with open('vlm_convstudy_gridI1', 'wb') as conv_file:
            pickle.dump(vlm_dct, conv_file)

        print('Next')
        del fc, ib, wt, wing, wingsys
    print('Ready')

