""" Loads configuration and analyses"""
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
from Q_prop2 import WingSys


def run(alpha, fp='cruise'):
    fc, ib, wt, wing, options = config()

    #update settings
#    wing.n = np.array([24, 24, 24, 24, 24])
#    wing.m = np.array([6, 21, 6, 21, 11])
    ib.grid_N = 10
    wt.grid = 10
    options.settings['N_slipstream_x'] = 50


    options.N_iter_max = 1

    # update flight conditions
    wing.alpha = alpha
    if fp == 'pull_up':
        fc.M += 0.04
        fc.alt += 3000*0.3048  # TODO: also update rho, mu

    # propellers to be taken into account
    propellers = [ib, wt]

    # initiate wing system
    wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)

    # run
    t0 = time.time()
    wingsys.analyse()

    print((time.time()-t0/60))
    print('Ready')
    return wingsys

run(2)
