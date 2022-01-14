""" Load the configuration and do a clean wing run"""
import numpy as np
import sys
import scipy.io
import numpy as np
from scipy import interpolate as si
import pandas as pd
import mat4py
from config5 import config
from Q_prop2 import WingSys


fc, ib, wt, wing, options = config()    # load configuration
print('done')


def run(alpha):
    vlm_dct = {}
    import time

    wing.alpha = alpha
    wing.n = np.array([12])
    wing.m = np.array([30])
    wing.opt = ['cos']
    propellers = []
    #options.turb = True
    options.visc = False

    wingsys = WingSys(fc, ib, wt, wing, options, propellers, plot=False)
    t0 = time.time()

    #wingsys.vlm.rho = 1.225
    #wingsys.vlm.Vinf = 30
    #wingsys.vlm.vlm()
    #print((time.time()-t0))/60
    wingsys.run_vlm()
    print(wingsys.vlm.res)
    print('Ready')
    return wingsys

res = run(6)

