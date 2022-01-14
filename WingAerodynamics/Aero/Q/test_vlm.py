import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import copy
import time

if 'BEM/' not in sys.path:
    sys.path.append('BEM/')

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from vlm import PyVLM  # vlm_or = original non-adapted vlm

def test_vlm(M,rho,a,b,c_r,c_t,le_t,alpha, airfoil,n):

    le_1 = np.array([0, 0])
    le_3 = np.array([le_t, b / 2])

    le = [le_1, le_3]
    chord = [c_r, c_t]

    n = [n]  # chordwise
    # n = [8,8]
    m = [30]  # spanwise
    opt = ['cos']

    vlm = PyVLM()

    airfoil = np.loadtxt('../Data/Airfoil/Wing/' + airfoil + '.dat', skiprows=1)
    vlm.airfoil.input_selig(airfoil[:, 0], airfoil[:, 1], airfoil)

    vlm.add_wing(le, chord, n, m, opt)
    #vlm.polar_dir = 'data_aero_wing/Roger_ex/'
    vlm.alpha = alpha
    vlm.Vinf = M*a
    vlm.rho = rho


    vlm.vlm()
    # vlm.Vpert['turb']=True
    # vlm.vlm_visc()
    vlm.strip_properties(False)
    print(vlm.res)
    return vlm


if __name__ == '__main__':
    M = 0.6
    rho = 0.54
    a = 344
    b = 34
    c_r = 4.2
    c_t =1.5
    le_t = 4.29
    alpha = 2.7
    airfoil = 'ATR72Smoothed' # NACA663-418 ATR72Smoothed
    n = 12
    alpha_lst =[]
    cl_lst = []
    for alpha in [2.7]:
        res = test_vlm(M,rho,a,b,c_r,c_t,le_t,alpha, airfoil, n)
        cl_lst.append(res.res['CL'])
        alpha_lst.append(alpha)


    import matplotlib.pyplot as plt
    plt.plot(alpha_lst,cl_lst)