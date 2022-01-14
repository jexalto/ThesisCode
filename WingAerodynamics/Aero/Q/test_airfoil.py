import sys
if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from airfoils import Airfoil
import numpy as np
import matplotlib.pyplot as plt


def thin_airfoil_polar(airfoil,M):
    af = Airfoil()
    airfoil = np.loadtxt('../Data/Airfoil/Wing/' + airfoil + '.dat', skiprows=1)
    af.input_selig(airfoil[:, 0], airfoil[:, 1], airfoil)

    a = np.linspace(-10,10,20)
    cl = af.cl_thin_airfoil(a, M)
    #plt.plot(a,cl )
    return a, cl

if __name__ == '__main__':
    airfoil = 'NACA2412'
    M=0
    a, cl = thin_airfoil_polar(airfoil, M)
    plt.plot(a, cl)