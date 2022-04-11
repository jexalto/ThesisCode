import numpy as np
import pandas as pd
import warnings
import numpy as np
from scipy.interpolate import Akima1DInterpolator
from scipy.optimize import root
import matplotlib.pyplot as plt

# dir = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/airfoils/N250/airfoils/'
dir = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/xfoil/xfoil_runner/data_output/'

import pandas as pd
import numpy as np

diameter = 0.237
radius = diameter/2  # [m]
rho_inf = 1.225  # [kg/m3]
mu_inf = 1.48e-5  # [Pa s]
a_inf = 343  # [m/s]
J = 0.8
v_inf = 40.
n = (J*v_inf/diameter)**(-1)

path = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/TUD_data/prowim/geometry/PROWIM.txt'

data = pd.read_csv(path, sep=',')
rR = data.iloc[:,1]
theta = data.iloc[:,2]
cR = data.iloc[:,3]

airfoilSecs_table = [
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    },
    {
        "Name": 'st1.dat',
        "span": 0.,
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "cl": 0.,
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
        'alphaMin': None,
        'alphaMax': None,
    }
]

for index, iAirfoilSec in enumerate(airfoilSecs_table):
    name = f'st{index+1}.txt'
    airfoildata = pd.read_csv(dir+name, skiprows=[11], header=5, delim_whitespace=True)
    alpha = airfoildata['alpha'][:]
    alpha = np.where(np.isnan(alpha), 0, alpha)
    cl = airfoildata['CL'][:]
    cl = np.where(cl=='         nan', 0, cl)
    cd = airfoildata['CD'][:]
    cd = np.where(cd=='         nan', 0, cd)
    cm = airfoildata['CM'][:]
    cm = np.where(cm=='         nan', 0, cm)
    rps = n
    v_rot = 2*np.pi*rR[index]*radius*rps
    v_tot = np.sqrt(v_inf**2+v_rot**2)
    Mach = v_tot/a_inf
    Re = rho_inf*v_tot*cR[index]*radius/mu_inf
# names=['alpha', 'CL', 'CD', 'CDp', 'CM', 'Top_Xtr', 'Bot_Xtr']
    name = 'Section '+str(index)+': st'+str(index+1)

    # if index!=0:
    #     iAirfoilSec['span']     = abs(rR[index*2]-rR[index*2-1])*radius
    iAirfoilSec['r']            = rR[index]
    iAirfoilSec['chord']        = cR[index]*radius
    iAirfoilSec['theta']        = theta[index]
    iAirfoilSec['alpha']        = alpha.tolist()
    iAirfoilSec['alphaMin']     = min(alpha.tolist())
    iAirfoilSec['alphaMax']     = max(alpha.tolist())
    iAirfoilSec['cl']           = cl.tolist()
    iAirfoilSec['cd']           = cd.tolist()
    iAirfoilSec['cm']           = cm.tolist()
    iAirfoilSec['name']         = name
None