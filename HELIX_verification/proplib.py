import pandas as pd
import numpy as np

diameter = 0.237
radius = diameter/2  # [m]
rho_inf = 1.225  # [kg/m3]
mu_inf = 1.48e-5  # [Pa s]
a_inf = 343  # [m/s]
J = 0.8
v_inf = 40.
n = J*v_inf/diameter

path = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/TUD_data/prowim/geometry/PROWIM.txt'

data = pd.read_csv(path, sep=',')
rR = data.iloc[:,1]
theta = data.iloc[:,2]
cR = data.iloc[:,3]

airfoilSecs = [
    {
        "Name": "Sec. 1: st1",
        "airfoil": "st1.dat",
        "r": 0.2 * radius,  # [m]
        "chord": 0.1463,  # [m]
        "theta": 0.,
        "re": None,  # []
        "mach": None,  # []
        "alphaMin": -6,  # [deg]
        "alphaMax": 16,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 2: st2",
        "airfoil": "NACA_5412.txt",
        "r": 0.45 * radius,  # [m]
        "chord": 0.1335,  # [m]
        "theta": 0.,
        "re": None,  # []
        "mach": None,  # []
        "alphaMin": -6,  # [deg]
        "alphaMax": 16,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 3: st3",
        "airfoil": "NACA_5405.txt",
        "r": 0.7 * radius,  # [m]
        "chord": 0.1210,  # [m]
        "theta": 0.,
        "re": 0.621e6,  # []
        "mach": 0.277,  # []
        "alphaMin": -5,  # [deg]
        "alphaMax": 7,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    },
    {
        "Name": "Sec. 4: st4",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
        "theta": 0.,
        "re": 0.671e6,  # []
        "mach": 0.396,  # []
        "alphaMin": -2,  # [deg]
        "alphaMax": 6.0,  # [deg]
        "alpha": None,  # [deg]
        "cl": None,  # []
        "cd": None,  # []
        "cm": None,  # []
        "alpha_L0": None,  # [rad]
        "Cl_alpha": None,  # [1/rad]
    }
]


for index, iAirfoilSec in enumerate(airfoilSecs):
    rps = n
    v_rot = 2*np.pi*rR[index]*radius*rps
    v_tot = np.sqrt(v_inf**2+v_rot**2)
    Mach = v_tot/a_inf
    Re = rho_inf*v_tot*cR[index]*radius/mu_inf

    name = 'Section '+str(index+1)+': st'+str(index+1)

    iAirfoilSec['Name']         = name
    iAirfoilSec['airfoil']      = 'st'+str(index+1)+'.dat'
    iAirfoilSec['r']            = rR[index]
    iAirfoilSec['chord']        = cR[index]*radius
    iAirfoilSec['theta']        = theta[index]
    iAirfoilSec['re']           = Re
    iAirfoilSec['mach']         = Mach
    iAirfoilSec['alphaMin']     = -8.
    iAirfoilSec['alphaMax']     = 4.

None