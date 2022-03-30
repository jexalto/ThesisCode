__author__ = "Bernardo Pacini"
__email__ = "bpacini@umich.edu"
__date__ = "May 14th, 2021 Fri"
__status__ = "Production"

# External Imports
import warnings
import numpy as np
from scipy.interpolate import Akima1DInterpolator
from scipy.optimize import root
import matplotlib.pyplot as plt

import sys
import os.path
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname('/home/jexalto/code/MDO_lab_env/packages/pyXLIGHT/pyxlight/'), os.path.pardir)))

import pyxlight.pyXLIGHT as pyXLIGHT
import niceplots

# Set Plotting Parameters
niceplots.setRCParams()
niceColors = niceplots.get_niceColors()

# Data Directory
dataDir = "./Airfoils/"

# =============================================================================
# General Parameters
# =============================================================================
radius = 1.8288  # [m]
omega = 75  # [rad/s]
rho_inf = 0.98437  # [kg/m3]
mu_inf = 1.84e-5  # [Pa s]
a_inf = 346.204  # [m/s]

# Define Airfoil Sections
airfoilSecs = [
    {
        "Name": "Sec. 1: NACA 5412",
        "airfoil": "NACA_5412.txt",
        "r": 0.2 * radius,  # [m]
        "chord": 0.1463,  # [m]
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
        "Name": "Sec. 2: NACA 5412",
        "airfoil": "NACA_5412.txt",
        "r": 0.45 * radius,  # [m]
        "chord": 0.1335,  # [m]
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
        "Name": "Sec. 3: NACA 5403",
        "airfoil": "NACA_5405.txt",
        "r": 0.7 * radius,  # [m]
        "chord": 0.1210,  # [m]
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
        "Name": "Sec. 4: NACA 1403",
        "airfoil": "NACA_1403.txt",
        "r": 1.0 * radius,  # [m]
        "chord": 0.09144,  # [m]
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
]


def airfoilAnalysis(plotting=True):
    """
    This function handles solving the provided airfoil sections for the rotor
    blade using XFoil, and then computes required parameters such as Cl_alpha
    and alpha_L0.
    NOTE: the methods for computing Cl_alpha and alpha_L0 do not always work,
    so the values should be checked and set manually if required.
    Parameters
    ----------
    plotting : bool, optional
        Flag to enable / disable plotting resulting airfoil data, by default
        True
    Returns
    -------
    airfoilSecs : dictionary
        Dictionary containing airfoil data including Name, alpha, cl, cd, cm
    """
    # =========================================================================
    # Solve Airfoils
    # =========================================================================
    for iSec in range(0, np.size(airfoilSecs)):
        solveAirfoil(airfoilSecs[iSec], omega, rho_inf, mu_inf, a_inf, nPerDeg=10, nIter=1000)

    # =========================================================================
    # Compute Zero-Lift Angle of Attack
    # =========================================================================
    for iSec in range(0, np.size(airfoilSecs)):
        alpha_L0, success = compute_alpha_L0(airfoilSecs[iSec]["alpha"], airfoilSecs[iSec]["cl"])

        if not success:
            warnings.warn(
                "Unable to find zero-lift angle of attack (alpha_L0) for airfoil section {}".format(iSec),
                RuntimeWarning,
            )
        else:
            airfoilSecs[iSec]["alpha_L0"] = alpha_L0

    # =========================================================================
    # Compute Coefficient of Lift vs. Alpha (Cl_alpha)
    # =========================================================================
    alpha_bnds = np.array([0.0, 3.0])
    for iSec in range(0, np.size(airfoilSecs)):
        airfoilSecs[iSec]["Cl_alpha"] = compute_cl_alpha(
            airfoilSecs[iSec]["alpha"], airfoilSecs[iSec]["cl"], alpha_bnds
        )

    # =========================================================================
    # Plot
    # =========================================================================
    if plotting:
        plot(airfoilSecs)

    return airfoilSecs


def solveAirfoil(sec, omega, rho_inf, mu_inf, a_inf, nPerDeg=1, nIter=10):
    """
    This function sets up an XFoil analysis with the airfoil dictionary
    provided and solves for the desired cl, cd, cm values in the desired
    angle of attack range. The dictionary provided is updated and returned.
    Parameters
    ----------
    sec : dictionary
        A dictionary of airfoil section properties: airfoil file, Reynolds
        number, Mach number, minimum angle of attack, maximum angle of attack.
        Entries angle of attack range, cl, cd, and cm are included but will be
        computed and written by this function.
    omega : float
        Angular velocity of rotor, in rad/s
    rho_inf : float
        Free-stream density, in kg/m3
    mu_inf : float
        Free-stream dnamic viscosity, in Pa s
    a_inf : float
        Free-stream sound speed, in m/s
    nPerDeg : integer
        The approximate number of intervals within each degree of angle of
        attack solved by XFoil, by default 10.
    Returns
    -------
    sec: dictionary
        An updated dictionary of airfoil section properties, filling in
        Reynolds number, Mach number, angle of attack range, cl, cd, and cm
        to the dictionary.
    """
    # Compute Reynolds Number and Mach Number
    u = omega * sec["r"]
    re = rho_inf * u * sec["chord"] / mu_inf
    sec["re"] = np.round(re / 1000, decimals=0) * 1000  # Round to nearest thousand

    mach = u / a_inf
    sec["mach"] = np.round(mach, decimals=3)  # Round to nearest thousandth

    # Initialize XFoil Analysis (with dummy analysis)
    airfoilsec = pyXLIGHT.xfoilAnalysis(dataDir + sec["airfoil"], re=sec["re"], mach=sec["mach"], iter=nIter)
    _, _, _, _ = airfoilsec.solveAlpha(0)

    # Setup Angle of Attack Range
    alphaMin = sec["alphaMin"]
    alphaMax = sec["alphaMax"]
    nAlpha = int(np.ceil((alphaMax - alphaMin) * nPerDeg) + 1)
    alphaVec = np.linspace(alphaMin, alphaMax, nAlpha)

    cl = np.empty(0)
    cd = np.empty(0)
    cm = np.empty(0)
    alpha = np.empty(0)

    # Iterate Over Angle of Attack Range
    for iAlpha in range(0, nAlpha):
        cl_loc, cd_loc, cm_loc, lexitflag = airfoilsec.solveAlpha(alphaVec[iAlpha])

        # If Converged, Save Result
        if not lexitflag:
            cl = np.append(cl, cl_loc)
            cd = np.append(cd, cd_loc)
            cm = np.append(cm, cm_loc)

            alpha = np.append(alpha, alphaVec[iAlpha])

    # Store in Dictionary
    sec["cl"] = cl
    sec["cd"] = cd
    sec["cm"] = cm
    sec["alpha"] = alpha

    # Cleanup XFoil Analysis
    del airfoilsec

    return sec


def compute_alpha_L0(alpha, cl):
    """
    This function computes the zero-lift angle of attack of an airfoil
    (alpha_L0) using the angle of attack (alpha) and corresponding coefficient
    of lift (cl).
    Parameters
    ----------
    alpha : ndarray (float)
        One-dimensional array of angle-of-attack data
    cl : ndarray (float)
        One-dimensional array of coefficient of lift data
    Returns
    -------
    alpha_L0 : float
        Zero-lift angle of attack (alpha_L0)
    success : boolean
        Success flag (0 = fail, 1 = success)
    """
    iZero = min(range(len(cl)), key=lambda i: np.abs(cl[i]))

    cl_alpha_interp = Akima1DInterpolator(alpha, cl)
    res = root(cl_alpha_interp, alpha[iZero])

    alpha_L0 = np.deg2rad(res.x)
    success = res.success

    return alpha_L0, success


def compute_cl_alpha(alpha, cl, alpha_bnds):
    """
    This function computes the lift-curve slope (cl_alpha) of an airfoil using
    the angle of attack (alpha) and corresponding coefficient of lift (cl).
    Parameters
    ----------
    alpha : ndarray (float)
        One-dimensional array of angle-of-attack data
    cl : ndarray (float)
        One-dimensional array of coefficient of lift data
    alpha_bnds : ndarray (float)
        One-dimensional, two-element array of angle of attack bounds for
        region used to compute cl_alpha
    Returns
    -------
    cl_alpha : float
        Lift-curve slope (cl_alpha)
    """
    # Find Locations in array closest to the desired angles
    iLow = min(range(len(alpha)), key=lambda i: np.abs(alpha[i] - alpha_bnds[0]))
    iHigh = min(range(len(alpha)), key=lambda i: np.abs(alpha[i] - alpha_bnds[1]))

    # Lower Point
    alpha_low = alpha[iLow]
    cl_low = cl[iLow]

    # Higher Point
    alpha_high = alpha[iHigh]
    cl_high = cl[iHigh]

    # Compute Cl_alpha Value
    cl_alpha = (cl_high - cl_low) / np.deg2rad(alpha_high - alpha_low)

    return cl_alpha


def plot(airfoilSecs):
    """
    This function plots airfoil performance data (cl vs. alpha, cd vs. alpha,
    cm vs. alpha, cd vs. cl).
    Parameters
    ----------
    airfoilSecs : dictionary
        Dictionary containing airfoil data including Name, alpha, cl, cd, cm
    """
    # Coefficient of Lift (Cl) vs. Angle of Attack (alpha)
    fig, ax = plt.subplots(figsize=(10, 8))
    for iSec in range(0, np.size(airfoilSecs)):
        ax.plot(
            airfoilSecs[iSec]["alpha"],
            airfoilSecs[iSec]["cl"],
            label=airfoilSecs[iSec]["Name"],
        )

    # Plot Line at Cl=0
    ax.plot([-10, 20], [0, 0], "--", c=niceColors["Grey"])

    ax.set_xlabel(r"$\alpha$ [$^\circ$]")
    ax.set_ylabel(r"$C_l$")
    ax.set_xlim(-6, 16)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.275))
    niceplots.adjust_spines(ax, outward=True)

    # Coefficient of Drag (Cd) vs. Angle of Attack (alpha)
    fig, ax = plt.subplots(figsize=(10, 8))
    for iSec in range(0, np.size(airfoilSecs)):
        ax.plot(
            airfoilSecs[iSec]["alpha"],
            airfoilSecs[iSec]["cd"],
            label=airfoilSecs[iSec]["Name"],
        )

    ax.set_xlabel(r"$\alpha$ [$^\circ$]")
    ax.set_ylabel(r"$C_d$")
    ax.set_xlim(-4.5, 16)
    ax.legend()
    niceplots.adjust_spines(ax, outward=True)

    # Coefficient of Moment (Cm) vs. Angle of Attack (alpha)
    fig, ax = plt.subplots(figsize=(10, 8))
    for iSec in range(0, np.size(airfoilSecs)):
        ax.plot(
            airfoilSecs[iSec]["alpha"],
            airfoilSecs[iSec]["cm"],
            label=airfoilSecs[iSec]["Name"],
        )

    ax.set_xlabel(r"$\alpha$ [$^\circ$]")
    ax.set_ylabel(r"$C_m$")
    ax.set_xlim(-4.5, 16)
    ax.legend()
    niceplots.adjust_spines(ax, outward=True)

    # Coefficient of Lift (Cl) vs Coefficient of Drag (Cd)
    fig, ax = plt.subplots(figsize=(10, 8))
    for iSec in range(0, np.size(airfoilSecs)):
        ax.plot(
            airfoilSecs[iSec]["cd"],
            airfoilSecs[iSec]["cl"],
            label=airfoilSecs[iSec]["Name"],
        )

    ax.set_xlabel(r"$C_d$")
    ax.set_ylabel(r"$C_l$")
    ax.legend()
    niceplots.adjust_spines(ax, outward=True)

    # Show Plots
    plt.show()


if __name__ == "__main__":
    airfoilAnalysis()