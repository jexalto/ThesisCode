__author__ = "Bernardo Pacini"
__email__ = "bpacini@umich.edu"
__date__ = "May 14th, 2021 Fri"
__status__ = "Production"

# Local Imports
from airfoils import airfoilAnalysis

# External Imports
import numpy as np
import json


# =============================================================================
# Define Data
# =============================================================================
# Radius
radius = 1.8288  # [m]

# Provided Radial Sections
nSec = 10
r_sec = np.linspace(0.2, 1.0, nSec)  # []

chord_r = np.array([0.2, 0.92, 1.0])  # []
chord = np.array([0.08, 0.06, 0.05]) * radius  # [m]

twist_r = np.array([0.2, 0.75, 1.0])  # []
twist = np.array([11.0, 0.0, -2.0])  # [deg]

thick_r = np.array([0.2, 0.45, 0.7, 1.0])  # []
thick = np.array([0.12, 0.12, 0.05, 0.03]) * np.interp(thick_r, chord_r, chord)  # [m]

alpha_0_r = np.array([0.2, 0.45, 0.7, 1.0])  # []
alpha_0 = np.deg2rad(np.array([9.45, 9.45, 6.91, 5.88]))  # [deg] Approximated from XFoil Results

alpha_L0_r = np.array([0.2, 0.45, 0.7, 1.0])  # []
# alpha_L0 = np.array([])  # [rad] Computed above with XFoil Data

Cl_alpha_r = np.array([0.2, 0.45, 0.7, 1.0])  # []
# Cl_alpha = np.array([])  # [1/rad] Computed above with XFoil Data

M_r = np.array([0.2, 0.45, 0.7, 1.0])  # []
M = np.array([50.0, 50.0, 50.0, 50.0])  # []


def generateRotor(fileName="rotor.json"):
    """
    This function generates the rotor geometry used on the NASA N+1 Quadcopter.
    Airfoil data is computed using XFoil (and the airfoilAnalysis function)
    and the data is then interpolated to user-specified radial locations on the
    rotor blade. A formatted JSON file is output ("rotor.json") with the data
    required for running HELIX.
    Parameters
    ----------
    fileName : string
        Name of file to be output as JSON geometry file
    """
    # =========================================================================
    # Run Airfoil Analysis to get Data
    # =========================================================================
    airfoilSecs = airfoilAnalysis(plotting=False)

    # =========================================================================
    # Compute Data
    # =========================================================================
    data = {}
    data["ref_point"] = np.array([0.0, r_sec[0] * radius, 0.0]).tolist()

    # Compute Chord
    data["chord"] = interpVal(r_sec, chord_r, chord).tolist()

    # Compute Twist
    data["twist"] = interpVal(r_sec, twist_r, twist).tolist()

    # Compute Thick
    data["thick"] = interpVal(r_sec, thick_r, thick).tolist()

    # Compute Alpha_0
    data["alpha_0"] = interpVal(r_sec, alpha_0_r, alpha_0).tolist()

    # Compute Alpha_L0
    alpha_L0 = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        alpha_L0 = np.append(alpha_L0, airfoilSecs[iSec]["alpha_L0"])

    data["alpha_L0"] = interpVal(r_sec, alpha_L0_r, alpha_L0).tolist()

    # Compute Cl_alpha
    Cl_alpha = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        Cl_alpha = np.append(Cl_alpha, airfoilSecs[iSec]["Cl_alpha"])

    data["Cl_alpha"] = interpVal(r_sec, Cl_alpha_r, Cl_alpha).tolist()

    # Compute M
    data["M"] = interpVal(r_sec, M_r, M).tolist()

    # Compute Span
    data["span"] = computeSpan(r_sec, radius).tolist()

    # =========================================================================
    # Write Data
    # =========================================================================
    writeJSON(data, fileName)


def interpVal(r, val_r, val):
    """
    This function interpolates a sectional quantity (val) at locations along
    the blade (val_r) to the desired locations (r).
    Parameters
    ----------
    r : ndarray (float)
        Desired radial locations
    val_r : ndarray(float)
        Provided radial locations
    val : ndarray (float)
        Quantity value at provided radial locations
    Returns
    -------
    val_out : ndarray (float)
        Quantity value at desired radial locations
    """
    val_out = np.interp(r, val_r, val)
    return val_out


def computeSpan(r, radius):
    """
    This function computes the span of each section of the rotor by finding the
    difference in radial locations between sections.
    Parameters
    ----------
    r : ndarray (float)
        Provided radial locations
    radius : float
        Total radius of rotor blade
    Returns
    -------
    span_out : ndarray (float)
        Span of each blade section
    """
    span_out = np.diff(r_sec * radius)
    return span_out


def writeJSON(data, fileName):
    """
    This function generates and writes a JSON version of the rotor parameters
    data dictionary to a .json file.
    Parameters
    ----------
    data : dictionary
        Dictionary of rotor sectional and spanwise data
    fileName : string
        Name of file to be output as JSON geometry file
    """
    with open(fileName, "w") as file:
        json.dump(data, file, indent=4)


if __name__ == "__main__":
    generateRotor()