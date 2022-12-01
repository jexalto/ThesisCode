__author__ = "Bernardo Pacini"
__email__ = "bpacini@umich.edu"
__date__ = "May 14th, 2021 Fri"
__status__ = "Production"

# Local Imports
# from asyncio.windows_events import NULL
from airfoils import airfoilAnalysis

# External Imports
import numpy as np
import json

from proplib import radius # this one is forn xfoil generated data
from airfoil_data_table import airfoilSecs_table as airfoilSecs # this one is for quinty's data
# =============================================================================
# Define Data
# =============================================================================
# Radius

# Provided Radial Sections
nSec = 20
r_sec = np.empty(0)
for iSec in range(0, np.size(airfoilSecs)):
    r_sec = np.append(r_sec, airfoilSecs[iSec]["r"])

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
    _ = airfoilAnalysis(plotting=False, screwXFOIL=True)

    f = open('data/airfoilsecs.json')
    airfoilSecs = json.load(f)
    f.close()
    # =========================================================================
    # Compute Data
    # =========================================================================
    data = {}
    # data["ref_point"] = np.array([0.0, r_sec[0] * radius, 0.0]).tolist()

    # Compute Chord
    chord = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        chord = np.append(chord, airfoilSecs[iSec]["chord"])

    data["chord"] = chord.tolist()

    # # Compute Twist
    theta = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        theta = np.append(theta, airfoilSecs[iSec]["theta"])

    data["theta"] = theta.tolist()

    # # Compute Thick
    thick = np.empty(0)
    thick_ = 0.
    for iSec in range(0, np.size(airfoilSecs)):
        thick = np.append(thick, thick_)

    data["thick"] = chord.tolist()

    # # Compute Alpha_0
    alpha_0 = np.empty(0)
    alpha_0_ = 8.*(2*np.pi)/360
    for iSec in range(0, np.size(airfoilSecs)):
        alpha_0 = np.append(alpha_0, alpha_0_)
    # Re=41000
    alpha0_data = np.array([5., 5., 5., 5., 11, 7., 7., 5., 5., 7., 8., 8., 9., 8.5, 7.5, 7.5, 7.5, 7.5, 7.5, 6., 7.])+2.
    data["alpha_0"] =  ((alpha0_data)/360*(2*np.pi)).tolist()

    # Compute Alpha_L0
    alpha_L0 = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        alpha_L0_section=airfoilSecs[iSec]["alpha_L0"]
        if not alpha_L0_section: # catch empty entries
            alpha_L0_section = -0.06363470452771325
        alpha_L0 = np.append(alpha_L0, alpha_L0_section)
    alpha_L0 = alpha_L0*0.85
    data["alpha_L0"] = alpha_L0.tolist()

    # Compute Cl_alpha
    Cl_alpha = np.empty(0)
    for iSec in range(0, np.size(airfoilSecs)):
        Cl_alpha_section=airfoilSecs[iSec]["Cl_alpha"]
        if not Cl_alpha_section: # catch empty entries
            Cl_alpha_section = 4.
        Cl_alpha = np.append(Cl_alpha, Cl_alpha_section)

    data["Cl_alpha"] = Cl_alpha.tolist()

    # Compute M
    data["M"] = interpVal(r_sec, M_r, M).tolist()

    # Compute Span
    data["span"] = computeSpan(r_sec, radius).tolist()

    # =========================================================================
    # Write Data
    # =========================================================================
    datadir = '/home/jexalto99/code/MDO_lab_env/ThesisCode/HELIX_verification/data/'
    writeJSON(data, datadir+fileName)


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