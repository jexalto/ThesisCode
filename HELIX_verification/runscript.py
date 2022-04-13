import helix.model

import helix.parameters.simparam_def as simparam_def

import helix.geometry.geometry_def as geo_def
import helix.geometry.geometry_def_parametric as geo_def_parametric

import helix.references.references_def as ref_def

# External
from mpi4py import MPI
import numpy as np

import matplotlib.pyplot as plt

import json

import niceplots

# from rotor import generateRotor

niceplots.setRCParams()

# generateRotor()

def main():
    """
    Main runner function for repeatedly running HELIX and plotting the results
    against PROWIM experimental data.
    """
    # =========================================================================
    # Import Experimental Data
    # =========================================================================
    with open("/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/data/prowim_data.json", "r") as file:
        dataPROWIM = json.load(file)
        file.close()

    # =========================================================================
    # Parameters
    # =========================================================================
    V_inf = 40.0
    diameter = 0.237
    nJ = 20

    JMin = np.min(dataPROWIM["J"])
    JMax = np.max(dataPROWIM["J"])

    # =========================================================================
    # Set Up Arrays
    # =========================================================================
    J = np.linspace(JMin, JMax, nJ)
    CT = np.zeros(nJ)

    # =========================================================================
    # Run Helix Cases
    # =========================================================================
    for i in range(0, nJ):
        # Compute Rotational Velocity
        omega = V_inf / (J[i] * diameter) * 2.0 * np.pi
        # Call Helix
        CT[i] = runHelix(omega)

    # =========================================================================
    # Plot Result
    # =========================================================================
    _, ax = plt.subplots(figsize=(10, 7))
    print(CT)

    # Plot Helix Data
    ax.plot(J, CT)

    # Plot Experimental Data
    ax.scatter(dataPROWIM["J"], dataPROWIM["CT"])

    ax.set_xlim([0.6, 1.05])
    ax.set_xlabel("Advance Ratio (J)")

    # ax.set_ylim([0.0, 0.16])
    ax.set_ylabel(r"Coefficient of Thrust ($C_T$)")

    # plt.tight_layout()
    niceplots.adjust_spines(ax, outward=True)

    plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/figures/CT_bs.png')


def runHelix(omega):
    """
    This function sets up and runs a HELIX case to compute the quantities of
    interest needed.
    Parameters
    ----------
    omega : float
        Rotational velocity of rotor (in rad/s)
    """
    # ---------------------- Assign Simulation Parameters ------------------- #
    simparam = set_simparam()

    # -------------------------- Set Reference Frames ----------------------- #
    references_def = set_references()

    # ----------------------------- Build Vehicle --------------------------- #
    geometry_def = geometry_definition(omega)

    # ------------------------------- Run HELIX ----------------------------- #
    model = helix.model.Model(thrust=True, torque=True, moment=False, loads=False)
    model.simparam_def = simparam
    model.geometry_def = geometry_def
    model.references_def = references_def
    model.comm = MPI.COMM_WORLD
    model.initialize_run()

    model.run()
    print(np.sum(model.f_geometry.geo_comp[0].a_ref))

    return model.f_geometry.rotor[0].ct[0]


def set_simparam():
    """
    This funcion sets the simulation parameters used in the HELIX case and
    returns a complete simulation parameters object
    Returns
    -------
    simparam : HELIX simulation paramters object (t_simparam_def)
        Object containing the HELIX simulation parameters used for running a
        HELIX case
    """
    simparam = simparam_def.t_simparam_def()
    simparam.basename = "PROWIM"

    simparam.dt = 0.01
    simparam.t_start = 0.0
    simparam.t_end = 0.02

    simparam.nt_rev = 30

    simparam.v_inf = np.array([0.0, 0.0, -40.0])
    simparam.rho_inf = 1.25

    return simparam


def set_references():
    """
    This function sets the reference frame definitions needed in the HELIX case
    and returns a reference frame definition object that is used to run the
    HELIX case.
    Returns
    -------
    references_def : HELIX reference frame definition object (t_references_def)
        Object containing the HELIX reference frame definition parameters used
        for running a HELIX case
    """
    # Initialize Reference Frame Defintions Holder
    references_def = ref_def.t_references_def()

    # Hub Frame
    Hub = ref_def.t_frame_def()
    Hub.Name = "Hub"
    Hub.Parent = "Root"
    Hub.origin = np.array([0.0, 0.0, 0.0])
    Hub.orientation = np.array([[1.0, 0.0, 0.0], [0.0, 1.0*np.cos(np.deg2rad(-0.2)), 1.0*np.sin(np.deg2rad(-0.2))], [0.0, -1.0*np.sin(np.deg2rad(-0.2)), 1.0*np.cos(np.deg2rad(-0.2))]])

    Hub.moving = False

    # Append to References
    references_def.append_frame(Hub)

    return references_def


def geometry_definition(omega):
    """
    This function sets the geometry definitions needed in the HELIX case case
    and returns a geometry definition object that is used to run the HELIX
    case.
    Returns
    -------
    geometry_definition : HELIX geometry definition object (t_geometry_def)
        Object containing the HELIX goemetry definition parameters used for
        running a HELIX case
    """
    # Load rotor data
    fileName = "/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/data/rotor.json"
    with open(fileName, "r") as file:
        airfoilData = json.load(file)
        file.close()

    # Initialize Geometry Component Definitions Holder
    geometry_def = geo_def.t_geometry_def()

    # ---------------------------- Blade Parameters -------------------------- #
    rotor = geo_def_parametric.t_geometry_def_parametric()
    rotor.CompName = "rotor"
    rotor.CompType = "rotor"
    rotor.RefName = "Hub"

    # Reference Parameters
    N_span = len(airfoilData["span"])
    rotor.ref_point = np.array([0.0, 0.0174195, 0.0]) # airfoilData["ref_point"]
    rotor.ref_chord_frac = 0.5

    # Symmetry Parameters
    rotor.symmetry = False
    rotor.symmetry_point = np.array([0.0, 0.0, 0.0])
    rotor.symmetry_normal = np.array([0.0, 1.0, 0.0])

    # Mirror Parameters
    rotor.mirror = False
    rotor.mirror_point = np.array([0.0, 0.0, 0.0])
    rotor.mirror_normal = np.array([0.0, 1.0, 0.0])

    # Initialize Rotor and Allocate Arrays
    rotor.initialize_parametric_geometry_definition(N_span)

    # Set Multiplicity
    rotor.multiple = True
    rotor.multiplicity = {
        "mult_type": "rotor",
        "n_blades": 4,
        "rot_axis": np.array([0.0, 0.0, 1.0]),
        "rot_rate": omega,
        "psi_0": 0.0,
        "hub_offset": 0.0,
        "n_dofs": 0,
    }
    # ------------------------ Blade Section Definition ---------------------- #
    for iSec in range(len(airfoilData["chord"])):
        rotor.sec[iSec].chord = airfoilData["chord"][iSec]
        rotor.sec[iSec].twist = airfoilData["theta"][iSec]
        rotor.sec[iSec].alpha_0 = airfoilData["alpha_0"][iSec]
        rotor.sec[iSec].alpha_L0 = airfoilData["alpha_L0"][iSec]
        rotor.sec[iSec].Cl_alpha = airfoilData["Cl_alpha"][iSec]
        rotor.sec[iSec].M = airfoilData["M"][iSec]

    # ------------------------- Blade Span Definition ----------------------- #
    for iSpan in range(len(airfoilData["span"])):
        rotor.span[iSpan].span = airfoilData["span"][iSpan]
        rotor.span[iSpan].sweep = 0.0
        rotor.span[iSpan].dihed = 0.0
        rotor.span[iSpan].N_elem_span = 2
        rotor.span[iSpan].span_type = 1

    # Append To Vehicle
    geometry_def.append_component(rotor)

    return geometry_def


if __name__ == "__main__":
    main()