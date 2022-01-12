# Helix Libraries
import helix.model

import helix.parameters.simparam_def as simparam_def

import helix.geometry.geometry_def as geo_def
import helix.geometry.geometry_def_parametric as geo_def_parametric

import helix.references.references_def as ref_def

import helix_pyf90.mod_solver_initialize as f90_solver_initialize
import helix_pyf90.mod_aerodynamic_coefficients as f90_aerodynamic_coefficients
# External
import numpy as np
import pandas as pd

###############################################################################
## Main - Run Script
###############################################################################
def run_helix(rpm, v_inf):
    # ---------------------- Assign Simulation Parameters ------------------- #
    simparam = set_simparam(v_inf)

    # -------------------------- Set Reference Frames ----------------------- #
    references_def = set_references()

    # ----------------------------- Build Vehicle --------------------------- #
    geometry_def = geometry_definition(rpm)

    # ------------------------------- Run HELIX ----------------------------- #
    model = helix.model.Model()
    model.simparam_def = simparam
    model.geometry_def = geometry_def
    model.references_def = references_def
    model.initialize_run()
    model.run()
    return model.f_geometry.rotor[0].thrust_mag[0]


###############################################################################
## Set Simulation Parameters
###############################################################################
def set_simparam(v_inf):
    simparam = simparam_def.t_simparam_def()
    simparam.basename = "MR8x45 Rotor"

    simparam.dt = 0.5
    simparam.t_start = 0.0
    simparam.t_end = 1.0

    simparam.nt_rev = 30

    # Flipping freestream to make Z-axis of rotation
    simparam.v_inf = v_inf
    simparam.rho_inf = 1.25

    return simparam


###############################################################################
## Generate Vehicle
###############################################################################
def geometry_definition(rpm):
    # Retrieve geometry
    data = pd.read_csv('/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/Ilinois_data/da4052_9x6.75_geom.txt', dtype=np.float16, sep=' ', skiprows=(1), header=None).values

    rR = data[:, 0]
    cR = data[:, 1]
    beta = data[:, 2]

    radius = 0.2286/2
    rhub = 0.15*radius

    # Initialize Geometry Component Definitions Holder
    geometry_def = geo_def.t_geometry_def()

    # ==========================================================================
    # Define Rotor 1
    # ==========================================================================
    # ---------------------------- Blade Parameters -------------------------- #
    rotor1 = geo_def_parametric.t_geometry_def_parametric()
    rotor1.CompName = "Rotor1"
    rotor1.CompType = "rotor"
    rotor1.RefName = "Hub"

    # Reference Parameters
    N_span = len(rR)+1
    rotor1.ref_point = np.array([0.0, 0.02023364, 0.0])
    rotor1.ref_chord_frac = 0.5

    # Symmetry Parameters
    rotor1.symmetry = False
    rotor1.symmetry_point = np.array([0.0, 0.0, 0.0])
    rotor1.symmetry_normal = np.array([0.0, 1.0, 0.0])

    # Mirror Parameters
    rotor1.mirror = False
    rotor1.mirror_point = np.array([0.0, 0.0, 0.0])
    rotor1.mirror_normal = np.array([0.0, 1.0, 0.0])

    # Initialize Rotor and Allocate Arrays
    rotor1.initialize_parametric_geometry_definition(N_span)

    rotor1.multiple = True
    rotor1.multiplicity = {
        "mult_type": "rotor",
        "n_blades": 2,
        "rot_axis": np.array([0.0, 0.0, 1.0]),
        "rot_rate": rpm / 60.0 * 2.0 * np.pi,
        "psi_0": 0.0,
        "hub_offset": 0.0,
        "n_dofs": 0,
    }

    # ------------------------ Blade Section Definition ---------------------- #
    for i in range(0, N_span-2):
        rotor1.sec[0].chord = (rR[i+1] - rR[i]) * radius
        rotor1.sec[0].twist = (beta[i+1] + beta[i])/2
        rotor1.sec[0].thick = 0.00216916
        rotor1.sec[0].alpha_0 = 0.3595378
        rotor1.sec[0].alpha_L0 = -0.03490658503
        rotor1.sec[0].Cl_alpha = 5.3407
        rotor1.sec[0].M = 50.0

        # Span 1  ------------------
        rotor1.span[0].span = 0.00126238
        rotor1.span[0].sweep = 0.0
        rotor1.span[0].dihed = 0.0
        rotor1.span[0].N_elem_span = 1
        rotor1.span[0].span_type = 1


    # Chord 1 ------------------
    rotor1.sec[0].chord = 0.01761998
    rotor1.sec[0].twist = 40.6848
    rotor1.sec[0].thick = 0.00216916
    rotor1.sec[0].alpha_0 = 0.3595378
    rotor1.sec[0].alpha_L0 = -0.03490658503
    rotor1.sec[0].Cl_alpha = 5.3407
    rotor1.sec[0].M = 50.0

    # Span 1  ------------------
    rotor1.span[0].span = 0.00126238
    rotor1.span[0].sweep = 0.0
    rotor1.span[0].dihed = 0.0
    rotor1.span[0].N_elem_span = 1
    rotor1.span[0].span_type = 1

    # Append To Vehicle
    geometry_def.append_component(rotor1)

    return geometry_def


# ###############################################################################
# ## Generate Vehicle
# ###############################################################################
def set_references():
    # Initialize Reference Frame Defintions Holder
    references_def = ref_def.t_references_def()

    # ==========================================================================
    # Hub Frame
    # ==========================================================================
    Hub = ref_def.t_frame_def()
    Hub.Name = "Hub"
    Hub.Parent = "Root"
    Hub.origin = np.array([0.0, 0.0, 0.0])
    Hub.orientation = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    Hub.moving = False

    # Append to References
    references_def.append_frame(Hub)

    return references_def
