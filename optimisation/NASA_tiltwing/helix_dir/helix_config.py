import helix.parameters.simparam_def as py_simparam_def
import helix.references.references_def as py_ref_def
import helix.geometry.geometry_def as geo_def
import helix.geometry.geometry_def_parametric as geo_def_parametric
import json

# Helix OpenMDAO Wrapper
import helix.openmdao.om_helix as om_helix

# External
import openmdao.api as om
import numpy as np

# rst simparam
def simparam_definition():
    simparam = py_simparam_def.t_simparam_def()
    simparam.basename = "exampleOpt"

    simparam.nt = 5
    simparam.t_start = 0.0
    simparam.t_end = 0.1

    simparam.nt_rev = 30

    simparam.v_inf = np.array([0.0, 0.0, -40.0])
    simparam.rho_inf = 1.2087

    return simparam


def references_definition():
    # Initialize Reference Frame Defintions Holder
    references_def = py_ref_def.t_references_def()

    # Hub Frame
    Hub = py_ref_def.t_frame_def()
    Hub.Name = "Hub"
    Hub.Parent = "Root"
    Hub.origin = np.array([0.0, 0.0, 0.0])
    Hub.orientation = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    Hub.moving = False

    # Append to References
    references_def.append_frame(Hub)

    # # Hub Frame
    # Hub1 = py_ref_def.t_frame_def()
    # Hub1.Name = "Hub1"
    # Hub1.Parent = "Root"
    # Hub1.origin = np.array([0.0, 0.0, 0.0])
    # Hub1.orientation = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    # Hub1.moving = False

    # # Append to References
    # references_def.append_frame(Hub1)

    return references_def


# rst ref (end)

# rst geodef
def geometry_definition():
    # Load rotor data
    f = open('ThesisCode/optimisation/NASA_tiltwing/helix_dir/prop_data/rotor.json')
    airfoilSecs = json.load(f)
    f.close()
    
    # Initialize Geometry Component Definitions Holder
    geometry_def = geo_def.t_geometry_def()

    # ==========================================================================
    # Define Rotor
    # ==========================================================================
    # ---------------------------- Blade Parameters -------------------------- #
    rotor = geo_def_parametric.t_geometry_def_parametric()
    rotor.CompName = "rotor"
    rotor.CompType = "rotor"
    rotor.RefName = "Hub"

    # Reference Parameters
    N_span = len(airfoilSecs["span"])
    rotor.ref_point = airfoilSecs["ref_point"]
    rotor.ref_chord_frac = 0.5

    # Symmetry Parameters
    rotor.symmetry = False
    rotor.symmetry_point = np.array([0.0, 0.0, 0.0])
    rotor.symmetry_normal = np.array([0.0, 1.0, 0.0])

    # Mirror Parameters
    rotor.mirror = True
    rotor.mirror_point = np.array([0.0, 0.0, 0.0])
    rotor.mirror_normal = np.array([0.0, 1.0, 0.0])

    # Initialize Rotor and Allocate Arrays
    rotor.initialize_parametric_geometry_definition(N_span)

    omega = 81.828 # rad/s
    nBlades = 5
    collective = np.deg2rad(37.0)

    rotor.multiple = True
    rotor.multiplicity = {
        "mult_type": "rotor",
        "n_blades": nBlades,
        "rot_axis": np.array([-1.0, 0.0, 0.0]),
        "rot_rate": omega,
        "psi_0": 0.0,
        "hub_offset": 0.0,
        "n_dofs": 1,
        "dofs": [
            {
                "hinge_type": "pitch",
                "hinge_offset": np.array([0.0, 0.0, 0.0]),
                "collective": collective,
                "cyclic_sine_amplitude": 0.0,
                "cyclic_cosine_amplitude": 0.0,
                "cyclic_sine_phase": 0.0,
                "cyclic_cosine_phase": 0.0,
            },
        ],
    }
    
    # ------------------------ Blade Section Definition ---------------------- #
    for iSec in range(len(airfoilSecs['chord'])):
        # Chord ------------------
        rotor.sec[iSec].chord = airfoilSecs['chord'][iSec]
        rotor.sec[iSec].twist = airfoilSecs['theta'][iSec]
        rotor.sec[iSec].thick = airfoilSecs['thick'][iSec]
        rotor.sec[iSec].alpha_0 = airfoilSecs['alpha_0'][iSec]
        rotor.sec[iSec].alpha_L0 = airfoilSecs['alpha_L0'][iSec]
        rotor.sec[iSec].Cl_alpha = airfoilSecs['Cl_alpha'][iSec]
        rotor.sec[iSec].M = airfoilSecs['M'][iSec]

    for iSpan in range(len(airfoilSecs["span"])):
        rotor.span[iSpan].span = airfoilSecs["span"][iSpan]
        rotor.span[iSpan].sweep = 0.0
        rotor.span[iSpan].dihed = 0.0
        rotor.span[iSpan].N_elem_span = 1
        rotor.span[iSpan].span_type = 1

    # Append To Vehicle
    geometry_def.append_component(rotor)

    return geometry_def
