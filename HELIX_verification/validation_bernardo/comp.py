__author__ = "Bernardo Pacini"
__email__ = "bpacini@umich.edu"
__date__ = "May 17th, 2021 Mon"
__status__ = "Production"

# External Imports
import json
import numpy as np

# Helix Libraries
import helix.parameters.simparam_def as helix_simparam_def
import helix.geometry.geometry_def_parametric as helix_geo_def_parametric
import helix.references.references_def as helix_ref_def

# Pulse Libraries
import pulse.parameters.simparam_def as pulse_simparam_def
import pulse.geometry.geometry_def_parametric as pulse_geo_def_parametric
import pulse.references.references_def as pulse_ref_def


def simparamShared(simparam):
    """
    This function sets simulation parameters that are common to both HELIX
    and PULSE analyses.
    Parameters
    ----------
    simparam : t_simparam_def (HELIX Python OR PULSE Python)
        Simulation parameters definition holder for PULSE Python wrapper OR
        HELIX Python wrapper (the function identifies which, if neither it will
        raise an error). The simparam object is modified in place.
    """
    if isinstance(simparam, helix_simparam_def.t_simparam_def):
        tool = "HELIX"
    elif isinstance(simparam, pulse_simparam_def.t_simparam_def):
        tool = "PULSE"
    else:
        raise TypeError("Simparam must be HELIX or PULSE simparam type")

    simparam.basename = "NASA_QuadE_1_RPM_" + tool

    simparam.v_inf = np.array([0.0, 0.0, 0.0])  # [m/s]
    simparam.rho_inf = 0.9848  # [kg/m3]


def geoDefShared(rotor):
    """
    This function sets geometry definition parameters that are common to both
    HELIX and PULSE analyses. Specifically, the function sets the geometry
    defintion for one parameterically defined rotor.
    Parameters
    ----------
    rotor : list (t_geometry_def_parametric (HELIX Python OR PULSE Python))
        A list of parametric geometry definition holders for HELIX Python
        wrapper OR PULSE Python wrapper (the function identifies which, if
        neither it will raise an error). The list is modified in place.
    """
    if isinstance(rotor[0], helix_geo_def_parametric.t_geometry_def_parametric):
        tool = "HELIX"
    elif isinstance(rotor[0], pulse_geo_def_parametric.t_geometry_def_parametric):
        tool = "PULSE"
    else:
        raise TypeError("Rotor must contain all HELIX or all PULSE parametric geometry component definitions.")

    # Temporary
    rpm = 716.197

    # Read Geometry JSON file
    fileName = "../rotor.json"
    with open(fileName, "r") as file:
        airfoilData = json.load(file)

    nSpan = np.size(airfoilData["span"])
    nBlades = 3

    # =========================================================================
    # Define Rotor
    # =========================================================================
    # --------------------------- Blade Parameters -------------------------- #
    rotor[0].CompName = "Rotor_1"
    rotor[0].CompType = "rotor"
    rotor[0].RefName = "Hub1"
    if tool == "PULSE":
        rotor[0].Compact = True

    # Reference Parameters
    rotor[0].ref_point = airfoilData["ref_point"]
    rotor[0].ref_chord_frac = 0.5

    # Symmetry Parameters
    rotor[0].symmetry = False
    rotor[0].symmetry_point = np.array([0.0, 0.0, 0.0])
    rotor[0].symmetry_normal = np.array([0.0, 1.0, 0.0])

    # Mirror Parameters
    rotor[0].mirror = False
    rotor[0].mirror_point = np.array([0.0, 0.0, 0.0])
    rotor[0].mirror_normal = np.array([0.0, 1.0, 0.0])

    # Initialize Rotor and Allocate Arrays
    if tool == "HELIX":
        rotor[0].initialize_parametric_geometry_definition(nSpan)
    elif tool == "PULSE":
        rotor[0].initialize_parametric_geometry_definition(nSpan, nBlades)

    rotor[0].multiple = True
    rotor[0].multiplicity = {
        "mult_type": "rotor",
        "n_blades": nBlades,
        "rot_axis": np.array([0.0, 0.0, 1.0]),
        "rot_rate": rpm / 60 * 2.0 * np.pi,
        "psi_0": np.pi,
        "hub_offset": 0.0,
        "n_dofs": 1,
        "dofs": [
            {
                "hinge_type": "pitch",
                "hinge_offset": np.array([0.0, 0.0, 0.0]),
                "collective": np.deg2rad(7.676),
                "cyclic_sine_amplitude": 0.0,
                "cyclic_cosine_amplitude": 0.0,
                "cyclic_sine_phase": 0.0,
                "cyclic_cosine_phase": 0.0,
            },
        ],
    }

    # ----------------------- Blade Section Definition ---------------------- #
    # Chord Sections
    for iSection in range(len(airfoilData["chord"])):
        rotor[0].sec[iSection].chord = airfoilData["chord"][iSection]
        rotor[0].sec[iSection].twist = airfoilData["twist"][iSection]

        if tool == "HELIX":
            rotor[0].sec[iSection].alpha_0 = airfoilData["alpha_0"][iSection]
            rotor[0].sec[iSection].alpha_L0 = airfoilData["alpha_L0"][iSection]
            rotor[0].sec[iSection].Cl_alpha = airfoilData["Cl_alpha"][iSection]
            rotor[0].sec[iSection].M = airfoilData["M"][iSection]

        elif tool == "PULSE":
            rotor[0].sec[iSection].thick = airfoilData["thick"][iSection]

    # Span Sections
    for iSpan in range(len(airfoilData["span"])):
        rotor[0].span[iSpan].span = airfoilData["span"][iSpan]
        rotor[0].span[iSpan].sweep = 0.0
        rotor[0].span[iSpan].dihed = 0.0
        rotor[0].span[iSpan].N_elem_span = 4
        rotor[0].span[iSpan].span_type = 1


def refDefShared(hub):
    """
    This function sets references parameters that are common to both HELIX and
    PULSE analyses. Specifically, the function sets the frame definitions of
    one reference frame.
    Parameters
    ----------
    hub : list (t_frame_def (HELIX Python OR PULSE Python))
        A list of reference frame definition holders for HELIX Python
        wrapper OR PULSE Python wrapper (the function identifies which, if
        neither it will raise an error). The list is modified in place.
    """
    if isinstance(hub[0], helix_ref_def.t_frame_def):
        tool = "HELIX"  # NOQA F841
    elif isinstance(hub[0], pulse_ref_def.t_frame_def):
        tool = "PULSE"  # NOQA F841
    else:
        raise TypeError("Hub must contain all HELIX or all PULSE reference frame definitions")

    # =========================================================================
    # Hub 1 Frame (Front Left)
    # =========================================================================
    hub[0].Name = "Hub1"
    hub[0].Parent = "Root"
    hub[0].origin = np.array([0.0, 0.0, 0.0])
    hub[0].orientation = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    hub[0].moving = False