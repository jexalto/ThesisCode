# rst imports
# Helix Libraries
import helix.model

import helix.parameters.simparam_def as simparam_def

import helix.geometry.geometry_def as geo_def
import helix.geometry.geometry_def_parametric as geo_def_parametric

import helix.references.references_def as ref_def

import helix_pyf90.mod_post_output_vtk as output_vtk

# External
from mpi4py import MPI
import numpy as np

import matplotlib.pyplot as plt

import json

# rst imports (end)

###############################################################################
## Main - Run Script
###############################################################################
# rst main
def main():
    # ---------------------- Assign Simulation Parameters ------------------- #
    simparam = set_simparam()

    # -------------------------- Set Reference Frames ----------------------- #
    references_def = set_references()

    # ----------------------------- Build Vehicle --------------------------- #
    geometry_def = geometry_definition()

    # ------------------------------- Run HELIX ----------------------------- #
    model = helix.model.Model(thrust=True, torque=True, moment=True, loads=True)
    model.simparam_def = simparam
    model.geometry_def = geometry_def
    model.references_def = references_def
    model.comm = MPI.COMM_WORLD
    model.initialize_run()

    model.run()
    # output_vtk.write_vtk(model.f_geometry,1,True,False,"test")

    thrust_vec = model.f_simparam.t_vec

    print(f'CT: {model.f_geometry.rotor[0].ct}')

    plt.plot(np.linspace(0, len(thrust_vec), len(thrust_vec)), thrust_vec)
    plt.grid()
    plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/figures/thrustvector.png')


# rst main (end)

###############################################################################
## Set Simulation Parameters
###############################################################################
# rst simparam
def set_simparam():
    simparam = simparam_def.t_simparam_def()
    simparam.basename = "MR8x45 Rotor"

    simparam.dt = 0.01
    simparam.t_start = 0.0
    simparam.t_end = 0.02

    simparam.nt_rev = 30

    simparam.v_inf = np.array([0., 0.0, -40.])
    simparam.rho_inf = 1.25

    return simparam


# rst simparam (end)

###############################################################################
## Set References
###############################################################################
# rst ref
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
    Hub.orientation = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    Hub.moving = False

    # Append to References
    references_def.append_frame(Hub)

    return references_def


# rst ref (end)

###############################################################################
## Generate Vehicle
###############################################################################
# rst geodef
def geometry_definition():
    # Load rotor data
    f = open('/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/data/rotor.json')
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
    rotor.ref_point = np.array([0.0, 0.02023364, 0.0])
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
    J = 0.7
    V = 40
    diameter = 0.237

    n = -(V/(J*diameter))*2*np.pi
    rotor.multiple = True
    rotor.multiplicity = {
        "mult_type": "rotor",
        "n_blades": 4,
        "rot_axis": np.array([0.0, 0.0, 1.0]),
        "rot_rate": n,
        "psi_0": 0.0,
        "hub_offset": 0.,
        "n_dofs": 0,
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
        rotor.sec[iSec].M = 50.0

    for iSpan in range(len(airfoilSecs["span"])):
        # input(airfoilSecs["span"][iSpan])
        rotor.span[iSpan].span = airfoilSecs["span"][iSpan]
        rotor.span[iSpan].sweep = 0.0
        rotor.span[iSpan].dihed = 0.0
        rotor.span[iSpan].N_elem_span = 4
        rotor.span[iSpan].span_type = 1

    # Append To Vehicle
    geometry_def.append_component(rotor)

    return geometry_def


# rst geodef (end)

###############################################################################
## Initializer
###############################################################################
if __name__ == "__main__":
    main()
