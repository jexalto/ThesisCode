"""
This script is an example of a black box setup using OpenAeroStruct.
We first construct and setup a problem, then manually change values in the
prob instance and run analysis at different points. Through this method,
we can manually explore the design space.

Although this example varies angle of attack (alpha), you could vary any input
into the problem, including wing geometry and flight conditions.
"""

import numpy as np
import pandas as pd

from openaerostruct.geometry.utils import generate_mesh

from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

import openmdao.api as om
from openaerostruct.utils.constants import grav_constant

from jet_correction_ib import jet_correction_ib

import matplotlib.pyplot as plt

# Create a dictionary to store options about the surface
mesh_dict = {
        # Wing definition
        "num_x": 5,  # number of chordwise points
        "num_y": 29,  # number of spanwise points --> NEEDS to be more than 11
        "wing_type": "rect",  # initial shape of the wing
        # either 'CRM' or 'rect'
        # 'CRM' can have different options
        # after it, such as 'CRM:alpha_2.75'
        # for the CRM shape at alpha=2.75
        "symmetry": True,  # if true, model one half of wing
        # reflected across the plane y = 0
        # Simple Geometric Variables
        "span": 10.0,  # full wingspan, even for symmetric cases
        "root_chord": 1.0,  # root chord
        "dihedral": 0.0,  # wing dihedral angle in degrees
        # positive is upward
        "sweep": 0.0,  # wing sweep angle in degrees
        # positive sweeps back
        "taper": 1.0,  # taper ratio; 1. is uniform chord
        "num_twist_cp": 5
    }

mesh = generate_mesh(mesh_dict) # twist_cp

# --- Correction factor calculation ---
# --- Include error for too little mesh points ---
# y = np.linspace(0, mesh_dict["span"]/2, int((mesh_dict["num_y"])) )
y = np.linspace(0, 10, 15 ) # np.linspace(0, mesh_dict["span"]/2, int((mesh_dict["num_y"])) )
mu = 0.95
crd = 1.0
loc = 5-0.25 # 0.5*mesh_dict["span"]/2 - 0.25
ib_loc = 5 # 0.5*mesh_dict["span"]/2

if mesh_dict["num_x"] == 1:
    n = 1
else:
    n = int((mesh_dict["num_x"]-3)/2 + 1)

G = jet_correction_ib(y, mu, crd, loc, ib_loc) # np.zeros((mesh_dict["num_y"]-1, mesh_dict["num_y"]-1)) # 
G = np.tile(G, (n, 1))
G = np.tile(G, (1, n))

quit()
# -------------------------------------

surface = {
    # Wing definition
    "name": "wing",  # name of the surface
    "symmetry": True,  # if true, model one half of wing
    # reflected across the plane y = 0
    "S_ref_type": "wetted",  # how we compute the wing area,
    # can be 'wetted' or 'projected'
    "fem_model_type": "tube",
    "thickness_cp": np.array([0.1, 0.2, 0.3]),
    # "twist_cp": twist_cp,
    "mesh": mesh,
    # Aerodynamic performance of the lifting surface at
    # an angle of attack of 0 (alpha=0).
    # These CL0 and CD0 values are added to the CL and CD
    # obtained from aerodynamic analysis of the surface to get
    # the total CL and CD.
    # These CL0 and CD0 values do not vary wrt alpha.
    "CL0": 0.0,  # CL of the surface at alpha=0
    "CD0": 0.0,  # CD of the surface at alpha=0
    # Airfoil properties for viscous drag calculation
    "k_lam": 0.05,  # percentage of chord with laminar
    # flow, used for viscous drag
    "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
    "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
    # thickness
    "with_viscous": False,
    "with_wave": False,  # if true, compute wave drag
    # Structural values are based on aluminum 7075
    "E": 70.0e9,  # [Pa] Young's modulus of the spar
    "G": 30.0e9,  # [Pa] shear modulus of the spar
    "yield": 500.0e6 / 2.5,  # [Pa] yield stress divided by 2.5 for limiting case
    "mrho": 3.0e3,  # [kg/m^3] material density
    "fem_origin": 0.35,  # normalized chordwise location of the spar
    "wing_weight_ratio": 2.0,
    "correction": G,
    "struct_weight_relief": False,  # True to add the weight of the structure to the loads on the structure
    "distributed_fuel_weight": False,
    # Constraints
    "exact_failure_constraint": False,  # if false, use KS function
}

# Create the problem and assign the model group
prob = om.Problem()

# Add problem information as an independent variables component
indep_var_comp = om.IndepVarComp()
indep_var_comp.add_output("v", val=100, units="m/s")
indep_var_comp.add_output("alpha", val=1.0, units="deg")
indep_var_comp.add_output("Mach_number", val=0.)
indep_var_comp.add_output("re", val=1.0e6, units="1/m")
indep_var_comp.add_output("rho", val=1.225, units="kg/m**3")
indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
indep_var_comp.add_output("R", val=11.165e6, units="m")
indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
indep_var_comp.add_output("load_factor", val=1.0)
indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")

prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

aerostruct_group = AerostructGeometry(surface=surface)

name = "wing"

# Add tmp_group to the problem with the name of the surface.
prob.model.add_subsystem(name, aerostruct_group)

point_name = "AS_point_0"

# Create the aero point group and add it to the model
AS_point = AerostructPoint(surfaces=[surface])

prob.model.add_subsystem(
    point_name,
    AS_point,
    promotes_inputs=[
        "v",
        "alpha",
        "Mach_number",
        "re",
        "rho",
        "CT",
        "R",
        "W0",
        "speed_of_sound",
        "empty_cg",
        "load_factor",
    ],
)

com_name = point_name + "." + name + "_perf"
prob.model.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

# Connect aerodyamic mesh to coupled group mesh
prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

# Connect performance calculation variables
prob.model.connect(name + ".radius", com_name + ".radius")
prob.model.connect(name + ".thickness", com_name + ".thickness")
prob.model.connect(name + ".nodes", com_name + ".nodes")
prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

# Set up the problem
prob.setup()

# Choose the angle of attack (alpha) values to sweep through
alphas = np.array([1.0]) # np.linspace(-5.0, 5.0, 11)

# Loopo through each alpha
for alpha in alphas:

    # Set the alpha in the problem and run analysis
    prob["alpha"] = alpha
    prob.run_model()

    print()
    print("Angle of attack:", prob["alpha"])
    print("CL:", prob["AS_point_0.wing_perf.CL"])
    print("CD:", prob["AS_point_0.wing_perf.CD"])

totlift = np.concatenate( (prob['AS_point_0.wing_perf.aero_funcs.liftcoeff.Cl'], np.flip(prob['AS_point_0.wing_perf.aero_funcs.liftcoeff.Cl'])) )
plt.plot(np.linspace(-mesh_dict["span"]/2, mesh_dict["span"]/2, int((mesh_dict["num_y"]-1))), totlift, label="OAS data")
plt.grid()
# plt.xlim((-mesh_dict["span"]/2, mesh_dict["span"]/2))
# plt.ylim((min(totlift)*0.98, max(totlift)*1.2))

data, target = np.array_split(np.loadtxt('/home/jexalto/code/MDO_lab_env/ThesisCode/OAS_validation/vlm_verification.txt', dtype=float), [-1], axis=0)

y_vlm_verification = data[:, 0]
cl_vlm_verification = data[:, 1]

plt.plot(y_vlm_verification, cl_vlm_verification, label="VLM verification data")
plt.legend()
plt.savefig("/home/jexalto/code/MDO_lab_env/ThesisCode/OAS_validation/Figures/liftdistribution.png")