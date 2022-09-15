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

import matplotlib.pyplot as plt

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

import openmdao.api as om

from correction_matrix import correction_, correction_location

def EOAS_system_(jet_radius, jet_loc_list, Vinf, r_min, span_max, filename, nx_input, prop_discr, span, nr_props):
    # Create a dictionary to store options about the surface
    mesh_dict = {
            # Wing definition
            "num_x": nx_input,  # number of chordwise points --> interesting, 5 points equates to two panels
            "num_y": 201,  # number of spanwise points --> NEEDS to be more than 11
            "wing_type": "rect",  # initial shape of the wing
            # either 'CRM' or 'rect'
            # 'CRM' can have different options
            # after it, such as 'CRM:alpha_2.75'
            # for the CRM shape at alpha=2.75
            "symmetry": False,  # if true, model one half of wing
            # reflected across the plane y = 0
            # Simple Geometric Variables
            "span": span,  # full wingspan, even for symmetric cases,
            "root_chord": 1.0,  # root chord
            "dihedral": 0.0,  # wing dihedral angle in degrees
            # positive is upward
            "sweep": 0.0,  # wing sweep angle in degrees
            # positive sweeps back
            "taper": 1.0,  # taper ratio; 1. is uniform chord
            "num_twist_cp": 5
        }
    
    panels_overset_wing = 1001
    panels_jet = 51

    steps = 100
    step = np.arange(0, steps, 1)

    def x2(steps):
        return 130-(0.05*steps-5)**2
    
    vel_distr_input = np.array([x2(step), x2(step)], order='F')
    radii_input = np.array(np.linspace(0.01, jet_radius, steps), order='F') # np.linspace(0.01, jet_radius, steps), 
    
    panels_span_VLM = mesh_dict['num_y']-1
    panels_chord_VLM = mesh_dict['num_x']-1
    nr_radii_input = len(step)
    
    correction, y_VLM, vel_vec = correction_(   panels_span_VLM, panels_chord_VLM,  panels_overset_wing, panels_jet, span, span_max, r_min, Vinf,
                                                vel_distr_input, radii_input, prop_discr, jet_loc_list, nr_props, nr_radii_input)
    print(vel_vec)
    y = y_VLM

    mesh = np.zeros((2, mesh_dict["num_y"], 3))
    mesh[0, :, 1] = y
    mesh[1, :, 1] = y
    mesh[1, :, 0] = np.ones((len(y)))
    mesh = generate_mesh(mesh_dict) # !¡

# you have to change the mesh file in fortran!

    surface = {
        # Wing definition
        "name": "wing",  # name of the surface
        "symmetry": False,  # if true, model one half of wing
        # reflected across the plane y = 0
        "S_ref_type": "wetted",  # how we compute the wing area,
        # can be 'wetted' or 'projected'
        "fem_model_type": "tube",
        "thickness_cp": np.array([0.1, 0.2, 0.3]),
        # "twist_cp": twist_cp,
        "mesh": mesh,
        # "span": 10.,
        "propeller": 1,
        "radii_shape": steps,
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
        "struct_weight_relief": False,  # True to add the weight of the structure to the loads on the structure
        "distributed_fuel_weight": False,
        # Constraints
        "exact_failure_constraint": False,  # if false, use KS function
    }

    # Create the problem and assign the model group
    prob = om.Problem()

    # Add problem information as an independent variables component
    indep_var_comp = om.IndepVarComp()
    indep_var_comp.add_output("v", val=Vinf, units="m/s")
    indep_var_comp.add_output("correction", val=correction)
    indep_var_comp.add_output("jet_loc", val=jet_loc_list, units='m')
    indep_var_comp.add_output("jet_radius", val=radii_input, units='m')
    indep_var_comp.add_output('velocity_distr', val=vel_vec, units='m/s')
    indep_var_comp.add_output("alpha", val=1.0, units="deg")
    indep_var_comp.add_output("Mach_number", val=0.84)
    indep_var_comp.add_output("re", val=1.0e6, units="1/m")
    indep_var_comp.add_output("rho", val=0.38, units="kg/m**3")
    indep_var_comp.add_output("CT", val=grav_constant * 17.0e-6, units="1/s")
    indep_var_comp.add_output("R", val=11.165e6, units="m")
    indep_var_comp.add_output("W0", val=0.4 * 3e5, units="kg")
    indep_var_comp.add_output("speed_of_sound", val=295.4, units="m/s")
    indep_var_comp.add_output("load_factor", val=1.0)
    indep_var_comp.add_output("empty_cg", val=np.zeros((3)), units="m")
    indep_var_comp.add_output("span", val=10., units="m")

    prob.model.add_subsystem("prob_vars", indep_var_comp, promotes=["*"])

    aerostruct_group = AerostructGeometry(surface=surface)

    name = "wing"

    # Add tmp_group to the problem with the name of the surface.
    prob.model.add_subsystem(name, aerostruct_group)

    point_name = "AS_point_0"        # Create the aero point group and add it to the model

    AS_point = AerostructPoint(surfaces=[surface])

    prob.model.add_subsystem(
        point_name,
        AS_point,
        promotes_inputs=[
            "v",
            "velocity_distr",
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
            'correction',
        ],
    )

    prob.model.connect('wing.mesh', 'AS_point_0.coupled.mesh')
    prob.model.connect('jet_loc', 'wing.geometry.mesh.jet_loc') #                               you need to reintroduce these!
    prob.model.connect('jet_radius', 'wing.geometry.mesh.jet_radius')

    prob.model.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
    prob.model.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

    # Connect aerodyamic mesh to coupled group mesh
    prob.model.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

    # Connect performance calculation variables
    com_name = point_name + "." + name + "_perf"
    prob.model.connect(name + ".radius", com_name + ".radius")
    prob.model.connect(name + ".thickness", com_name + ".thickness")
    prob.model.connect(name + ".nodes", com_name + ".nodes")
    prob.model.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
    prob.model.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
    prob.model.connect(name + ".t_over_c", com_name + ".t_over_c")

    # Set up the problem
    prob.setup()

    # Choose the angle of attack (alpha) values to sweep through
    alphas = np.array([2.0]) # np.linspace(-5.0, 5.0, 11)

    # Loopo through each alpha
    for alpha in alphas:

        # Set the alpha in the problem and run analysis
        prob["alpha"] = alpha
        prob.run_model()

        om.n2(prob)

        print()
        print("Angle of attack:", prob["alpha"])
        print("CL:", prob["AS_point_0.wing_perf.CL"])
        print("CD:", prob["AS_point_0.wing_perf.CD"])

    # --- Plotting ---
    data_cfd = pd.read_csv('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/data/cfd.txt')
    y_cfd = data_cfd['z[m]'].to_numpy()
    y_cfd *= -1
    beta_cfd = y_cfd/jet_radius
    cl_cfd = data_cfd[' Cl'].to_numpy()

    # print(prob.get_val('AS_point_0.wing_perf.CL'))
    cl_dist = prob.get_val('AS_point_0.wing_perf.aero_funcs.Cl')
    mesh_error = sum(y_VLM-prob['wing.mesh'][0, :, 1])
    print(f'mesh error: {mesh_error}')

    y_jet = jet_loc_list

    plt.clf()
    y = mesh[1, :, 1]
    y_ = np.zeros((len(y)-1))
    for i in range(len(y)-1):
        y_[i] = (y[i+1]+y[i])/2
    plt.plot(y_, cl_dist*1.06, label="Correction applied")
    # plt.plot(y_cfd, cl_cfd, label='CFD data')

    plt.scatter(y_jet, np.ones(len(jet_loc_list))*0.125, marker='x', color='b', label='Propeller')
    plt.scatter(y_jet-jet_radius, np.ones(len(jet_loc_list))*0.125, marker='|', color='b')
    plt.scatter(y_jet+jet_radius, np.ones(len(jet_loc_list))*0.125, marker='|', color='b')
    # plt.plot([y_jet-jet_radius, y_jet, y_jet+jet_radius], np.ones((3, 1))*0.125, color='b')
    plt.legend()
    plt.grid()
    plt.xlabel("Wingspan [m]")
    plt.ylabel("Lift coefficient [Cl]")
    # plt.title(f'VLM Panels={mesh_dict["num_y"]-1}-{nx_input-1}, Radius{jet_radius}')
    
    name_ = f'liftdistribution_r{jet_radius}_d{jet_loc_list[0]}.png'
    plt.savefig(filename+name_)       # this one is for the movie file
    # plt.savefig(filename+f"panels/liftdistribution_r{jet_radius}_d{jet_loc}_p{panels_overset_wing}.png")
    plt.clf()
    print(name_)
    return filename+name_