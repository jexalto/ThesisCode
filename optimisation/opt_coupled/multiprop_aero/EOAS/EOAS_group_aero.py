from bezier import Surface
import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint

class EOAS(om.Group):

    def initialize(self):
        self.options.declare('panels_span_VLM', default=300)
        self.options.declare('panels_chord_VLM', default=1)
        self.options.declare('span_0', default=10.)
        self.options.declare('radii_shape', default=20)

    def setup(self):
        span_0 = self.options['span_0']
        radii_shape = self.options['radii_shape']

        mesh_dict = {
            # Wing definition
            "num_x": self.options['panels_chord_VLM']+1,  # number of chordwise points
            "num_y": self.options['panels_span_VLM']+1,  # number of spanwise points --> NEEDS to be more than 11
            "wing_type": "rect",  # initial shape of the wing
            # either 'CRM' or 'rect'
            # 'CRM' can have different options
            # after it, such as 'CRM:alpha_2.75'
            # for the CRM shape at alpha=2.75
            "symmetry": False,  # if true, model one half of wing
            # reflected across the plane y = 0
            # Simple Geometric Variables
            "span": span_0,  # full wingspan, even for symmetric cases
            "chord": 0.15,
            # "dihedral": 0.0,  # wing dihedral angle in degrees
            # positive is upward
            # "sweep": 0.0,  # wing sweep angle in degrees
            # positive sweeps back
            # "taper": 1.0,  # taper ratio; 1. is uniform chord
            "num_twist_cp": 10
        }

        mesh = generate_mesh(mesh_dict)

        spanwisesections = 7

        surface = {
            # Wing definition
            "name": "wing",  # name of the surface
            "symmetry": False,  # if true, model one half of wing
            # "type": "aero",
            # reflected across the plane y = 0
            "S_ref_type": "projected",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "fem_model_type": "tube",
            "thickness_cp": np.array([0.001, 0.002, 0.003, 0.002, 0.001])*10,
            # "data_x_upper": upper_x,
            # "data_x_lower": lower_x,
            # "data_y_upper": upper_y,
            # "data_y_lower": lower_y,
            # "thickness_cp": np.array([0.1, 0.2, 0.3]),
            "spar_thickness_cp": np.array([0.004, 0.005, 0.005, 0.008, 0.008, 0.01]),  # [m]
            "skin_thickness_cp": np.array([0.005, 0.01, 0.015, 0.020, 0.025, 0.026]),
            "original_wingbox_airfoil_t_over_c": 0.12,
            # "strength_factor_for_upper_skin": 1.0, 
            "mesh": mesh,
            "span": span_0,
            "chord_cp": np.ones(spanwisesections)*0.24,  # Define chord using 3 B-spline cp's
            "twist_cp": np.zeros(spanwisesections),
            "propeller": 2,
            "n_point_masses": 2,
            "radii_shape": radii_shape,
            "electric": 0,
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            "CL0": 0.0,  # CL of the surface at alpha=0
            "CD0": 0.001,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.15,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c_cp": np.array([0.25]),  # thickness over chord ratio (NACA0015)
            "c_max_t": 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            "with_viscous": True,
            "electric": True,
            "with_wave": False,  # if true, compute wave drag
            # Structural values are based on aluminum 7075
            "E": 70.0e9,  # [Pa] Young's modulus of the spar
            "G": 30.0e9,  # [Pa] shear modulus of the spar
            "yield": 500.0e6 / 2.5,  # [Pa] yield stress divided by 2.5 for limiting case
            "mrho": 3.0e3,  # [kg/m^3] material density
            "fem_origin": 0.35,  # normalized chordwise location of the spar
            "wing_weight_ratio": 2.0,
            "struct_weight_relief": True,  # True to add the weight of the structure to the loads on the structure
            "distributed_fuel_weight": False,
            # Constraints
            "exact_failure_constraint": False,  # if false, use KS function
        }

        name = "wing"

        geom_group = Geometry(surface=surface)
        self.add_subsystem(name, geom_group)

        point_name = "AS_point_0"

        # Create the aero point group and add it to the model
        AS_point = AeroPoint(surfaces=[surface])

        self.add_subsystem(
            point_name,
            AS_point,
            promotes_inputs=[
                "v",
                "velocity_distr",
                "alpha",
                "Mach_number",
                "re",
                "rho",
                'correction',
            ],
        )

        self.connect(name + ".mesh", point_name + "." + name + ".def_mesh")

        # Perform the connections with the modified names within the
        # 'aero_states' group.
        self.connect(name + ".mesh", point_name + ".aero_states." + name + "_def_mesh")

        self.connect(name + ".t_over_c", point_name + "." + name + "_perf." + "t_over_c")