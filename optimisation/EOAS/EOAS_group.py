from bezier import Surface
import numpy as np
import openmdao.api as om

from parametersinput import parameters
from openaerostruct.utils.constants import grav_constant
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class EOAS(om.Group):

    def initialize(self):
        self.options.declare('panels_span_VLM', default=100)
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
            "root_chord": 1.0,  # root chord
            # "dihedral": 0.0,  # wing dihedral angle in degrees
            # positive is upward
            # "sweep": 0.0,  # wing sweep angle in degrees
            # positive sweeps back
            "taper": 1.0,  # taper ratio; 1. is uniform chord
            "num_twist_cp": 5
        }

        mesh = generate_mesh(mesh_dict) # twist_cp
        chord_cp = np.ones((1))
        twist_cp = np.zeros((5))

        data_x_upper = np.array([0.000000,  0.005000, 0.007500, 0.012500, 0.025000, 0.050000, 0.075000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.750000, 0.800000, 0.850000, 0.900000, 0.950000, 1.000000])

        data_y_upper = np.array([0.000000, 0.011930, 0.014360, 0.018150, 0.025080, 0.034770, 0.042020, 0.047990, 0.057320, 0.064230, 0.069260, 0.072700, 0.074630, 0.074870, 0.073130, 0.069780, 0.065170, 0.059560, 0.053110, 0.046000, 0.038470, 0.030840, 0.023210, 0.015580, 0.007950, 0.000320])

        data_y_lower = np.array([0.000000, -0.011930, -0.014360, -0.018150, -0.025080, -0.034770, -0.042020, -0.047990, -0.057320, -0.064230, -0.069260, -0.072700, -0.074630, -0.074870, -0.073130, -0.069780, -0.065170, -0.059560, -0.053110, -0.046000, -0.038470, -0.030840, -0.023210, -0.015580, -0.007950, -0.000320])

        surface = {
            # Wing definition
            "name": "wing",  # name of the surface
            "symmetry": False,  # if true, model one half of wing
            # reflected across the plane y = 0
            "S_ref_type": "wetted",  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            "fem_model_type": "wingbox",
            "spar_thickness_cp": np.array([0.004, 0.005, 0.005, 0.008, 0.008, 0.01]),  # [m]
            "skin_thickness_cp": np.array([0.005, 0.01, 0.015, 0.020, 0.025, 0.026]),
            "original_wingbox_airfoil_t_over_c": 0.12,
            "data_x_upper":data_x_upper,
            "data_x_lower":data_x_upper,
            "data_y_upper":data_y_upper,
            "data_y_lower":data_y_lower,
            "strength_factor_for_upper_skin": 1.0,
            # "twist_cp": twist_cp,
            "mesh": mesh,
            "span": span_0,
            "chord_cp": chord_cp,
            "propeller": 1,
            # "n_point_masses": 1,
            "radii_shape": radii_shape,
            "twist_cp": twist_cp,
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
            "with_viscous": True,
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
        
        aerostruct_group = AerostructGeometry(surface=surface)

        name = "wing"

        # Add tmp_group to the problem with the name of the surface.
        self.add_subsystem(name, aerostruct_group)

        point_name = "AS_point_0"

        # Create the aero point group and add it to the model
        AS_point = AerostructPoint(surfaces=[surface])

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
                "CT",
                "R",
                "W0",
                "speed_of_sound",
                "empty_cg",
                "load_factor",
                'correction',
            ],
        )

        self.connect('wing.mesh', 'AS_point_0.coupled.mesh')

        self.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
        self.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        self.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        point_name = 'AS_point_0'
        com_name = point_name + "." + name + "_perf"
        # self.connect(name + ".radius", com_name + ".radius")
        # self.connect(name + ".nodes", com_name + ".nodes")
        self.connect(name + ".t_over_c", com_name + ".t_over_c")

        self.connect(name + ".nodes", com_name + ".nodes")
        self.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        self.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")

        # Connect wingbox properties to von Mises stress calcs
        self.connect(name + ".Qz", com_name + ".Qz")
        self.connect(name + ".J", com_name + ".J")
        self.connect(name + ".A_enc", com_name + ".A_enc")
        self.connect(name + ".htop", com_name + ".htop")
        self.connect(name + ".hbottom", com_name + ".hbottom")
        self.connect(name + ".hfront", com_name + ".hfront")
        self.connect(name + ".hrear", com_name + ".hrear")

        self.connect(name + ".spar_thickness", com_name + ".spar_thickness")
        # self.connect(name + ".t_over_c", com_name + ".t_over_c")
