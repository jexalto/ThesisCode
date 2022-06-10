from bezier import Surface
import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

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

        upper_x = np.array([    0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4,
                                0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,],
            dtype="complex128",
        )
        lower_x = np.array([    0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43,
                                0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
            ],
            dtype="complex128",
        )
        upper_y = np.array([    0.0447, 0.046, 0.0472, 0.0484, 0.0495, 0.0505, 0.0514, 0.0523, 0.0531, 0.0538, 0.0545, 0.0551, 0.0557, 0.0563, 0.0568, 0.0573, 0.0577, 0.0581, 0.0585, 0.0588, 0.0591, 0.0593, 0.0595, 0.0597,
                                0.0599, 0.06, 0.0601, 0.0602, 0.0602, 0.0602, 0.0602, 0.0602, 0.0601, 0.06, 0.0599, 0.0598, 0.0596, 0.0594, 0.0592, 0.0589, 0.0586, 0.0583, 0.058, 0.0576, 0.0572, 0.0568, 0.0563, 0.0558, 0.0553, 0.0547, 0.0541,
            ],
            dtype="complex128",
        )
        lower_y = np.array([    -0.0447, -0.046, -0.0473, -0.0485, -0.0496, -0.0506, -0.0515, -0.0524, -0.0532, -0.054, -0.0547, -0.0554, -0.056, -0.0565, -0.057, -0.0575, -0.0579, -0.0583, -0.0586, -0.0589, -0.0592, -0.0594,
                                -0.0595, -0.0596, -0.0597, -0.0598, -0.0598, -0.0598, -0.0598, -0.0597, -0.0596, -0.0594, -0.0592, -0.0589, -0.0586, -0.0582, -0.0578, -0.0573, -0.0567, -0.0561, -0.0554, -0.0546, -0.0538, -0.0529, -0.0519, -0.0509, -0.0497, -0.0485, -0.0472, -0.0458, -0.0444,
            ],
            dtype="complex128",
        )   

        mesh = generate_mesh(mesh_dict) # twist_cp
        # chord_cp = np.ones((10))
        twist_cp = np.zeros((5))

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
            # "spar_thickness_cp": np.array([0.004, 0.005, 0.005, 0.008, 0.008, 0.01]),  # [m]
            # "skin_thickness_cp": np.array([0.005, 0.01, 0.015, 0.020, 0.025, 0.026]),
            # "original_wingbox_airfoil_t_over_c": 0.12,
            # "strength_factor_for_upper_skin": 1.0, 
            "mesh": mesh,
            "span": span_0,
            "chord_cp": np.ones(5)*0.15,  # Define chord using 3 B-spline cp's
            "twist_cp": np.zeros(5),
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
            "CD0": 0.0,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            "k_lam": 0.15,  # percentage of chord with laminar
            # flow, used for viscous drag
            "t_over_c_cp": np.array([0.15]),  # thickness over chord ratio (NACA0015)
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
        com_name = point_name + "." + name + "_perf"
        self.connect(name + ".radius", com_name + ".radius")
        self.connect(name + ".thickness", com_name + ".thickness")
        self.connect(name + ".nodes", com_name + ".nodes")
        self.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        self.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        self.connect(name + ".t_over_c", com_name + ".t_over_c")
        com_name = "AS_point_0" + "." + name + "_perf."

        # self.connect(name + ".Qz", com_name + "Qz")
        # self.connect(name + ".J", com_name + "J")
        # self.connect(name + ".A_enc", com_name + "A_enc")
        # self.connect(name + ".htop", com_name + "htop")
        # self.connect(name + ".hbottom", com_name + "hbottom")
        # self.connect(name + ".hfront", com_name + "hfront")
        # self.connect(name + ".hrear", com_name + "hrear")

        # self.connect(name + ".spar_thickness", com_name + "spar_thickness")
