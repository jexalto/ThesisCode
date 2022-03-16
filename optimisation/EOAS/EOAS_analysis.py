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

class surface(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('panels_span', default=71, desc='number of spanwise panels on the VLM')
        self.options.declare('panels_chord', default=71, desc='number of chordwise panels on the VLM')

    def setup(self):
        panels_span = self.options('panels_span')
        panels_chord = self.options('panels_chord')

        self.add_input('mesh', val=np.zeros((panels_chord, panels_span, 3)))
        self.add_input('correction', val=np.zeros((panels_span, panels_span)))
        self.add_input('correction_loc', val=np.zeros((panels_span)))

        self.add_input('chord', val=1.0, units='m')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('jet_radius', val=1.0, units='m')
        self.add_input('vinf', val=1.0, units='m/s')
        self.add_input('vjet', val=1.0, units='m/s')
        self.add_input('aoa', val=1.0, units='deg')

        self.add_input("Mach_number", val=0.84)
        self.add_input("re", val=1.0e6, units="1/m")
        self.add_input("rho", val=0.38, units="kg/m**3")
        self.add_input("CT", val=grav_constant * 17.0e-6, units="1/s")
        self.add_input("R", val=11.165e6, units="m")
        self.add_input("W0", val=0.4 * 3e5, units="kg")
        self.add_input("speed_of_sound", val=295.4, units="m/s")
        self.add_input("load_factor", val=1.)
        self.add_input("empty_cg", val=np.zeros((3)), units="m")

        self.add_output('mesh', val=np.zeros((panels_chord, panels_span, 3)))
        self.add_output('correction', val=np.zeros((panels_span, panels_span)))
        self.add_output('correction_loc', val=np.zeros((panels_span)))


    def compute(self, inputs, outputs):
        outputs['mesh'] = inputs['mesh']
        outputs['correction'] = inputs['correction']
        outputs['correction_loc'] = inputs['correction_loc']

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
            "mesh": inputs['mesh'],
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
                "vjet",
                "r0",
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
        self.connect(name + ".local_stiff_transformed", point_name + ".coupled." + name + ".local_stiff_transformed")
        self.connect(name + ".nodes", point_name + ".coupled." + name + ".nodes")

        # Connect aerodyamic mesh to coupled group mesh
        self.connect(name + ".mesh", point_name + ".coupled." + name + ".mesh")

        # Connect performance calculation variables
        self.connect(name + ".radius", com_name + ".radius")
        self.connect(name + ".thickness", com_name + ".thickness")
        self.connect(name + ".nodes", com_name + ".nodes")
        self.connect(name + ".cg_location", point_name + "." + "total_perf." + name + "_cg_location")
        self.connect(name + ".structural_mass", point_name + "." + "total_perf." + name + "_structural_mass")
        self.connect(name + ".t_over_c", com_name + ".t_over_c")

        # Set up the problem
        self.setup()

        om.n2()

        # Choose the angle of attack (alpha) values to sweep through
        alpha = np.array(2.0)

            # Set the alpha in the problem and run analysis
        prob["alpha"] = alpha
        prob.run_model()

        print("\nAngle of attack:", prob["alpha"])
        print("CL:", prob["AS_point_0.wing_perf.CL"])
        print("CD:", prob["AS_point_0.wing_perf.CD"])

        CL = prob["AS_point_0.wing_perf.CL"]
        CD = prob["AS_point_0.wing_perf.CD"]

        cl = prob.get_val('AS_point_0.wing_perf.aero_funcs.Cl')
        cd = prob.get_val('AS_point_0.wing_perf.aero_funcs.Cd')

        Wstruc = prob["wing.structural_mass"][0]
        Wfuel = prob["AS_point_0.fuelburn"][0]
