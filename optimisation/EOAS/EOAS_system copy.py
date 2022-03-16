import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant
from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

class EOAS(om.Group):
    def initialize(self):

    def setup(self):
        self.options.declare('panels_VLM', default=71, desc='number of spanwise panels on the VLM')
        self.options.declare('nx', default=2, desc='number of chordwise panels')
        
        # panels_VLM = self.options['panels_VLM']
        # nx = self.options['nx']

        # self.add_input('mesh', shape=(nx, panels_VLM), val=np.zeros((nx, panels_VLM)))
        # self.add_input('correction_matrix', shape=(panels_VLM, panels_VLM), val=np.zeros((panels_VLM, panels_VLM)))
        # self.add_input('chord', val=1.0, units='m')
        # self.add_input('jet_loc', val=1.0, units='m')
        # self.add_input('jet_radius', val=1.0, units='m')
        # self.add_input('vinf', val=1.0, units='m/s')
        # self.add_input('vjet', val=1.0, units='m/s')
        # self.add_input('aoa', val=1.0, units='deg')

        # self.add_input("Mach_number", val=0.84)
        # self.add_input("re", val=1.0e6, units="1/m")
        # self.add_input("rho", val=0.38, units="kg/m**3")
        # self.add_input("CT", val=grav_constant * 17.0e-6, units="1/s")
        # self.add_input("R", val=11.165e6, units="m")
        # self.add_input("W0", val=0.4 * 3e5, units="kg")
        # self.add_input("speed_of_sound", val=295.4, units="m/s")
        # self.add_input("load_factor", val=1.)
        # self.add_input("empty_cg", val=np.zeros((3)), units="m")

        # self.add_output('CL', val=1.0)
        # self.add_output('CD', val=1.0)
        # self.add_output('cl', shape=(1, panels_VLM), val=np.zeros((1, panels_VLM)))
        # self.add_output('cd', shape=(1, panels_VLM), val=np.zeros((1, panels_VLM)))
        # self.add_output('Wfuel', val=1.0)
        # self.add_output('Wwing', val=1.0)

    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        span = inputs['span']
        jet_loc = inputs['jet_loc']
        jet_radius = inputs['jet_radius']
        Vinf = inputs['vinf']
        Vjet = inputs['vjet']
        mesh = inputs['mesh']
        correction = inputs['correction']
        alpha = inputs['aoa']
        chord = 1.#inputs['chord']

        panels_VLM = self.options['panels_VLM']
        panels_chord = self.options['nx']

        self.add_subsystem('correction_loc_array', )
        
        y = mesh[:, 0, 0]
        correction_loc = np.zeros((panels_VLM))

        for index in range(panels_VLM):
            correction_loc[index] = (y[index+1]+y[index])/2

        correction_loc = np.where(correction_loc>(jet_loc+jet_radius), 0, correction_loc)
        correction_loc = np.where(correction_loc<(jet_loc-jet_radius), 0, correction_loc)
        correction_loc = np.where(correction_loc==0, 0, 1)

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
            "correction": correction,
            "correction_loc": correction_loc,
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
        indep_var_comp.add_output("vjet", val=Vjet, units="m/s")
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
        # The connections are if statements, this value=that value if connected
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

        outputs['CL'], outputs['CD'] = CL, CD
        outputs['cl'], outputs['cd'] = cl, cd
        outputs['Wstruc'], outputs['Wfuel'] = Wstruc, Wfuel