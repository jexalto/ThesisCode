import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix.openmdao.om_helix import HELIX_Group
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from helix_dir.helix_config import geometry
from rethorst_dir.rethorst_group import Rethorst
from constraints.cons_proploc import proploc
from obj_function import obj_func

class master(om.Problem):

    def setup(self):
        self.add_subsystem('parameters', subsys=parameters())

        simparam_def = simparam_definition()
        references_def = references_definition()
        geometry_def = geometry_definition()
        self.add_subsystem('helix', subsys=HELIX_Group(
                                                        simparam_def=simparam_def,
                                                        references_def=references_def,
                                                        geometry_def=geometry_def,
                                                        thrust_calc=True,
                                                        torque_calc=True,
                                                        moment_calc=True,
                                                        loads_calc=True,
                                                        velocity_distribution_calc=True,
                                                        ),
        )

        self.add_subsystem('rethorst', subsys=Rethorst(panels_VLM=301))
        
        self.add_subsystem('EOAS', subsys=EOAS())

        self.add_subsystem('obj_fun', subsys=obj_func())
        
        self.add_subsystem('cons_proploc', subsys=proploc())

        self.connect('rethorst.velocity_distribution', 'helix.rotorcomp_0_velocitydistribution')
        self.connect('rethorst.radii', 'helix.rotorcomp_0_radii')
        
        self.connect('EOAS.span', 'rethorst.span')

        self.connect('EOAS.correction_matrix', 'rethorst.correction_matrix')
        self.connect('EOAS.velocity_vector', 'rethorst.velocity_vector')
        
        self.connect('jet_loc', 'EOAS.wing.geometry.mesh.jet_loc')
        self.connect('radii', 'EOAS.wing.geometry.mesh.jet_radius')

prob = om.Problem()
model = prob.model

prob.model = master()

# --- Connecting parameters ----
model.connect('EOAS.v', 'parameters.vinf')
model.connect('EOAS.alpha', 'parameters.alpha')
model.connect('EOAS.Mach_number', 'parameters.Mach_number')
model.connect('EOAS.re', 'parameters.re')
model.connect('EOAS.rho', 'parameters.rho')
model.connect('EOAS.CT', 'parameters.CT')
model.connect('EOAS.R', 'parameters.R')
model.connect('EOAS.W0', 'parameters.W0')
model.connect('EOAS.speed_of_sound', 'parameters.speed_of_sound')
model.connect('EOAS.load_factor', 'parameters.load_factor')
model.connect('EOAS.empty_cg', 'parameters.empty_cg')

# --- Adding design variables ---
parametric = geometry_definition()
for iParametric in range(0, np.size(parametric)):
    name = "helix.geodef_parametric_{:d}_".format(iParametric)
    
    N_chord = np.size(parametric.parametric_def[iParametric].sec)
    N_span = np.size(parametric.parametric_def[iParametric].span)
    span = parametric.parametric_def[iParametric].span
    sweep = parametric.parametric_def[iParametric].sweep
    dihed = parametric.parametric_def[iParametric].dihed
    chord = parametric.parametric_def[iParametric].chord
    twist = parametric.parametric_def[iParametric].twist
    alpha_0 = parametric.parametric_def[iParametric].alpha_0
    alpha_L0 = parametric.parametric_def[iParametric].alpha_L0
    Cl_alpha = parametric.parametric_def[iParametric].Cl_alpha
    M = parametric.parametric_def[iParametric].M

    model.add_design_var(name + "span", shape=N_span, val=span)
    model.add_design_var(name + "sweep", shape=N_span, val=sweep)
    model.add_design_var(name + "dihed", shape=N_span, val=dihed)

    # Section Variables (hard-coding shape for now)
    model.add_design_var(name + "chord", shape=N_chord, val=chord)
    model.add_design_var(name + "twist", shape=N_chord, val=twist)
    model.add_design_var(name + "alpha_0", shape=N_chord, val=alpha_0)
    model.add_design_var(name + "alpha_L0", shape=N_chord, val=alpha_L0)
    model.add_design_var(name + "Cl_alpha", shape=N_chord, val=Cl_alpha)
    model.add_design_var(name + "M", shape=N_chord, val=M)

model.add_design_var('rethorst.span', lower=10, upper=30)
model.add_design_var('rethorst.jet_loc', lower=0.5, upper=4.)
model.add_design_var('EOAS.twist', lower=0.5, upper=4.)
model.add_design_var('EOAS.chord', lower=0.5, upper=4.)

model.add_constraint('constraint1', upper=0.)

model.add_objective('Wfuel', scaler=1.0)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
# prob.driver.options['maxiter'] = 100
prob.driver.options['tol'] = 1e-8

prob.setup()

prob.set_val('master.span', 10.)
prob.set_val('master.chord', 10.)
prob.set_val('master.jet_loc', 10.)
prob.set_val('master.jet_radius', 10.)
prob.set_val('master.Vinf', 10.)
prob.set_val('master.Vjet', 10.)
prob.set_val('master.aoa', 10.)

prob.run_model()