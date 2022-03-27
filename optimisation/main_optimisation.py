import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix.openmdao.om_helix import HELIX_Group
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from rethorst_dir.rethorst_group import Rethorst
from constraints.cons_proploc import proploc
# from obj_function import obj_func

class master(om.Group):

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
        
        N_elem_span = 0
        for iParametric in range(0, np.size(geometry_def)):
            N_span = np.size(geometry_def.parametric_def[iParametric].span)
            for iSpan in range(0, N_span):
                N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span
        
        self.add_subsystem('rethorst', subsys=Rethorst(vel_distr_shape=N_elem_span, panels_span_VLM=200))
        
        self.add_subsystem('EOAS', subsys=EOAS(panels_span_VLM=200))

        # self.add_subsystem('obj_fun', subsys=obj_func())
        
        # self.add_subsystem('cons_proploc', subsys=proploc())

        self.connect('helix.rotorcomp_0_velocity_distribution', 'rethorst.velocity_vector')
        self.connect('helix.rotorcomp_0_radii', 'rethorst.radii')

        self.connect('rethorst.correction_matrix', 'EOAS.correction')
        self.connect('rethorst.wing_veldistr', 'EOAS.velocity_distr')

        # --- Connecting parameters ----
        self.connect('parameters.vinf', 'EOAS.v')
        self.connect('parameters.alpha', 'EOAS.alpha')
        self.connect('parameters.Mach_number', 'EOAS.Mach_number')
        self.connect('parameters.re', 'EOAS.re')
        self.connect('parameters.rho', 'EOAS.rho')
        self.connect('parameters.CT', 'EOAS.CT')
        self.connect('parameters.R', 'EOAS.R')
        self.connect('parameters.W0', 'EOAS.W0')
        self.connect('parameters.speed_of_sound', 'EOAS.speed_of_sound')
        self.connect('parameters.load_factor', 'EOAS.load_factor')
        self.connect('parameters.empty_cg', 'EOAS.empty_cg')

        self.connect('parameters.vinf', 'rethorst.Vinf')
        
        self.add_objective('EOAS.AS_point_0.fuelburn')

prob = om.Problem()
model = prob.model

prob.model = master()

# --- Adding design variables ---
parametric = geometry_definition()

if True:
    for iParametric in range(0, np.size(parametric)):
        name = "helix.geodef_parametric_{:d}_".format(iParametric)
        
        N_chord = np.size(parametric.parametric_def[iParametric].sec)
        N_span = np.size(parametric.parametric_def[iParametric].span)
        
        chord = np.zeros(N_chord)
        twist = np.zeros(N_chord)
        alpha_0 = np.zeros(N_chord)
        alpha_L0 = np.zeros(N_chord)
        Cl_alpha = np.zeros(N_chord)
        M = np.zeros(N_chord)
        for iChord in range(0, N_chord):
            chord[iChord] = parametric.parametric_def[iParametric].sec[iChord].chord
            twist[iChord] = parametric.parametric_def[iParametric].sec[iChord].twist
            alpha_0[iChord] = parametric.parametric_def[iParametric].sec[iChord].alpha_0
            alpha_L0[iChord] = parametric.parametric_def[iParametric].sec[iChord].alpha_L0
            Cl_alpha[iChord] = parametric.parametric_def[iParametric].sec[iChord].Cl_alpha
            M[iChord] = parametric.parametric_def[iParametric].sec[iChord].M
        
        span = np.zeros(N_span)
        sweep = np.zeros(N_span)
        dihed = np.zeros(N_span)
        for iSpan in range(0, N_span):
            span[iSpan] = parametric.parametric_def[iParametric].span[iSpan].span
            sweep[iSpan] = parametric.parametric_def[iParametric].span[iSpan].sweep
            dihed[iSpan] = parametric.parametric_def[iParametric].span[iSpan].dihed

        # Section Variables (hard-coding shape for now)
        # model.add_design_var(name + "span",lower=0.1, upper=1)
        # model.add_design_var(name + "chord",lower=0.01, upper=0.1)
        # model.add_design_var(name + "twist",lower=-50, upper=50)
        # model.add_design_var(name + "alpha_0",lower=0.1, upper=1)
        # model.add_design_var(name + "alpha_L0",lower=-0.5, upper=0.5)
        # model.add_design_var(name + "Cl_alpha",lower=np.pi, upper=2.5*np.pi)

model.add_design_var('rethorst.span', lower=10, upper=30) # couple this to EOAS.span
model.add_design_var('rethorst.jet_loc', lower=0.5, upper=4.)
# model.add_design_var('EOAS.wing.geometry.span', lower=10, upper=30)
# model.add_design_var('EOAS.wing.geometry.mesh.jet_loc', lower=-5., upper=5.)
# model.add_design_var('EOAS.wing.geometry.twist_cp', lower=0.5, upper=4.)
# model.add_design_var('EOAS.wing.geometry.chord', lower=0.5, upper=4.)

prob.setup()

prob.set_val(name + "span", val=span)
prob.set_val(name + "chord",val=chord)
prob.set_val(name + "twist",val=twist)
prob.set_val(name + "alpha_0",val=alpha_0)
prob.set_val(name + "alpha_L0",val=alpha_L0)
prob.set_val(name + "Cl_alpha",val=Cl_alpha)
prob.set_val('rethorst.span',val=10.)
prob.set_val('rethorst.jet_loc',val=0.)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
# prob.driver.options['maxiter'] = 100
prob.driver.options['tol'] = 1e-9

prob.run_driver()
# prob.run_model()

# --- Plotting of results ---
print("The fuel burn value is", prob["EOAS.AS_point_0.fuelburn"][0], "[kg]")

cl_dist = prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl')
y = prob['EOAS.wing.mesh'][0, :, 1]
y_ = np.zeros((len(y)-1))

for index in range(len(y)-1):
    y_[index] = (y[index+1]+y[index])/2

plt.plot(y_, cl_dist, label='Cl distribution')
plt.grid()
plt.legend()
plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/figures/cl_distr.png')