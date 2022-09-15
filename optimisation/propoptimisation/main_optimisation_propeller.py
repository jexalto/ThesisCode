import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile
import scipy as sp
from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix.openmdao.om_helix import HELIX_Group
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from helix_dir.propweight import propweight
from rethorst_dir.rethorst_group import Rethorst
from obj_function import obj_function
from constraints import constraints
from helix_dir.linear_radius import linear_radius
from printsys import printsys

import json

import niceplots

class master(om.Group):

    def setup(self):
        self.add_subsystem('parameters', subsys=parameters())

        self.add_subsystem('linear_radius', subsys=linear_radius())

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
                                                        power_calc=True,
                                                        velocity_distribution_calc=True,
                                                        ),
        )
        
        N_elem_span = 0
        for iParametric in range(0, np.size(geometry_def)):
            N_span = np.size(geometry_def.parametric_def[iParametric].span)
            for iSpan in range(0, N_span):
                N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span

        self.add_subsystem('obj_function', subsys=obj_function())

        self.add_subsystem('constraints', subsys=constraints())

        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radius',                               'linear_radius.radius')
        self.connect('linear_radius.propspan_sectional',                'helix.geodef_parametric_0_span')

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('helix.rotorcomp_0_power',                         'obj_function.power')

        # ==================================
        # ===== Connecting Constraints =====
        # ==================================
        self.connect('helix.rotorcomp_0_thrust',                        'constraints.thrust')
        
    def configure(self):
        geometry_def = geometry_definition()
        parametric = geometry_def

        if True:
            for iParametric in range(0, np.size(parametric)):
                name = "helix.geodef_parametric_{:d}_".format(iParametric)
                # self.add_design_var(name + "span",lower=0.001, upper=.1)
                # self.add_design_var(name + "chord",lower=0.01, upper=0.2)
                self.add_design_var(name + "twist",lower=15, upper=75)
                self.add_design_var(name + "rot_rate",lower=-1700, upper=-900)
                # self.add_design_var(name + "alpha_0",lower=0.1, upper=1)
                # self.add_design_var(name + "alpha_L0",lower=-0.5, upper=0.5)
                # self.add_design_var(name + "Cl_alpha",lower=np.pi, upper=2.5*np.pi)
        
        self.add_design_var('parameters.radius', lower=0.08, upper=0.2)

        self.add_objective("obj_function.objective", scaler=0.1)

        self.add_constraint("constraints.constraint_thrust", equals=28)

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

prob.setup(mode='fwd')

prob.set_val(name + "span", val=span)
prob.set_val(name + "chord",val=chord)
prob.set_val(name + "twist",val=twist)
prob.set_val(name + "alpha_0",val=alpha_0)
prob.set_val(name + "alpha_L0",val=alpha_L0)
prob.set_val(name + "Cl_alpha",val=Cl_alpha)
# prob.set_val('prop_weight.power', val=10000000)
# prob.set_val('EOAS.correction', val=0)
# prob.set_val('EOAS.velocity_distr', val=40.)

prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings = {
    "Major feasibility tolerance": 1.0e-5,
    "Major optimality tolerance": 1.0e-10,
    "Minor feasibility tolerance": 1.0e-10,
    "Verify level": -1,
    "Function precision": 1.0e-10,
    # "Major iterations limit": 50,
    "Nonderivative linesearch": None,
    "Print file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/propoptimisation/00_results/snopt_output/opt_SNOPT_print.txt",
    "Summary file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/propoptimisation/00_results/snopt_output/opt_SNOPT_summary.txt",
}

span_orig = prob.get_val("helix.geodef_parametric_0_span")
chord_orig  = prob.get_val("helix.geodef_parametric_0_chord")
twist_orig  = prob.get_val("helix.geodef_parametric_0_twist")
prob.run_model()
vel_distr_orig = np.copy(prob.get_val("helix.rotorcomp_0_velocity_distribution"))

# prob.run_model()
# prob.model.approx_totals()
prob.run_driver()
# prob.check_partials(compact_print=True, show_only_incorrect=False, includes=['*helix*'], form='central', step=1e-8) # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
# prob.check_totals(compact_print=True,  form='central')

# ===========================
# === Printing of results ===
# ===========================
print('Prop Radius: ', np.round(prob.get_val('parameters.radius'), 2))
print('Rotation Rate: ', np.round(prob.get_val('helix.geodef_parametric_0_rot_rate'), 2))
print('Power: ', np.round(prob.get_val("helix.rotorcomp_0_power"), 4))
print("Span: ", np.round(prob.get_val("helix.geodef_parametric_0_span"), 4))
print("Chord: ", np.round(prob.get_val("helix.geodef_parametric_0_chord"), 4))
print("Twist: ", np.round(prob.get_val("helix.geodef_parametric_0_twist"), 4))
print("Velocity distribution: ", prob.get_val("helix.rotorcomp_0_velocity_distribution"))

# ===========================
# === Plotting of results ===
# ===========================
span = prob.get_val("helix.geodef_parametric_0_span")
chord = prob.get_val("helix.geodef_parametric_0_chord")
twist = prob.get_val("helix.geodef_parametric_0_twist")
vel_distr = prob.get_val("helix.rotorcomp_0_velocity_distribution")

span_           = np.zeros(len(span)+1)
span_orig_      = np.zeros(len(span)+1)

for iSpan in range(len(span)):
    span_[iSpan+1] = span_[iSpan]+span[iSpan]
    span_orig_[iSpan+1] = span_orig_[iSpan]+span_orig[iSpan]

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(span_, chord, label='Optimised')
ax.plot(span_orig_, chord_orig, label='Original')
ax.set_xlabel(r'Spanwise location $y$')
ax.set_ylabel(r'$m$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/propoptimisation/00_results/figures/prop_results/prop_chordopt.png')
plt.show()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(span_, twist, label='Optimised')
ax.plot(span_orig_, twist_orig, label='Original')
ax.set_xlabel(r'Spanwise location $y$')
ax.set_ylabel(r'Twist $\deg$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/propoptimisation/00_results/figures/prop_results/prop_twistopt.png')
plt.show()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(np.linspace(0, 1, len(vel_distr)), vel_distr, label='Optimised')
ax.plot(np.linspace(0, 1, len(vel_distr)), vel_distr_orig, label='Original')
ax.set_xlabel(r'Spanwise Location $y$')
ax.set_ylabel(r'Exhaust Velocity $m/s$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/propoptimisation/00_results/figures/prop_results/prop_veldistr.png')
plt.show()