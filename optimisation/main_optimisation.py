import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix.openmdao.om_helix import HELIX_Group
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from helix_dir.propweight import propweight
from rethorst_dir.rethorst_group import Rethorst
from obj_function import obj_function
from constraints import constraints

import niceplots

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
        
        span_max = 0.748*2
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=60))

        self.add_subsystem('prop_weight', subsys=propweight())

        self.add_subsystem('EOAS', subsys=EOAS(panels_span_VLM=60, span_0=0.748, radii_shape=N_elem_span+1))

        self.add_subsystem('obj_function', subsys=obj_function())

        self.add_subsystem('constraints', subsys=constraints())
        
        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('helix.rotorcomp_0_velocity_distribution',         'rethorst.velocity_vector')
        self.connect('helix.rotorcomp_0_radii',                         'rethorst.radii')
        # self.connect('helix.power',                                     'prop_weight.power')

        # self.connect('prop_weight.prop_weight',                         'EOAS.AS_point_0.coupled.wing.point_masses')
        self.connect('helix.rotorcomp_0_radii',                         'EOAS.wing.geometry.mesh.jet_radius')
        self.connect('rethorst.correction_matrix',                      'EOAS.correction')
        self.connect('rethorst.wing_veldistr',                          'EOAS.velocity_distr')

        self.connect('parameters.twist',                                'EOAS.wing.twist_cp')
        self.connect('parameters.chord',                                'EOAS.wing.geometry.chord_cp')
        self.connect('parameters.jet_loc',                              'EOAS.wing.geometry.mesh.jet_loc')
        # self.connect('parameters.jet_loc',                              'EOAS.AS_point_0.coupled.wing.point_mass_locations')
        self.connect('parameters.span',                                 'EOAS.wing.geometry.span')
        self.connect('parameters.vinf',                                 'EOAS.v')
        self.connect('parameters.alpha',                                'EOAS.alpha')
        self.connect('parameters.Mach_number',                          'EOAS.Mach_number')
        self.connect('parameters.re',                                   'EOAS.re')
        self.connect('parameters.rho',                                  'EOAS.rho')
        self.connect('parameters.CT',                                   'EOAS.CT')
        self.connect('parameters.R',                                    'EOAS.R')
        self.connect('parameters.W0',                                   'EOAS.W0')
        self.connect('parameters.speed_of_sound',                       'EOAS.speed_of_sound')
        self.connect('parameters.load_factor',                          'EOAS.load_factor')
        self.connect('parameters.empty_cg',                             'EOAS.empty_cg')

        self.connect('parameters.jet_loc',                              'rethorst.jet_loc')
        self.connect('parameters.span',                                 'rethorst.span')
        self.connect('parameters.vinf',                                 'rethorst.vinf')

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('EOAS.AS_point_0.wing_perf.L',                     'obj_function.lift')
        self.connect('EOAS.AS_point_0.wing_perf.D',                     'obj_function.drag')

        # ==================================
        # ===== Connecting Constraints =====
        # ==================================
        # self.connect('EOAS.AS_point_0.wing_perf.L',                     'constraints.lift')
        # self.connect('EOAS.AS_point_0.total_perf.fuelburn',             'constraints.fuel_weight')
        # self.connect('EOAS.AS_point_0.total_perf.wing_structural_mass', 'constraints.struc_weight')
        # self.connect('helix.rotorcomp_0_thrust',                        'constraints.thrust')
        # self.connect('EOAS.AS_point_0.wing_perf.L',                     'constraints.drag')
        self.connect('parameters.jet_loc',                              'constraints.jet_loc')
        self.connect('parameters.span',                                 'constraints.span')
        
        parametric = geometry_def

        if True:
            for iParametric in range(0, np.size(parametric)):
                name = "helix.geodef_parametric_{:d}_".format(iParametric)
                # self.add_design_var(name + "span",lower=0.001, upper=.05)
                # self.add_design_var(name + "chord",lower=0.01, upper=0.1)
                # self.add_design_var(name + "twist",lower=-50, upper=50)
                # self.add_design_var(name + "alpha_0",lower=0.1, upper=1)
                # self.add_design_var(name + "alpha_L0",lower=-0.5, upper=0.5)
                # self.add_design_var(name + "Cl_alpha",lower=np.pi, upper=2.5*np.pi)
        
        self.add_design_var('parameters.span', lower=1.*0.5, upper=span_max, scaler=1/0.748)
        self.add_design_var('parameters.jet_loc', lower=-0.5*0.748, upper=0.5*0.748, scaler=1/0.374)
        self.add_design_var('parameters.chord', lower=0.2*0.48, upper=4.*0.48, scaler=1/0.48)
        # self.add_design_var('parameters.twist', lower=-3, upper=3, scaler=1.)
        
        self.add_objective("obj_function.objective", scaler=1./0.02931652)
        
        self.add_constraint('EOAS.AS_point_0.wing_perf.L', equals=4.05033303)
        # self.add_constraint('constraints.constraint_thrust_drag', equals=0.)
        # self.add_constraint('constraints.constraint_lift_weight', equals=0.)
        self.add_constraint('constraints.constraint_jetloc', upper=0.)


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
print(span)
prob.set_val(name + "span", val=span)
prob.set_val(name + "chord",val=chord)
prob.set_val(name + "twist",val=twist)
prob.set_val(name + "alpha_0",val=alpha_0)
prob.set_val(name + "alpha_L0",val=alpha_L0)
prob.set_val(name + "Cl_alpha",val=Cl_alpha)
prob.set_val('prop_weight.power', val=10000000)

prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings = {
    "Major feasibility tolerance": 1.0e-5,
    "Major optimality tolerance": 1.0e-4,
    "Minor feasibility tolerance": 1.0e-5,
    "Verify level": -1,
    "Function precision": 1.0e-8,
    "Major iterations limit": 50,
    "Nonderivative linesearch": None,
    "Print file": "/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/snopt_output/opt_SNOPT_print.txt",
    "Summary file": "/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/snopt_output/opt_SNOPT_summary.txt",
}

prob.run_driver()
# prob.run_model()

# prob.check_partials(compact_print=True, show_only_incorrect=False, includes=['*EOAS*'], form='central') # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
# prob.check_totals(compact_print=True,  form='central')

# ===========================
# === Printing of results ===
# ===========================

print('wing span: \t\t', prob.get_val('rethorst.span'))
print('propeller location:\t', prob.get_val('rethorst.jet_loc'))
print('chord: \t\t', prob.get_val('parameters.chord'))
print('twist: \t\t', prob.get_val('parameters.twist'))
print('L/D:\t\t', prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.L')/prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.D'))
print('CL: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.CL'))
print('CD: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.CD'))
print('Propeller radius: \t', prob.get_val('helix.geodef_parametric_0_span'))

print("The fuel burn value is", prob["EOAS.AS_point_0.fuelburn"][0], "[kg]")
print("The lift generated is:\t", prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.L'), " N")

# ===========================
# === Plotting of results ===
# ===========================

cl_dist = prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl')
y = prob['EOAS.wing.mesh'][0, :, 1]
y_ = np.zeros((len(y)-1))

for index in range(len(y)-1):
    y_[index] = (y[index+1]+y[index])/2

_, ax = plt.subplots(figsize=(10, 7))

# Plot Helix Data
ax.plot(y_, cl_dist, label='Cl distribution')

# ax.set_xlim([0.6, 1.05])
ax.set_xlabel(r'Spanwise location $y$')
# ax.set_ylim([0.0, 0.16])
ax.set_ylabel(r'$C_L$')

ax.grid()

# plt.tight_layout()
niceplots.adjust_spines(ax, outward=True)

plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/figures/cl_distr.png')

_, ax = plt.subplots(figsize=(10, 7))

# Plot Helix Data
ax.plot(prob.get_val('helix.rotorcomp_0_radii'), prob.get_val('helix.rotorcomp_0_velocity_distribution'))

# ax.set_xlim([0.6, 1.05])
ax.set_xlabel(r'Propeller radius $R$')
# ax.set_ylim([0.0, 0.16])
ax.set_ylabel(r'Slipstream velocity $V$')

ax.grid()

# plt.tight_layout()
niceplots.adjust_spines(ax, outward=True)

plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/figures/velocitydistribution.png')

_, ax = plt.subplots(figsize=(10, 7))

# Plot Helix Data
veldistr = prob.get_val('rethorst.wing_veldistr')
ax.plot(np.linspace(-0.748/2, 0.748/2, len(veldistr)), veldistr, label='Prop on')
freestream = np.ones(len(veldistr)+2)*40
freestream[0] = 0
freestream[-1] = 0
x = np.concatenate(([-0.748/2], np.linspace(-0.748/2, 0.748/2, len(veldistr)), [0.748/2]))
ax.plot(x, freestream, label='Prop off')

# ax.set_xlim([0.6, 1.05])
ax.set_xlabel(r'Wing span')
ax.set_ylim([0.0, 45])
ax.set_ylabel(r'Velocity distribution $V$')

ax.grid()
ax.legend()

# plt.tight_layout()
niceplots.adjust_spines(ax, outward=True)

plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/figures/wingveldistr.png')
# openmdao iprof_totals <iprof.0>