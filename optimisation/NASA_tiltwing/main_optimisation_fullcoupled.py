import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile
from scipy import interpolate

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix.openmdao.om_helix import HELIX_Group
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from helix_dir.propweight import propweight
from rethorst_dir.rethorst_group import Rethorst
from EOAS.inducedAoA import propinflow
from helix_dir.linear_radius import linear_radius
from helix_dir.helixcoupling import helixcoupler

from obj_function import obj_function
from constraints import constraints

import json

import niceplots

class master(om.Group):

    def setup(self):
        nx = 2
        ny = 201
        self.span = 13.49
        span_max = 13.5
        self.nr_blades = 5
        self.nr_props = 2

        simparam_def = simparam_definition()
        references_def = references_definition()
        geometry_def = geometry_definition(self.nr_blades)
        N_elem_span = 0
        for iParametric in range(0, np.size(geometry_def)):
            N_span = np.size(geometry_def.parametric_def[iParametric].span)
            for iSpan in range(0, N_span):
                N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span
        
        self.add_subsystem('parameters', subsys=parameters())
        self.add_subsystem('linear_radius0', subsys=linear_radius())
        self.add_subsystem('linear_radius1', subsys=linear_radius())
        self.add_subsystem('helix0', subsys=HELIX_Group(
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
        self.add_subsystem('helix1', subsys=HELIX_Group(
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
        self.add_subsystem('helix_coupler', subsys=helixcoupler(nr_propellers=self.nr_props, nr_blades=self.nr_blades,vel_distr_shape=N_elem_span))
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))
        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=self.span, radii_shape=N_elem_span))

        # self.add_subsystem('prop_weight', subsys=propweight())
        # self.add_subsystem('propinflow', subsys=propinflow(nr_props=nr_props, ny=ny, nx=nx, propdist_chord=0.1))

        self.add_subsystem('constraints', subsys=constraints())
        self.add_subsystem('obj_function', subsys=obj_function())
        
        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radius0',                              'linear_radius0.radius')
        self.connect('parameters.radius1',                              'linear_radius1.radius')
        self.connect('linear_radius0.propspan_sectional',               'helix0.geodef_parametric_0_span')
        self.connect('linear_radius1.propspan_sectional',               'helix1.geodef_parametric_0_span')
        self.connect('helix0.rotorcomp_0_velocity_distribution',        'helix_coupler.vel_distr_0')
        self.connect('helix1.rotorcomp_0_velocity_distribution',        'helix_coupler.vel_distr_1')
        self.connect('helix0.rotorcomp_0_radii',                        'helix_coupler.radii_0')
        self.connect('helix1.rotorcomp_0_radii',                        'helix_coupler.radii_1')

        self.connect('helix_coupler.vel_distr_tot',                     'rethorst.velocity_vector')
        self.connect('helix_coupler.radii_tot',                         'rethorst.radii')
        # self.connect('helix.rotorcomp_0_power',                         'prop_weight.power')

        self.connect('parameters.point_masses',                         'EOAS.AS_point_0.coupled.wing.point_masses')
        self.connect('helix_coupler.radii_tot',                         'EOAS.wing.geometry.mesh.jet_radius')
        self.connect('rethorst.correction_matrix',                      'EOAS.correction')
        self.connect('rethorst.wing_veldistr',                          'EOAS.velocity_distr')

        self.connect('parameters.twist',                                'EOAS.wing.twist_cp')
        self.connect('parameters.chord',                                'EOAS.wing.geometry.chord_cp')
        self.connect('parameters.jet_loc',                              'EOAS.wing.geometry.mesh.jet_loc')
        self.connect('parameters.jet_loc',                              'EOAS.AS_point_0.coupled.wing.point_mass_locations')
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

        # self.connect('EOAS.AS_point_0.coupled.aero_states.horseshoe_circulations',  'propinflow.circulation')
        # self.connect('EOAS.wing.mesh',                                              'propinflow.mesh')
        # self.connect('parameters.jet_loc',                                          'propinflow.jet_loc')
        # self.connect('propinflow.propinflow0',                                      'helix0.simparamdef_v_inf')
        # self.connect('propinflow.propinflow1',                                      'helix1.simparamdef_v_inf')

        # self.nonlinear_solver = om.NonlinearBlockGS(rtol=1e-10, atol=1e-2)
        # self.linear_solver = om.LinearBlockGS()

        # self.nonlinear_solver.options['iprint'] = 2
        # # self.nonlinear_solver.options['maxiter'] = 10
        # # self.nonlinear_solver.options['solve_subsystems'] = True
        # self.nonlinear_solver.linesearch = om.ArmijoGoldsteinLS()
        # self.nonlinear_solver.linesearch.options['maxiter'] = 10
        # self.nonlinear_solver.linesearch.options['iprint'] = 2

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('helix0.rotorcomp_0_power',                         'obj_function.power0')
        self.connect('helix1.rotorcomp_0_power',                         'obj_function.power1')

        # ==================================
        # ===== Connecting Constraints =====
        # ==================================
        self.connect('EOAS.AS_point_0.L_equals_W',                          'constraints.L_W')
        self.connect('helix0.rotorcomp_0_thrust',                           'constraints.thrust0')
        self.connect('helix1.rotorcomp_0_thrust',                           'constraints.thrust1')
        self.connect('EOAS.AS_point_0.wing_perf.CD',                        'constraints.CD')
        self.connect('parameters.rho',                                      'constraints.rho')
        self.connect('parameters.vinf',                                     'constraints.V')
        self.connect('EOAS.AS_point_0.coupled.wing.S_ref',                  'constraints.surface')
        
    def configure(self):
        geometry_def = geometry_definition(self.nr_blades)
        parametric = geometry_def

        if True:
            for iProp in range(self.nr_props):
                for iParametric in range(0, np.size(parametric)):
                    name = "helix{:d}.geodef_parametric_{:d}_".format(iProp, iParametric)
                    self.add_design_var(name + "rot_rate",          lower=20, upper=200, scaler=1/80)
                    self.add_design_var(name + "twist",             lower=10, upper=70, scaler=1/40)

        self.add_design_var('parameters.radius0',                    lower=0.55, upper=2., scaler=1/0.8936)
        self.add_design_var('parameters.radius1',                    lower=0.55, upper=2., scaler=1/0.8936)

        # self.add_design_var('parameters.jet_loc',                   lower=[-self.span/2, 1.], upper=[-1., self.span/2])
        self.add_design_var('parameters.chord',                     lower=0.2, upper=3.0, scaler=1.)
        self.add_design_var('parameters.twist',                     lower=-3.5, upper=5, scaler=1.)

        self.add_objective("obj_function.objective",                scaler=1/50938.53744861)

        self.add_constraint('constraints.constraint_lift_weight',   upper=0.)
        self.add_constraint('constraints.constraint_thrust_drag',   lower=0.)
        self.add_constraint('EOAS.wing.structural_mass',            lower=0.)
        self.add_constraint('EOAS.AS_point_0.wing_perf.failure',    upper=0.)


prob = om.Problem()
model = prob.model
prob.model = master()
wing_discr = 5
nr_blades = 5
# --- Adding design variables ---
parametric = geometry_definition(nr_blades)

prob.setup(mode='fwd')

nr_props = prob.model.nr_props
if True:
    for iProp in range(nr_props):
        for iParametric in range(0, np.size(parametric)):
            name = "helix{:d}.geodef_parametric_{:d}_".format(iProp, iParametric)
            
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

            prob.set_val(name + "chord",val=chord)
            prob.set_val(name + "twist",val=twist)
            prob.set_val(name + "alpha_0",val=alpha_0)
            prob.set_val(name + "alpha_L0",val=alpha_L0)
            prob.set_val(name + "Cl_alpha",val=Cl_alpha)


prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings = {
    "Major feasibility tolerance": 1.0e-5,
    "Major optimality tolerance": 1.0e-7,
    "Minor feasibility tolerance": 1.0e-5,
    "Verify level": -1,
    "Function precision": 1.0e-6,
    # "Major iterations limit": 50,
    "Nonderivative linesearch": None,
    "Print file": f"00_results/snopt_output/opt_SNOPT_print_{wing_discr}sections.txt",
    "Summary file": f"00_results/snopt_output/opt_SNOPT_summary_{wing_discr}sections.txt",
}

prob.driver.options['debug_print'] = ['desvars', 'objs']

span_orig_prop = prob.get_val("helix0.geodef_parametric_0_span")
chord_orig_prop  = prob.get_val("helix0.geodef_parametric_0_chord")
twist_orig_prop  = prob.get_val("helix0.geodef_parametric_0_twist")
# 
prob.run_model()
# prob.model.approx_totals()
# prob.model.list_inputs(includes=['*helix0.geodef_parametric_0_span*', '*helix1.geodef_parametric_0_span*'])
# prob.run_driver()
# prob.check_partials(compact_print=True, show_only_incorrect=True, includes=['rethorst'], form='central', step=1e-8) # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
prob.check_totals(compact_print=True,  form='central')

# ===========================
# === Printing of results ===
# ===========================
print('chord: \t\t', prob.get_val('parameters.chord'))
print('jet loc: \t\t', prob.get_val('parameters.jet_loc'))
print('twist: \t\t', prob.get_val('parameters.twist'))
print('L/D:\t\t', prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.L')/prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.D'))
print('L: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.L'))
print('D: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.D'))
print('CD: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.CD'))
print("The fuel burn value:\t", prob["EOAS.AS_point_0.fuelburn"][0], "[kg]")
print("Structural mass:\t", prob.get_val('EOAS.wing.structural_mass'), " kg")
print()
print('Power: ', prob.get_val("obj_function.objective"))

cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
y = prob['EOAS.wing.mesh'][0, :, 1]
y_ = np.zeros((len(y)-1))

for index in range(len(y)-1):
    y_[index] = (y[index+1]+y[index])/2
_, ax = plt.subplots(figsize=(10, 7))
ax.plot(y_, cl_opt, label='Optimised')
ax.set_xlabel(r'Spanwise location $y$')
ax.set_ylabel(r'$C_L$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('cl_distr_chordproploc.png')

if False:
    # ===========================
    # === Plotting of results ===
    # ===========================
    chord = np.copy(prob.get_val('parameters.chord'))
    CDopt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.CD'))
    chord_opt = np.copy(prob.get_val('parameters.chord'))
    twist_opt = np.copy(prob.get_val('parameters.twist'))
    cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))

    # ========== Reset ==========
    prob.set_val('parameters.jet_loc', val=[-0.3, 0.3])
    prob.set_val('parameters.chord', val=np.ones((5))*0.15)
    prob.set_val('parameters.twist', val=np.zeros(5))

    chord_orig = prob.get_val('parameters.chord')
    twist_orig = prob.get_val('parameters.twist')
    chord_b = prob.get_val('parameters.chord')
    cl_dist_b = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
    CD = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.CD'))

    chord_len = np.linspace(0, len(chord), len(chord))
    cl_dist_len = np.linspace(0, len(chord), len(cl_opt))
    chord_wing = interpolate.interp1d(chord_len, chord)
    cl_c_dist = cl_opt*chord_wing(cl_dist_len)

    y = prob['EOAS.wing.mesh'][0, :, 1]
    span_chord = np.linspace(-y[0]/2, y[0]/2, len(chord_opt))
    y_ = np.zeros((len(y)-1))

    title = f'CD0: {CD}, CDopt: {CDopt}'
    for index in range(len(y)-1):
        y_[index] = (y[index+1]+y[index])/2

    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(y_, cl_dist_b, label='Original')
    ax.plot(y_, cl_opt, label='Optimised')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'$C_L$')
    ax.set_title(title)
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('cl_distr_chordproploc.png')

    chord_len = np.linspace(0, len(chord), len(chord))
    cl_dist_len = np.linspace(0, len(chord), len(cl_opt))

    chord_wing_ = interpolate.interp1d(chord_len, chord_b)
    cl_c_dist_ = cl_dist_b*chord_wing_(cl_dist_len)

    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(y_, cl_c_dist_, label='Original')
    ax.plot(y_, cl_c_dist, label='Optimised')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'$C_L\cdot c$')
    ax.set_title(title)
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/fullycoupled/00_results/figures/clc_distr_chordproploc.png')

    data = [{"Chord (m)": chord_orig,
            "Twist (deg)": twist_orig},
            {"Chord (m)": chord_opt,
            "Twist (deg)": twist_opt}]

    niceplots.setRCParams()

    f, axarr = niceplots.stacked_plots(
        "Span (m)",
        span_chord,
        data,
        figsize=(10, 6),
        line_scaler=0.5,
        filename='/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/fullycoupled/00_results/figures/chord_opt.png',
        dpi=400,
        # labels=[['Original', 'Optimised'], ['Original', 'Optimised']]
    )

    # span = prob.get_val("helix.geodef_parametric_0_span")
    # chord = prob.get_val("helix.geodef_parametric_0_chord")
    # twist = prob.get_val("helix.geodef_parametric_0_twist")

    span_           = np.zeros(len(span)+1)
    span_orig_      = np.zeros(len(span)+1)

    for iSpan in range(len(span)):
        span_[iSpan+1] = span_[iSpan]+span[iSpan]
        span_orig_[iSpan+1] = span_orig_[iSpan]+span_orig_prop[iSpan]

    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(span_, chord, label='Optimised')
    ax.plot(span_orig_, chord_orig_prop, label='Original')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'$m$')
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/prop_results/chord_opt.png')

    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(span_, twist, label='Optimised')
    ax.plot(span_orig_, twist_orig_prop, label='Original')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'Twist $\deg$')
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/prop_results/twist_opt.png')