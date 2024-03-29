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

from obj_function import obj_function
from constraints import constraints

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
        
        span_max = 0.748*5
        nx = 3
        ny = 81
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))

        # self.add_subsystem('prop_weight', subsys=propweight())

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=0.748*2*0.96, radii_shape=N_elem_span+1))

        self.add_subsystem('propinflow', subsys=propinflow(ny=ny, nx=nx, propdist_chord=0.1))

        self.add_subsystem('constraints', subsys=constraints())

        self.add_subsystem('obj_function', subsys=obj_function())
        
        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radius',                               'linear_radius.radius')
        self.connect('linear_radius.propspan_sectional',                'helix.geodef_parametric_0_span')
        self.connect('helix.rotorcomp_0_velocity_distribution',         'rethorst.velocity_vector')
        self.connect('helix.rotorcomp_0_radii',                         'rethorst.radii')
        # self.connect('helix.rotorcomp_0_power',                         'prop_weight.power')

        self.connect('parameters.point_masses',                         'EOAS.AS_point_0.coupled.wing.point_masses')
        self.connect('helix.rotorcomp_0_radii',                         'EOAS.wing.geometry.mesh.jet_radius')
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

        self.connect('parameters.vinf',                                             'propinflow.vinf')
        self.connect('EOAS.AS_point_0.coupled.aero_states.horseshoe_circulations',  'propinflow.circulation')
        self.connect('EOAS.wing.mesh',                                              'propinflow.mesh')
        self.connect('parameters.jet_loc',                                          'propinflow.jet_loc')
        self.connect('propinflow.propinflow',                                       'helix.simparamdef_v_inf')

        self.nonlinear_solver = om.NewtonSolver(solve_subsystems=True, atol=1e-6, rtol=1e-99)
        self.linear_solver = om.ScipyKrylov()
        # self.linear_solver.precon = om.LinearBlockJac()
        # self.linear_solver.precon.options['maxiter'] = 10
        self.linear_solver.options['maxiter'] = 15
        self.nonlinear_solver.options['maxiter'] = 30
        self.nonlinear_solver.linesearch = om.ArmijoGoldsteinLS()
        self.nonlinear_solver.linesearch.options['maxiter'] = 50
        self.nonlinear_solver.linesearch.options['iprint'] = 2

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('helix.rotorcomp_0_power',                         'obj_function.power')
        # self.connect('EOAS.AS_point_0.wing_perf.D',                     'obj_function.power')
        # self.connect('EOAS.wing.structural_mass',                       'obj_function.power')

        # ==================================
        # ===== Connecting Constraints =====
        # ==================================
        self.connect('EOAS.AS_point_0.L_equals_W',                          'constraints.L_W')
        self.connect('helix.rotorcomp_0_thrust',                            'constraints.thrust')
        self.connect('EOAS.AS_point_0.wing_perf.CD',                        'constraints.CD')
        self.connect('EOAS.AS_point_0.coupled.wing.S_ref',                  'constraints.CDv')
        self.connect('parameters.rho',                                      'constraints.rho')
        self.connect('parameters.vinf',                                     'constraints.V')
        self.connect('EOAS.AS_point_0.coupled.wing.S_ref',                  'constraints.surface')
        
    def configure(self):
        geometry_def = geometry_definition()
        parametric = geometry_def

        if True:
            for iParametric in range(0, np.size(parametric)):
                name = "helix.geodef_parametric_{:d}_".format(iParametric)
                self.add_design_var(name + "rot_rate",              lower=-1700, upper=-1200, scaler=1/1000)
                self.add_design_var(name + "twist",                 lower=20, upper=70, scaler=1/10)
                # self.add_design_var(name + "span",lower=0.003, upper=.05)
                # self.add_design_var(name + "chord",lower=0.05, upper=0.1)
                # self.add_design_var(name + "alpha_0",lower=0.1, upper=1)
                # self.add_design_var(name + "alpha_L0",lower=-0.5, upper=0.5)
                # self.add_design_var(name + "Cl_alpha",lower=np.pi, upper=2.5*np.pi)

        self.add_design_var('parameters.radius',                    lower=0.08, upper=0.2, scaler=10)

        # self.add_design_var('parameters.span', lower=1.*0.5, upper=3., scaler=1/0.748)
        # self.add_design_var('parameters.jet_loc',                   lower=[-0.748*0.96/2, 0.1], upper=[-0.1, 0.748*0.96/2])
        self.add_design_var('parameters.chord',                     lower=0.015, upper=2.0, scaler=10.)
        self.add_design_var('parameters.twist',                     lower=-5, upper=5, scaler=1.)

        self.add_objective("obj_function.objective",                scaler=1/1000)

        self.add_constraint('constraints.constraint_lift_weight',   upper=0., scaler=10)
        self.add_constraint('constraints.constraint_thrust_drag',   lower=0., scaler=10)
        # self.add_constraint("constraints.constraint_thrust",        equals=28)
        self.add_constraint('EOAS.wing.structural_mass',            lower=0.)
        self.add_constraint('EOAS.AS_point_0.wing_perf.failure',    upper=0.)


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
prob.set_solver_print(level=2)
# om.n2(prob)
prob.set_val(name + "span", val=span)
prob.set_val(name + "chord",val=chord)
prob.set_val(name + "twist",val=twist)
prob.set_val(name + "alpha_0",val=alpha_0)
prob.set_val(name + "alpha_L0",val=alpha_L0)
prob.set_val(name + "Cl_alpha",val=Cl_alpha)
# prob.set_val('EOAS.correction', val=0)
# prob.set_val('EOAS.velocity_distr', val=40.)
# prob.set_val('prop_weight.power', val=[12000, 12000])

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
    "Print file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/snopt_output/opt_SNOPT_print_full.txt",
    "Summary file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/snopt_output/opt_SNOPT_summary_full.txt",
}

span_orig_prop = prob.get_val("helix.geodef_parametric_0_span")
chord_orig_prop  = prob.get_val("helix.geodef_parametric_0_chord")
twist_orig_prop  = prob.get_val("helix.geodef_parametric_0_twist")

# recorder = om.SqliteRecorder("/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/data/cases.sql")
# prob.driver.add_recorder(recorder)
# prob.driver.recording_options["record_objectives"] = True
# prob.driver.recording_options["record_constraints"] = True
# prob.driver.recording_options["record_desvars"] = True
# prob.driver.recording_options["record_inputs"] = True
# prob.driver.recording_options['includes'] = ['helix.rotorcomp_0_thrust']

prob.run_model()
# prob.model.approx_totals()
# prob.run_driver()
# prob.check_partials(compact_print=True, show_only_incorrect=False, includes=['*propinflow*'], form='central', step=1e-8) # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
prob.check_totals(compact_print=True,  form='central')
prob.cleanup()

# ===========================
# === Printing of results ===
# ===========================
print('chord: \t\t', prob.get_val('parameters.chord'))
print('jet loc: \t\t', prob.get_val('parameters.jet_loc'))
print('radius: \t\t', prob.get_val('parameters.radius'))
print('twist: \t\t', prob.get_val('parameters.twist'))
print('L/D:\t\t', prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.L')/prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.D'))
print('L: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.L'))
print('D: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.D'))
print('CD: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.CD'))
print("The fuel burn value:\t", prob["EOAS.AS_point_0.fuelburn"][0], "[kg]")
print("Structural mass:\t", prob.get_val('EOAS.wing.structural_mass'), " kg")
print()
print('Power: ', prob.get_val("obj_function.objective"))
print("Span: ", prob.get_val("helix.geodef_parametric_0_span"))
print("Chord: ", prob.get_val("helix.geodef_parametric_0_chord"))
print("Velocity distribution: ", prob.get_val("helix.rotorcomp_0_velocity_distribution"))

if True:
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
    plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/cl_distr_chordproploc.png')

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
    plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/clc_distr_chordproploc.png')

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
        filename='/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/chord_opt.png',
        dpi=400,
        # labels=[['Original', 'Optimised'], ['Original', 'Optimised']]
    )

    span = prob.get_val("helix.geodef_parametric_0_span")
    chord = prob.get_val("helix.geodef_parametric_0_chord")
    twist = prob.get_val("helix.geodef_parametric_0_twist")

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