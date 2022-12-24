import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile
from scipy import interpolate

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
        # =================================
        # ====== Defining parameters ======
        # =================================
        span_max = 0.748*3
        nx = 2
        ny = 141

        simparam_def = simparam_definition()
        references_def = references_definition()
        geometry_def = geometry_definition()
        N_elem_span = 0
        for iParametric in range(0, np.size(geometry_def)):
            N_span = np.size(geometry_def.parametric_def[iParametric].span)
            for iSpan in range(0, N_span):
                N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span

        # =================================
        # ======= Adding Subsystems =======
        # =================================
        self.add_subsystem('parameters', subsys=parameters())

        self.add_subsystem('linear_radius', subsys=linear_radius())
        
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
        
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=0.748*2, radii_shape=N_elem_span+1))

        self.add_subsystem('constraints', subsys=constraints())

        self.add_subsystem('obj_function', subsys=obj_function())
        
        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radius',                               'linear_radius.radius')
        self.connect('linear_radius.propspan_sectional',                'helix.geodef_parametric_0_span')
        self.connect('helix.rotorcomp_0_velocity_distribution',         'rethorst.velocity_vector')
        self.connect('helix.rotorcomp_0_radii',                         'rethorst.radii')
        self.connect('helix.rotorcomp_0_radii',                         'linear_radius.radii')

        self.connect('parameters.point_masses',                         'EOAS.AS_point_0.coupled.wing.point_masses')
        self.connect('helix.rotorcomp_0_radii',                         'EOAS.wing.geometry.mesh.jet_radius')
        self.connect('rethorst.correction_matrix',                      'EOAS.correction')
        self.connect('rethorst.wing_veldistr',                          'EOAS.velocity_distr')

        self.connect('parameters.twist',                                'linear_radius.twist')
        self.connect('linear_radius.twist_list',                        'EOAS.wing.twist_cp')
        self.connect('parameters.chord',                                'linear_radius.chord')
        self.connect('linear_radius.chord_list',                        'EOAS.wing.geometry.chord_cp')
        self.connect('parameters.jet_loc',                              'linear_radius.jet_loc')
        self.connect('linear_radius.jet_loc_list',                      'EOAS.wing.geometry.mesh.jet_loc')
        self.connect('linear_radius.jet_loc_list',                      'EOAS.AS_point_0.coupled.wing.point_mass_locations')
        self.connect('parameters.span',                                 'EOAS.wing.geometry.span')
        self.connect('EOAS.wing.mesh',                                  'linear_radius.mesh')
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

        self.connect('linear_radius.jet_loc_list',                      'rethorst.jet_loc')
        self.connect('parameters.span',                                 'rethorst.span')
        self.connect('parameters.vinf',                                 'rethorst.vinf')

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('helix.rotorcomp_0_power',                         'obj_function.power')

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
                self.add_design_var(name + "rot_rate",              lower=-1700, upper=-800, scaler=1/1000)
                self.add_design_var(name + "twist",lower=10, upper=70, scaler=1/100)

        # self.add_design_var('parameters.radius',                    lower=0.08, upper=0.2, scaler=1.)

        self.add_design_var('parameters.chord',                     lower=0.01, upper=2.0, scaler=10)
        self.add_design_var('parameters.twist',                     lower=-5, upper=10, scaler=100)

        self.add_objective("obj_function.objective",                scaler=1/1000)

        self.add_constraint('constraints.constraint_lift_weight',   equals=0., scaler=10)
        self.add_constraint('constraints.constraint_thrust_drag',   equals=0., scaler=1)
        self.add_constraint('EOAS.wing.structural_mass',            lower=2., scaler=1)
        self.add_constraint('EOAS.AS_point_0.wing_perf.failure',    upper=0., scaler=1000)


prob = om.Problem()
model = prob.model
prob.model = master()

# ===========================
# ======= Adding DVs ========
# ===========================
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

prob.setup(mode='rev')

prob.set_val(name + "span", val=span)
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
    "Nonderivative linesearch": None,
    "Print file": "00_results/snopt_output/opt_SNOPT_print.txt",
    "Summary file": "00_results/snopt_output/opt_SNOPT_summary.txt",
}

span_orig_prop = prob.get_val("helix.geodef_parametric_0_span")
chord_orig_prop  = prob.get_val("helix.geodef_parametric_0_chord")
twist_orig_prop  = prob.get_val("helix.geodef_parametric_0_twist")

prob.driver.options['debug_print'] = ['desvars', 'objs']
# om.n2(prob, outfile='coupled.html')

# prob.run_model()
# prob.model.approx_totals()
prob.run_driver()
# prob.check_partials(compact_print=True, show_only_incorrect=False, includes=['*rethorst*'], form='central', step=1e-8) # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
# prob.check_totals(compact_print=True, form='central')
# prob.cleanup()

# ===========================
# === Printing of results ===
# ===========================
cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
y = prob['EOAS.wing.mesh'][0, :, 1]
with open('00_results/data/meshresults/mesh.txt', 'w') as file:
    np.savetxt(file, y, fmt='%.5f')

y_ = np.zeros((len(y)-1))
for index in range(len(y)-1):
    y_[index] = (y[index+1]+y[index])/2

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

# ===========================
# === Plotting of results ===
# ===========================
if True:
    chordsections = 7
    niceplots.setRCParams()
    
    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(y_, cl_opt, label='Optimised')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'$C_L$')
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('00_results/figures/cl_distr.png')

    chord = np.ones(chordsections)*0.24 # np.copy(prob.get_val('parameters.chord'))
    CDopt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.CD'))
    chord_opt = np.copy(prob.get_val('parameters.chord'))
    twist_opt = np.copy(prob.get_val('parameters.twist'))
    cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))

    # ========== Reset ==========
    prob.set_val('parameters.jet_loc', val=[-0.2])
    prob.set_val('parameters.chord', val=np.ones((4))*0.24)
    prob.set_val('parameters.twist', val=np.zeros(4))

    chord_orig = np.ones(chordsections)*0.24 # prob.get_val('parameters.chord')
    twist_orig = prob.get_val('parameters.twist')
    chord_b = prob.get_val('parameters.chord')
    chord_b = chord_b.append(chord_b[1])
    chord_b = chord_b.append(chord_b[0])
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
    plt.savefig('00_results/figures/cl_distr_chordproploc.png')

    chord_len = np.linspace(0, 5, 5)
    cl_dist_len = np.linspace(0, 5, len(cl_opt))

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
    plt.savefig('00_results/figures/clc_distr_chordproploc.png')

    data = [{"Chord (m)": chord_orig,
            "Twist (deg)": twist_orig},
            {"Chord (m)": chord_opt,
            "Twist (deg)": twist_opt}]

    # f, axarr = niceplots.stacked_plots(
    #     "Span (m)",
    #     span_chord,
    #     data,
    #     figsize=(10, 6),
    #     line_scaler=0.5,
    #     filename='00_results/figures/chord_opt.png',
    #     dpi=400
    # )

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
    ax.plot(span_, chord_orig_prop, label='Original')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'$m$')
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('00_results/figures/prop_results/chord_opt.png')

    _, ax = plt.subplots(figsize=(10, 7))
    ax.plot(span_, twist, label='Optimised')
    ax.plot(span_, twist_orig_prop, label='Original')
    ax.set_xlabel(r'Spanwise location $y$')
    ax.set_ylabel(r'Twist $\deg$')
    ax.legend()
    ax.grid()
    niceplots.adjust_spines(ax, outward=True)
    plt.savefig('00_results/figures/prop_results/twist_opt.png')