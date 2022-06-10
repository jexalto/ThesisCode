import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile
from scipy import interpolate

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group_wing import EOAS
from obj_function_wing import obj_function
from constraints import constraints

import niceplots

class master(om.Group):

    def setup(self):
        self.add_subsystem('parameters', subsys=parameters())

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=2+2, panels_span_VLM=100, span_0=0.748))

        self.add_subsystem('obj_function', subsys=obj_function())
        
        # =================================
        # ===== Connecting Subsystems =====
        # =================================

        self.connect('parameters.twist',                                'EOAS.wing.twist_cp')
        self.connect('parameters.chord',                                'EOAS.wing.geometry.chord_cp')
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

        # ================================      
        # ===== Connecting Objective =====      
        # ================================      
        self.connect('EOAS.AS_point_0.wing_perf.D',                     'obj_function.drag')

        # ==================================
        # ===== Connecting Constraints =====
        # ==================================
        # self.connect('EOAS.AS_point_0.wing_perf.L',                     'constraints.lift')
        # self.connect('EOAS.AS_point_0.total_perf.fuelburn',             'constraints.fuel_weight')
        # self.connect('EOAS.AS_point_0.total_perf.wing_structural_mass', 'constraints.struc_weight')
        
    def configure(self):
        
        self.add_design_var('parameters.chord', lower=0.01, upper=2.0, scaler=100.)
        self.add_design_var('parameters.twist', lower=-5, upper=5, scaler=1.)
        
        # self.add_objective("EOAS.AS_point_0.wing_perf.D", scaler=1.)
        # self.add_objective("EOAS.wing.structural_mass")
        self.add_objective("obj_function.objective", scaler=0.1)

        self.add_constraint('EOAS.AS_point_0.wing_perf.L', equals=15.11)

prob = om.Problem()
model = prob.model

prob.model = master()

case_file = '/home/jexalto/code/MDO_lab_env/ThesisCode/optimisation/00_results/data/cases.sql'

prob.setup(mode='fwd')

prob.set_val('EOAS.correction', val=0)
prob.set_val('EOAS.velocity_distr', val=40.)

prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings = {
    "Major feasibility tolerance": 1.0e-5,
    "Major optimality tolerance": 1.0e-5,
    "Minor feasibility tolerance": 1.0e-5,
    "Verify level": -1,
    "Function precision": 1.0e-6,
    # "Major iterations limit": 50,
    "Nonderivative linesearch": None,
    "Print file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/snopt_output/opt_SNOPT_print.txt",
    "Summary file": "/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/snopt_output/opt_SNOPT_summary.txt",
}

# prob.run_model()
prob.model.approx_totals()
prob.run_driver()
# prob.check_partials(compact_print=True, show_only_incorrect=True, includes=['*EOAS.wing.geometry*'], form='central', step=1e-8) # excludes=['*parameters*, *helix*, *EOAS*, *rethorst*']
# prob.check_totals(compact_print=True,  form='central')


# ===========================
# === Printing of results ===
# ===========================
print('twist: \t\t', prob.get_val('parameters.twist'))
print('chord: \t\t', prob.get_val('parameters.chord'))
print('L/D:\t\t', prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.L')/prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.D'))
print('L: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.L'))
print('D: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.D'))
print('CD: \t\t', prob.get_val('EOAS.AS_point_0.wing_perf.CD'))
print("The fuel burn value:\t", prob["EOAS.AS_point_0.fuelburn"][0], "[kg]")
print("Structural mass:\t", prob.get_val('EOAS.wing.structural_mass'), " kg")

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
prob.run_model()

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
# ax.set_title(title)
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/wing_results/cl_distr_chordproploc.png')

chord_len = np.linspace(0, len(chord), len(chord))
cl_dist_len = np.linspace(0, len(chord), len(cl_opt))

chord_wing_ = interpolate.interp1d(chord_len, chord_b)
cl_c_dist_ = cl_dist_b*chord_wing_(cl_dist_len)

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(y_, cl_c_dist_, label='Original')
ax.plot(y_, cl_c_dist, label='Optimised')
ax.set_xlabel(r'Spanwise location $y$')
ax.set_ylabel(r'$C_L\cdot c$')
# ax.set_title(title)
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/wing_results/clc_distr_chordproploc.png')

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
    filename='/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/figures/wing_results/chord_opt.png',
    dpi=400
    # labels=[['Original', 'Optimised'], ['Original', 'Optimised']]
)
