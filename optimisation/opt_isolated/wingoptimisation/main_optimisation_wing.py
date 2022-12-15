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

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=2, panels_span_VLM=40, span_0=0.748*2))

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
        
        self.add_design_var('parameters.chord', lower=0.01, upper=6.0, scaler=1.)
        self.add_design_var('parameters.twist', lower=-5, upper=5, scaler=1.)
        
        # self.add_objective("EOAS.AS_point_0.wing_perf.D", scaler=1.)
        self.add_constraint("EOAS.wing.structural_mass", lower=0.)
        self.add_objective("obj_function.objective", scaler=10)

        self.add_constraint('EOAS.AS_point_0.wing_perf.L', equals=54, scaler=1/54) # works on 16, 1

prob = om.Problem()
model = prob.model

prob.model = master()

prob.setup(mode='rev')

prob.set_val('EOAS.correction', val=0)
prob.set_val('EOAS.velocity_distr', val=40.)

prob.driver = om.pyOptSparseDriver()
prob.driver.options['optimizer'] = 'SNOPT'
prob.driver.opt_settings = {
    "Major feasibility tolerance": 1.0e-5,
    "Major optimality tolerance": 1.0e-8,
    "Minor feasibility tolerance": 1.0e-5,
    "Verify level": -1,
    "Function precision": 1.0e-6,
    # "Major iterations limit": 50,
    "Nonderivative linesearch": None,
    "Print file": "00_results/snopt_output/opt_SNOPT_print.txt",
    "Summary file": "00_results/snopt_output/opt_SNOPT_summary.txt",
}

prob.driver.options['debug_print'] = ['desvars', 'objs']

# prob.run_model()
# prob.model.approx_totals()
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
y = prob['EOAS.wing.mesh'][0, :, 1]

wing_discr=5
# ========== Reset ==========
prob.set_val('parameters.jet_loc', val=[-0.3, 0.3])
prob.set_val('parameters.chord', val=np.ones((wing_discr))*0.15)
prob.set_val('parameters.twist', val=np.zeros(wing_discr))
prob.run_model()

chord_orig = prob.get_val('parameters.chord')
twist_orig = np.zeros(wing_discr)
chord_b = prob.get_val('parameters.chord')
cl_dist_b = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
CD = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.CD'))

chord_len = np.linspace(0, len(chord), len(chord))
cl_dist_len = np.linspace(0, len(chord), len(cl_opt))
chord_wing = interpolate.interp1d(chord_len, chord)
cl_c_dist = cl_opt*chord_wing(cl_dist_len)

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
plt.savefig('00_results/figures/wing_results/cl_distr_chordproploc.png')

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
plt.savefig('00_results/figures/wing_results/clc_distr_chordproploc.png')

x0 = 0; a = 0.2  # x center, half width                                       
y0 = 0; b = 0.748  # y center, half height                                      
x = np.linspace(-0.748, 0.748, len(cl_c_dist))  # x values of interest
ellipse = a*np.sqrt(1-(x)**2/b**2)

span_0 = np.linspace(-y[0], y[0], len(chord_opt))
span_1 = np.linspace(-y[0], y[0], len(twist_orig))
span_2 = np.linspace(-y[0], y[0], len(cl_c_dist))

span_chord = [span_0, span_1, x]

data = [{"Chord\n(m)": chord_orig,
        "Twist\n(deg)": twist_orig,
        "Lift Distribution\n"r"($C_L\cdot c$)": cl_c_dist_},
        {"Chord\n(m)": chord_opt,
        "Twist\n(deg)": twist_opt,
        "Lift Distribution\n"r"($C_L\cdot c$)": cl_c_dist}]

niceplots.setRCParams()

f, axarr = niceplots.stacked_plots(
    "Span (m)",
    span_chord,
    data,
    figsize=(16, 10),
    line_scaler=0.5,
    filename='00_results/figures/wing_results/wing_stacked_plots.png',
    dpi=600,
    legend=['Original', 'Optimized'],
    multiple_x=True,
)
