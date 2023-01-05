import numpy as np
import openmdao.api as om
import matplotlib.pyplot as plt
from openmdao.devtools import iprofile
from scipy import interpolate

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from helix_dir.helix_config import simparam_definition, geometry_definition, references_definition
from rethorst_dir.rethorst_group import Rethorst
from helix_dir.linear_radius import linear_radius

import niceplots

class master(om.Group):

    def setup(self):
        # =================================
        # ====== Defining parameters ======
        # =================================
        span_max = 0.748*2
        nx = 5
        ny = 201

        # simparam_def = simparam_definition()
        # references_def = references_definition()
        # geometry_def = geometry_definition()
        # N_elem_span = 0
        # for iParametric in range(0, np.size(geometry_def)):
        #     N_span = np.size(geometry_def.parametric_def[iParametric].span)
        #     for iSpan in range(0, N_span):
        #         N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span
        N_elem_span =19
        # =================================
        # ======= Adding Subsystems =======
        # =================================
        self.add_subsystem('parameters', subsys=parameters())

        self.add_subsystem('linear_radius', subsys=linear_radius(nr_sections = 41))
        
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=0.748*2, radii_shape=N_elem_span+1))

        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radii',                                'EOAS.wing.geometry.mesh.jet_radius')
        self.connect('parameters.vel_distr',                            'rethorst.velocity_vector')
        self.connect('parameters.radii',                                'rethorst.radii')
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


prob = om.Problem()
model = prob.model
prob.model = master()

# ===========================
# ======= Adding DVs ========
# ===========================
prob.setup(mode='rev')

prob.run_model()
prob.cleanup()

# ===========================
# === Printing of results ===
# ===========================
cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
velocity = np.copy(prob.get_val('rethorst.wing_veldistr'))
velocity_ = np.zeros(len(velocity)+2)
velocity_[1:-1] = velocity
y = prob['EOAS.wing.mesh'][0, :, 1]
with open('00_results/data/meshresults/mesh.txt', 'w') as file:
    np.savetxt(file, y, fmt='%.5f')

y_ = np.zeros((len(y)-1))
for index in range(len(y)-1):
    y_[index] = (y[index+1]+y[index])/2

class master(om.Group):

    def setup(self):
        # =================================
        # ====== Defining parameters ======
        # =================================
        span_max = 0.748*3
        nx = 5
        ny = 201

        # simparam_def = simparam_definition()
        # references_def = references_definition()
        # geometry_def = geometry_definition()
        # N_elem_span = 0
        # for iParametric in range(0, np.size(geometry_def)):
        #     N_span = np.size(geometry_def.parametric_def[iParametric].span)
        #     for iSpan in range(0, N_span):
        #         N_elem_span += geometry_def.parametric_def[iParametric].span[iSpan].N_elem_span
        N_elem_span =19
        # =================================
        # ======= Adding Subsystems =======
        # =================================
        self.add_subsystem('parameters', subsys=parameters())

        self.add_subsystem('linear_radius', subsys=linear_radius(nr_sections = 41))
        
        self.add_subsystem('rethorst', subsys=Rethorst(span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))

        self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=0.748*2, radii_shape=N_elem_span+1))

        # =================================
        # ===== Connecting Subsystems =====
        # =================================
        self.connect('parameters.radii',                                'EOAS.wing.geometry.mesh.jet_radius')
        self.connect('parameters.vel_distr',                            'rethorst.velocity_vector')
        self.connect('parameters.radii',                                'rethorst.radii')
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


prob = om.Problem()
model = prob.model
prob.model = master()

# ===========================
# ======= Adding DVs ========
# ===========================
prob.setup(mode='rev')

prob.set_val('EOAS.correction', val=0.)
prob.run_model()
cl_opt_ = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))

# ===========================
# === Plotting of results ===
# ===========================
chordsections = 7
niceplots.setRCParams()

_, (ax,ax2) = plt.subplots(2, 1, figsize=(16, 12))
ax.plot(y_, cl_opt, label='Lift distribution')
ax.plot(y_, cl_opt_, label='No Correction')
ax.set_ylabel(r'$C_L$')
ax.legend()
# ax.grid()
niceplots.adjust_spines(ax, outward=True)
ax2.plot(np.linspace(y_[0], y_[-1], len(velocity_)), velocity_, label='Velocity distribution')
ax2.set_xlabel(r'Spanwise location $y$')
ax2.set_ylabel(r'Velocity $m/s$')
ax2.legend()
# ax.grid()
niceplots.adjust_spines(ax2, outward=True)
plt.savefig('00_results/figures/cl_distr.png')