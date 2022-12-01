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
from rethorst_dir.rethorst_group import Rethorst
from EOAS.inducedAoA import propinflow
from helix_dir.helixcoupling import helixcoupler
from helix_dir.linear_radius import linear_radius
from obj_function import obj_function
from constraints import constraints

import json

import niceplots

def wingprop(J):

    class master(om.Group):

        def setup(self):
            self.add_subsystem('parameters', subsys=parameters())

            simparam_def = simparam_definition()
            references_def = references_definition()
            geometry_def = geometry_definition(J)
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
            
            span_max = 5
            nx = 3
            ny = 151
            self.add_subsystem('helix_coupler', subsys=helixcoupler(nr_propellers=1, nr_blades=4,vel_distr_shape=N_elem_span))

            self.add_subsystem('rethorst', subsys=Rethorst(nr_props=1, span_max=span_max, vel_distr_shape=N_elem_span, panels_span_VLM=ny-1, panels_chord_VLM=nx-1))

            self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, span_0=0.748, radii_shape=N_elem_span+1))

            self.add_subsystem('propinflow', subsys=propinflow(ny=ny, nx=nx, propdist_chord=0.2))
            
            self.add_subsystem('linear_radius', subsys=linear_radius())
            
            # =================================
            # ===== Connecting Subsystems =====
            # =================================
            self.connect('parameters.radius',                              'linear_radius.radius')
            self.connect('linear_radius.propspan_sectional',               'helix.geodef_parametric_0_span')
            self.connect('helix.rotorcomp_0_velocity_distribution',        'helix_coupler.vel_distr_0')
            self.connect('helix.rotorcomp_0_radii',                        'helix_coupler.radii_0')
            
            self.connect('helix_coupler.vel_distr_tot',                     'rethorst.velocity_vector')
            self.connect('helix_coupler.radii_tot',                         'rethorst.radii')

            self.connect('parameters.jet_loc',                              'rethorst.jet_loc')
            self.connect('parameters.span',                                 'rethorst.span')
            self.connect('parameters.vinf',                                 'rethorst.vinf')

            self.connect('helix_coupler.radii_tot',                         'EOAS.wing.geometry.mesh.jet_radius')
            self.connect('rethorst.correction_matrix',                      'EOAS.correction')
            self.connect('rethorst.wing_veldistr',                          'EOAS.velocity_distr')

            self.connect('parameters.twist',                                'EOAS.wing.twist_cp')
            self.connect('parameters.chord',                                'EOAS.wing.geometry.chord_cp')
            self.connect('parameters.jet_loc',                              'EOAS.wing.geometry.mesh.jet_loc')
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

            self.connect('EOAS.AS_point_0.coupled.aero_states.horseshoe_circulations',  'propinflow.circulation')
            self.connect('EOAS.wing.mesh',                                              'propinflow.mesh')
            self.connect('parameters.jet_loc',                                          'propinflow.jet_loc')

    prob = om.Problem()
    model = prob.model

    prob.model = master()

    # --- Adding design variables ---
    parametric = geometry_definition(J)

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
     
    steps=6
    alphas = np.linspace(-8.0, 10.0, steps)

    CL_numerical = np.zeros(steps)
    CD_numerical = np.zeros(steps)
    count = 0
    for alpha in alphas:

        # Set the alpha in the problem and run analysis
        prob['parameters.alpha'] = alpha
        print(f'{alpha:_^20}')
        prob.run_model()

        print()
        print("Angle of attack:", prob["parameters.alpha"])
        print("CL:", prob["EOAS.AS_point_0.wing_perf.CL"])
        print("CD:", prob["EOAS.AS_point_0.wing_perf.CD"])

        CL_numerical[count] = prob["EOAS.AS_point_0.wing_perf.CL"][0]
        CD_numerical[count] = prob["EOAS.AS_point_0.wing_perf.CD"][0]
        
        count+=1
    
        cl_opt = np.copy(prob.get_val('EOAS.AS_point_0.wing_perf.aero_funcs.Cl'))
        print(CL_numerical, cl_opt)
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
        plt.savefig(f'00_results/cl_distr_chordproploc{alpha}.png')
    return alphas, CL_numerical, CD_numerical

# ===========================
# === Plotting of results ===
# ===========================