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
from obj_function import obj_function
from constraints import constraints

import json

import niceplots

def wing():

    class master(om.Group):

        def setup(self):
            self.add_subsystem('parameters', subsys=parameters())

            simparam_def = simparam_definition()
            references_def = references_definition()
            geometry_def = geometry_definition(J=1.)
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
            
            span_max = 0.748*2
            nx = 21
            ny = 41

            self.add_subsystem('EOAS', subsys=EOAS(panels_chord_VLM=nx-1, panels_span_VLM=ny-1, radii_shape=N_elem_span+1))

            self.add_subsystem('propinflow', subsys=propinflow(ny=ny, nx=nx, propdist_chord=0.2))
            
            # =================================
            # ===== Connecting Subsystems =====
            # =================================
            self.connect('helix.rotorcomp_0_radii',                         'EOAS.wing.geometry.mesh.jet_radius')

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

    prob.setup()

    prob.set_val('EOAS.correction', val=0)
    prob.set_val('EOAS.velocity_distr', val=40)

    alphas = np.linspace(-8.0, 10.0, 10)

    CL_numerical = np.zeros(10)
    CD_numerical = np.zeros(10)
    count = 0
    for alpha in alphas:

        # Set the alpha in the problem and run analysis
        prob['parameters.alpha'] = alpha
        prob.run_model()

        print()
        print("Angle of attack:", prob["parameters.alpha"])
        print("CL:", prob["EOAS.AS_point_0.wing_perf.CL"])
        print("CD:", prob["EOAS.AS_point_0.wing_perf.CD"])

        CL_numerical[count] = prob["EOAS.AS_point_0.wing_perf.CL"][0]
        CD_numerical[count] = prob["EOAS.AS_point_0.wing_perf.CD"][0]
        
        count+=1
    return alphas, CL_numerical, CD_numerical

# ===========================
# === Plotting of results ===
# ===========================