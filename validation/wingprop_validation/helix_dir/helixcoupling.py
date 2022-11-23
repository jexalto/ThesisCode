'''
Author: J. Exalto

Description: This script separates the information about the two propellers that HELIX returns.
'''

import numpy as np
import openmdao.api as om

class helixcoupler(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('nr_propellers', default=2, desc='Nr of propellers')
        self.options.declare('nr_blades', default=5, desc='Nr of propellers')
        self.options.declare('vel_distr_shape', default=20, desc='Vel distr discretisation')

    def setup(self):
        nr_propellers = self.options['nr_propellers']
        nr_blades = self.options['nr_blades']
        vel_distr_shape = self.options['vel_distr_shape']
        for iProp in range(nr_propellers):
            self.add_input('vel_distr_'+str(iProp), shape_by_conn=True, units='m/s')
            self.add_input('radii_'+str(iProp), shape_by_conn=True, units='m')

            cols = np.arange(0, vel_distr_shape, 1)
            rows = np.arange(vel_distr_shape*iProp, vel_distr_shape*(iProp+1), 1)

            self.declare_partials('vel_distr_tot', 'vel_distr_'+str(iProp), rows=rows, cols=cols, val=1)
            self.declare_partials('radii_tot', 'radii_'+str(iProp), rows=rows, cols=cols, val=1)

        self.add_output('vel_distr_tot', val=np.zeros((nr_propellers, vel_distr_shape)), units='m/s')
        self.add_output('radii_tot', val=np.zeros((nr_propellers, vel_distr_shape)), units='m')

    def compute(self, inputs, outputs):
        nr_propellers = self.options['nr_propellers']
        nr_blades = self.options['nr_blades']
        vel_distr_shape = self.options['vel_distr_shape']

        vel_distr_tot   = np.zeros((nr_propellers, vel_distr_shape))
        radii_tot       = np.zeros((nr_propellers, vel_distr_shape))
        
        for iProp in range(nr_propellers):
            print('//////////////////////////////////////')
            print(inputs["vel_distr_"+str(iProp)])
            vel_distr_tot[iProp, :]  = inputs["vel_distr_"+str(iProp)]
            radii_tot[iProp, :]      = inputs["radii_"+str(iProp)]

        outputs['vel_distr_tot'] = vel_distr_tot
        outputs['radii_tot'] = radii_tot
