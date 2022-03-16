import numpy as np
import openmdao.api as om

class proploc(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self);
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('jet_radius', val=1.0, units='m')
        self.add_input('span', val=1.0, units='m')

        self.add_output('constraint1', val=0., units='-')

    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        span = inputs['span']
        jet_loc = inputs['jet_loc']
        jet_radius = inputs['jet_radius']

        constraint = ((abs(jet_loc) + jet_radius))/(span/2*0.95) - 1

        outputs['cosntraint1'] = constraint
