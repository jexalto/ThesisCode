import numpy as np
import openmdao.api as om

class obj_function(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        # self.add_input('fuelburn', val=1.)
        self.add_input('lift', val=1., units='N')
        self.add_input('drag', val=1., units='N')

        self.add_output('objective', val=1.)

        # self.declare_partials('objective', 'fuelburn')
        self.declare_partials('objective', 'lift')
        self.declare_partials('objective', 'drag')

    def compute(self, inputs, outputs):
        lift = inputs['lift']
        drag = inputs['drag']
        
        outputs['objective'] = drag

    def compute_partials(self, inputs, partials):
        lift = inputs['lift']
        # drag = inputs['drag']

        partials['objective', 'lift'] = 0.
        partials['objective', 'drag'] = 1.