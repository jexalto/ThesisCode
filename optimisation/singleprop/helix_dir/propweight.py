import numpy as np
import openmdao.api as om

class propweight(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        self.add_input('power', val=1., units='W')
        self.add_output('prop_weight', val=1., units='kg')

        self.declare_partials('prop_weight', 'power', val=11920)

    def compute(self, inputs, outputs):
        power = inputs['power']
        weight = power/(11.92*1000)

        outputs['prop_weight'] = weight