import numpy as np
import openmdao.api as om

class propweight(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        self.add_input('power', shape_by_conn=True, units='W')
        self.add_output('prop_weight', shape=2, units='kg')

        self.declare_partials('prop_weight', 'power', val=1/11920)

    def compute(self, inputs, outputs):
        weight = np.zeros(2)
        power = inputs['power'][0]
        weight[0] = power/(11.92*1000)
        weight[1] = weight[0]
        outputs['prop_weight'] = weight