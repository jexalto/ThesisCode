import numpy as np
import openmdao.api as om

class obj_function(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        # self.add_input('fuelburn', val=1.)
        # self.add_input('lift', val=1., units='N')
        # self.add_input('drag', val=1., units='N')
        self.add_input('power', shape_by_conn=True)

        self.add_output('objective', val=1.)

        # self.declare_partials('objective', 'fuelburn')
        # self.declare_partials('objective', 'lift')
        # self.declare_partials('objective', 'drag')
        self.declare_partials('objective', 'power')

    def compute(self, inputs, outputs):
        power = inputs['power']

        print('Strucweight: ', np.round(power, 6))
        outputs['objective'] = power[0]

    def compute_partials(self, inputs, partials):

        # partials['objective', 'lift'] = 0.6
        partials['objective', 'power'] = 1.