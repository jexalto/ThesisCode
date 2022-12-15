import numpy as np
import openmdao.api as om

class obj_function(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        # self.add_input('fuelburn', val=1.)
        # self.add_input('lift', val=1., units='N')
        # self.add_input('drag', val=1., units='N')
        self.add_input('drag', shape_by_conn=True)
        self.add_input('velocity', shape_by_conn=True)

        self.add_output('objective', val=1.)

        # self.declare_partials('objective', 'fuelburn')
        # self.declare_partials('objective', 'lift')
        # self.declare_partials('objective', 'drag')
        self.declare_partials('objective', 'drag')
        self.declare_partials('objective', 'velocity', val=0.)

    def compute(self, inputs, outputs):
        drag = inputs['drag']
        velocity = inputs['velocity']

        print('Velocity: ', velocity)
        print('Drag: ', drag)
        outputs['objective'] = drag

    def compute_partials(self, inputs, partials):

        # partials['objective', 'lift'] = 0.6
        partials['objective', 'drag'] = 1.