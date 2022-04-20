import numpy as np
import openmdao.api as om

class obj_function(om.ExplicitComponent):
    def initialize(self):
        pass

    def setup(self):
        # self.add_input('fuelburn', val=1.)
        self.add_input('lift', val=1., units='N')
        self.add_input('drag', val=1., units='N')
        self.add_input('fuelburn', val=1.)#, units='N')
        self.add_input('weight_struc', val=1.)

        self.add_output('objective', val=1.)

        # self.declare_partials('objective', 'fuelburn')
        self.declare_partials('objective', 'lift')
        self.declare_partials('objective', 'drag')
        self.declare_partials('objective', 'fuelburn')
        self.declare_partials('objective', 'weight_struc')

    def compute(self, inputs, outputs):
        lift = inputs['lift']
        drag = inputs['drag']
        fuelburn = inputs['fuelburn']
        weight_struc = inputs['weight_struc']*1000
        total_weight = fuelburn+weight_struc
        print(f'fuelburn:\t{fuelburn} kg')
        print(f'weight_struc:\t{weight_struc} kg')
        print(f'total_weight:\t{total_weight} kg')
        outputs['objective'] = total_weight

    def compute_partials(self, inputs, partials):
        lift = inputs['lift']
        # drag = inputs['drag']

        partials['objective', 'lift'] = 0.
        partials['objective', 'drag'] = 0.
        partials['objective', 'fuelburn'] = 1.
        partials['objective', 'weight_struc'] = 1000