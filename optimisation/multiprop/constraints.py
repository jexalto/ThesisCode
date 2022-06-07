import numpy as np
import openmdao.api as om

class constraints(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('Waw', default=2000, desc='fuselage weight')

    def setup(self):
        self.add_input('lift', val=1.0, units='N')
        self.add_input('fuel_weight', val=1.0, units='N')
        self.add_input('struc_weight', val=1.0, units='N')
        self.add_input('thrust', shape_by_conn=True, units='N')
        self.add_input('drag', val=1.0, units='N')
        self.add_input('jet_loc', val=1.0, units='m')
        self.add_input('span', val=1.0, units='m')

        self.add_output('constraint_lift_weight', val=0.)
        self.add_output('constraint_thrust_drag', val=0.)
        self.add_output('constraint_jetloc', val=0.)
        
        self.declare_partials('constraint_lift_weight', 'lift', method='fd')
        self.declare_partials('constraint_lift_weight', 'fuel_weight', method='fd')
        self.declare_partials('constraint_lift_weight', 'struc_weight', method='fd')
        self.declare_partials('constraint_thrust_drag', 'thrust', val=1.)
        self.declare_partials('constraint_thrust_drag', 'drag', method='fd')
        self.declare_partials('constraint_jetloc', 'jet_loc', method='fd')
        self.declare_partials('constraint_jetloc', 'span', method='fd')

    def compute(self, inputs, outputs):
        lift = inputs['lift']
        fuel_weight = inputs['fuel_weight']
        struc_weight = inputs['struc_weight']
        Waw = self.options['Waw']
        thrust = inputs['thrust'][2][0]
        drag = inputs['drag']
        span = inputs['span']
        jet_loc = inputs['jet_loc']

        constraint_thrust_drag = thrust #1.-thrust/drag
        constraint_lift_weight = 1.-lift/(fuel_weight+struc_weight+Waw)
        constraint_jetloc = abs(jet_loc)/span-1.

        outputs['constraint_thrust_drag'] = constraint_thrust_drag
        outputs['constraint_lift_weight'] = constraint_lift_weight
        outputs['constraint_jetloc'] = constraint_jetloc

    # def compute_partials(self, inputs, partials):
    #     pass
