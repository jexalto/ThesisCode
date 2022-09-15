import numpy as np
import openmdao.api as om

class constraints(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('Waw', default=2000, desc='fuselage weight')

    def setup(self):
        self.add_input('L_W', val=1.0, units='N')
        self.add_input('thrust', shape_by_conn=True, units='N')
        self.add_input('drag', val=1.0, units='N')

        self.add_output('constraint_lift_weight', val=0.)
        self.add_output('constraint_thrust_drag', val=0.)
        self.add_output('constraint_thrust', val=0.)
        
        self.declare_partials('constraint_thrust', 'thrust', cols=[10], rows=[0], val=1.)
        self.declare_partials('constraint_lift_weight', 'L_W', val=1.)
        self.declare_partials('constraint_lift_weight', 'thrust', val=0)
        self.declare_partials('constraint_lift_weight', 'drag', val=0)
        self.declare_partials('constraint_thrust_drag', 'L_W', val=0)
        self.declare_partials('constraint_thrust_drag', 'thrust', rows=[0], cols=[10])
        self.declare_partials('constraint_thrust_drag', 'drag')

    def compute(self, inputs, outputs):
        LW = inputs['L_W']
        thrust = inputs['thrust'][2][0]
        drag = inputs['drag']*100

        print(f'Lift equals Weight: {LW}')
        print('Thrust equals Drag: ', thrust-drag, 'Thrust: ', thrust, 'Drag: ', drag[0])
        outputs['constraint_lift_weight'] = LW
        outputs['constraint_thrust_drag'] = thrust-drag
        outputs['constraint_thrust'] = thrust
        print('Thrust: ', thrust)

    def compute_partials(self, inputs, partials):
        thrust = inputs['thrust'][2][0]
        drag = inputs['drag']

        partials['constraint_thrust_drag', 'drag'] = -100 #thrust/drag**2
        partials['constraint_thrust_drag', 'thrust'] = 1
