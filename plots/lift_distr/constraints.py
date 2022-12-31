import numpy as np
import openmdao.api as om

class constraints(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('Waw', default=2000, desc='fuselage weight')

    def setup(self):
        self.add_input('L_W', val=1.0, units='N')
        self.add_input('thrust', shape_by_conn=True, units='N')
        self.add_input('CD', val=1.0, units='N')
        self.add_input('CDv', val=0.)
        self.add_input('rho', val=0., units='kg/m**3')
        self.add_input('V', val=0., units='m/s')
        self.add_input('surface', val=0., units='m**2')

        self.add_output('constraint_lift_weight', val=0.)
        self.add_output('constraint_thrust_drag', val=0.)
        self.add_output('constraint_thrust', val=0.)
        
        self.declare_partials('constraint_thrust', 'thrust', cols=[10], rows=[0], val=1.)
        self.declare_partials('constraint_lift_weight', 'L_W', val=1.)
        self.declare_partials('constraint_lift_weight', 'thrust', val=0)
        self.declare_partials('constraint_lift_weight', 'CD', val=0)
        self.declare_partials('constraint_lift_weight', 'CDv', val=0)
        self.declare_partials('constraint_lift_weight', 'rho', val=0)
        self.declare_partials('constraint_lift_weight', 'V', val=0)
        self.declare_partials('constraint_lift_weight', 'surface', val=0)
        self.declare_partials('constraint_thrust_drag', 'L_W', val=0)
        self.declare_partials('constraint_thrust_drag', 'thrust', rows=[0], cols=[10])
        self.declare_partials('constraint_thrust_drag', 'CD')
        self.declare_partials('constraint_thrust_drag', 'CDv')
        self.declare_partials('constraint_thrust_drag', 'rho')
        self.declare_partials('constraint_thrust_drag', 'V')
        self.declare_partials('constraint_thrust_drag', 'surface')

    def compute(self, inputs, outputs):
        LW = inputs['L_W']
        thrust = inputs['thrust'][2][0]*2
        CD = inputs['CD']
        rho = inputs["rho"]
        V = inputs["V"]
        S = inputs["surface"]

        drag = 0.5*(CD)*rho*V**2*S
        
        print()
        print('=====================')
        print('==== Constraints ====')
        print('=====================')
        print(f'Lift equals Weight: {LW}') # This is 1-lift/weight
        print('Thrust equals Drag: ', 1-drag/thrust, 'Thrust: ', thrust, 'Drag: ', drag, '\n')
        
        outputs['constraint_lift_weight'] = LW
        outputs['constraint_thrust_drag'] = 1-drag/thrust
        outputs['constraint_thrust'] = -thrust

    def compute_partials(self, inputs, partials):
        thrust = inputs['thrust'][2][0]
        CD = inputs['CD']
        CDv = inputs["CDv"]
        rho = inputs["rho"]
        V = inputs["V"]
        S = inputs["surface"]

        partials['constraint_thrust_drag', 'thrust'] = 2*(0.5*rho*V**2*S*(CD+CDv))/(2*thrust)**2
        partials['constraint_thrust_drag', 'CD'] = -(0.5*rho*V**2*S)/(2*thrust)
        partials['constraint_thrust_drag', 'CDv'] = -(0.5*rho*V**2*S)/(2*thrust)
        partials['constraint_thrust_drag', 'rho'] = -(0.5*V**2*S*(CD+CDv))/(2*thrust)
        partials['constraint_thrust_drag', 'V'] = -(rho*V*S*(CD+CDv))/(2*thrust)
        partials['constraint_thrust_drag', 'surface'] = -(0.5*rho*V**2*(CD+CDv))/(2*thrust)
