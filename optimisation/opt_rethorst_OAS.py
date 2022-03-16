import numpy as np
import openmdao.api as om

from openaerostruct.utils.constants import grav_constant

from parametersinput import parameters
from EOAS.EOAS_group import EOAS
from rethorst.rethorst_system import Rethorst
from constraints.cons_proploc import proploc
from obj_function import obj_func

class master(om.group):

    def setup(self):
        main = self.add_subsystem('main', om.Group(), promotes_inputs=[())])

        main.add_subsystem('parameters', subsys=parameters(),           promotes_outputs=[( '*')])

        main.add_subsystem('rethorst', subsys=Rethorst(panels_VLM=71),  promotes_inputs=[(  'rethorst.span', 'rethorst.jet_loc', 'rethorst.jet_radius',
                                                                                            'rethorst.vinf', 'rethorst.vjet')],
                                                                        promotes_outputs=[( 'rethorst.correction_matrix', 'rethorst.mesh')])
        
        main.add_subsystem('EOAS', subsys=EOAS(),                       promotes_outputs=[( 'EOAS.CL', 'EOAS.CD', 'EOAS.cl', 'EOAS.cd',
                                                                                            'EOAS.Wstruc', 'EOAS.Wfuel')])

        self.add_subsystem('obj_fun', subsys=obj_func,                  promotes_inputs=[(  'Wfuel', 'Wstruc', 'CL', 'CD', 'cl', 'cd')],
                                                                        promotes_outputs=[( 'f')])
        
        self.add_subsystem('cons_proploc', subsys=proploc(),            promotes_input=[(   'jet_loc', 'jet_radius', 'span')],
                                                                        promotes_output=[(  'constraint1')])

        
        self.connect('rethorst.correction_matrix', 'EOAS.correction_matrix')
        self.connect('rethorst.mesh', 'EOAS.mesh')

        self.add_design_var('span', lower=10, upper=30)
        self.add_design_var('chord', lower=0.5, upper=4)
        self.add_design_var('jet_loc', lower=0, upper=10)
        self.add_design_var('jet_radius', lower=0.5, upper=.5)
        self.add_design_var('Vinf', lower=0, upper=28.5)
        self.add_design_var('Vjet', lower=0, upper=28.5)
        self.add_design_var('aoa', lower=0, upper=28.5)

        self.add_constraint('constraint1', upper=0.)

        self.add_objective('Wfuel', scaler=1.0)

prob = om.Problem()
model = prob.model

prob.model = master()

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
# prob.driver.options['maxiter'] = 100
prob.driver.options['tol'] = 1e-8

prob.setup()

prob.set_val('master.span', 10.)
prob.set_val('master.chord', 10.)
prob.set_val('master.jet_loc', 10.)
prob.set_val('master.jet_radius', 10.)
prob.set_val('master.Vinf', 10.)
prob.set_val('master.Vjet', 10.)
prob.set_val('master.aoa', 10.)

prob.run_model()