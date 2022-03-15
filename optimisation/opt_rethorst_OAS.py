import numpy as np
import openmdao.api as om

from EOAS_system import EOAS
from rethorst_system import Rethorst

class Circuit(om.group):

    def setup(self):
        self.add_subsystem('rethorst', subsys=Rethorst(panels_VLM=71),  promotes_inputs=[(  'rethorst.span', 'rethorst.jet_loc', 'rethorst.jet_radius',
                                                                                            'rethorst.Vinf', 'rethorst.Vjet')],
                                                                        promotes_outputs=[( 'rethorst.correction_matrix', 'rethorst.mesh')])
        
        self.add_subsystem('EOAS', subsys=EOAS(panels_VLM=71),          promotes_inputs=[(  'EOAS.mesh', 'EOAS.correction_matrix', 'EOAS.chord', 
                                                                                            'EOAS.jet_loc', 'EOAS.jet_radius',
                                                                                            'EOAS.Vinf', 'EOAS.Vjet',
                                                                                            'EOAS.aoa')],
                                                                        promotes_outputs=[( 'EOAS.CL', 'EOAS.CD', 'EOAS.cl', 'EOAS.cd',
                                                                                            'EOAS.Wstruc', 'EOAS.Wfuel')])

        self.connect('transfer.va', 'dv2.v1')

        self.add_design_var('span', lower=10, upper=30)
        self.add_design_var('chord', lower=0.5, upper=4)
        self.add_design_var('jet_loc', lower=0, upper=10)
        self.add_design_var('jet_radius', lower=0.5, upper=.5)
        self.add_design_var('Vinf', lower=0, upper=28.5)
        self.add_design_var('Vjet', lower=0, upper=28.5)
        self.add_design_var('aoa', lower=0, upper=28.5)

        self.add_constraint('dinc', lower=28.5, upper=28.5, scaler=1.0)

        self.add_objective('delta_v', scaler=1.0)

prob = om.Problem()
model = prob.model

model.add_subsystem('circuit', Circuit())

prob.setup()

prob.set_val('circuit.I_in', 0.1)
prob.set_val('circuit.Vg', 0.)

# set some initial guesses
prob.set_val('circuit.n1.V', 10.)
prob.set_val('circuit.n2.V', 1e-3)

prob.run_model()