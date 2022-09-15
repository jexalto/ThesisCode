import numpy as np 
import openmdao.api as om
from openmdao.test_suite.components.sellar import SellarDis1withDerivatives, SellarDis2withDerivatives

prob = om.Problem()
model = prob.model

model.add_subsystem('d1', SellarDis1withDerivatives(), promotes=['x', 'z', 'y1', 'y2'])
model.add_subsystem('d2', SellarDis2withDerivatives(), promotes=['z', 'y1', 'y2'])

model.add_subsystem('obj_cmp', om.ExecComp('obj = x**2 + z[1] + y1 + exp(-y2)',
                                           z=np.array([0.0, 0.0]), x=0.0),
                    promotes=['obj', 'x', 'z', 'y1', 'y2'])

model.add_subsystem('con_cmp1', om.ExecComp('con1 = 3.16 - y1'), promotes=['con1', 'y1'])
model.add_subsystem('con_cmp2', om.ExecComp('con2 = y2 - 24.0'), promotes=['con2', 'y2'])

model.nonlinear_solver = om.NonlinearBlockGS()

model.linear_solver = om.PETScKrylov()
model.linear_solver.options['iprint'] = 2
model.linear_solver.options['atol'] = 1.0e-20

prob.setup()

prob.set_val('x', 1.)
prob.set_val('z', np.array([5.0, 2.0]))

prob.run_model()

wrt = ['z']
of = ['obj']

J = prob.compute_totals(of=of, wrt=wrt, return_format='flat_dict')