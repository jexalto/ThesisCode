import openmdao.api as om
import numpy as np
import matplotlib.pyplot as plt

import niceplots

cr = om.CaseReader('/home/jexalto99/code/MDO_lab_env/ThesisCode/optimisation/multiprop/00_results/data/cases.sql')

# get the last driver case
driver_cases = cr.list_cases('driver')

cases = cr.get_cases()

iterations = len(cases)

iter_list = np.arange(0, iterations, 1)

param_HELIX_radius          = np.zeros(iterations)
param_HELIX_J               = np.zeros(iterations)
# param_HELIX_twist           = np.zeros(iterations)
# param_OAS_twist             = np.zeros(iterations)
# param_OAS_chord             = np.zeros(iterations)
cons_LW                     = np.zeros(iterations)
cons_TD                     = np.zeros(iterations)
cons_Wstruc                 = np.zeros(iterations)
cons_failure                = np.zeros(iterations)
objective                   = np.zeros(iterations)

for iCase, case in enumerate(cases):
    param_HELIX_radius[iCase]      = case['parameters.radius']
    param_HELIX_J[iCase]           = case['helix.geodef_parametric_0_rot_rate']/(2*np.pi)*60
    # param_HELIX_twist[iCase]       = case['helix.geodef_parametric_0_rot_twist']
    # param_OAS_twist[iCase]         = case['parameters.twist']
    # param_OAS_chord[iCase]         = case['parameters.twist']

    cons_LW[iCase]                 = case['constraints.constraint_lift_weight']
    cons_TD[iCase]                 = case['constraints.constraint_thrust_drag']
    cons_Wstruc[iCase]             = case['EOAS.wing.structural_mass']
    cons_failure[iCase]            = case['EOAS.AS_point_0.wing_perf.failure']

    objective[iCase]               = case['obj_function.objective']

niceplots.setRCParams()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(iter_list, param_HELIX_radius, label='Radius')
ax.set_xlabel(r'Iteration $-$')
ax.set_ylabel(r'Propeller Radius $m$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.show()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(iter_list, param_HELIX_J, label='RPM')
ax.set_xlabel(r'Iteration $-$')
ax.set_ylabel(r'Rotational Speed $RPM$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.show()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(iter_list, cons_LW, label='LW')
ax.plot(iter_list, cons_TD, label='TD')
ax.plot(iter_list, cons_Wstruc, label='Wstruc')
# ax.plot(iter_list, cons_failure, label='Failure')
ax.set_xlabel(r'Iteration $-$')
ax.set_ylabel(r'Constraint Value $-$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.show()

_, ax = plt.subplots(figsize=(10, 7))
ax.plot(iter_list, objective, label='objective')
# ax.plot(iter_list, cons_failure, label='Failure')
ax.set_xlabel(r'Iteration $-$')
ax.set_ylabel(r'objective Value $-$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.show()