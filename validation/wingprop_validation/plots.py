import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import niceplots

from wingprop_analysis import wingprop
from wing_analysis import wing


file = '/home/mdolabuser/mount/ThesisCode/HELIX_verification/TUD_data/prowim/SecIIIC_ModelII/SecIIIC_Fig24_CLCD_ModelII_conventional_ReD640k.txt'

data = pd.read_csv(file, delimiter=',', skiprows=22)

# n=0
# index1 = n*19
# index2 = (n+1)*19
# aoa = data['AoA'][index1:index2]
# CL_Jinf = data['CL'][index1:index2]
# CD_Jinf = data['CD'][index1:index2]
# J_inf = data['J'][index1+1]

# alphas, CL_num_Jinf, CD_num_Jinf = wing()

n=1
index1 = n*19
index2 = (n+1)*19
aoa = data['AoA'][index1:index2]
CL_J1 = data['CL'][index1:index2]
CD_J1 = data['CD'][index1:index2]
J_1 = data['J'][index1+1]

alphas, CL_num_J1, CD_num_J1 = wingprop(J=1.)

n=4
index1 = n*19
index2 = (n+1)*19
aoa = data['AoA'][index1:index2]
CL_J0796 = data['CL'][index1:index2]
CD_J0796 = data['CD'][index1:index2]
J_0796 = data['J'][index1+1]

alphas, CL_num_J0796, CD_num_J0796 = wingprop(J=0.7962)

niceplots.setRCParams()

_, ax = plt.subplots(figsize=(10, 7))
ax.scatter(alphas, CL_num_Jinf, label='Numerical, prop-off')
ax.scatter(alphas, CL_num_J1, label=f'Numerical, J=1.0')
ax.scatter(alphas, CL_num_J0796, label=f'Numerical, J=0.7962')
# ax.scatter(aoa, CL_Jinf, label=f'Experimental, prop-off')
ax.set_xlabel("Angle of Attack (deg)")
ax.set_ylabel(r"Lift Coefficient ($C_L$)")
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig(f'/home/jexalto99/code/MDO_lab_env/ThesisCode/validation/wingprop_validation/00_results/figures/numericaldata.png')

# _, ax = plt.subplots(figsize=(10, 7))
# ax.plot(alphas, CL_num_J1, label=f'Numerical, J=1.0')
# ax.scatter(aoa, CL_J1, label=f'Experimental, J=1.0')
# ax.set_xlabel("Angle of Attack (deg)")
# ax.set_ylabel(r"Lift Coefficient ($C_L$)")
# ax.legend()
# ax.grid()
# niceplots.adjust_spines(ax, outward=True)
# plt.savefig(f'/home/jexalto99/code/MDO_lab_env/ThesisCode/validation/wingprop_validation/00_results/figures/wingprop_validation_J1.0.png')

# _, ax = plt.subplots(figsize=(10, 7))
# ax.plot(alphas, CL_num_J0796, label=f'Numerical, J=0.7962')
# ax.scatter(aoa, CL_J0796, label=f'Experimental, J=0.7962')
# ax.set_xlabel("Angle of Attack (deg)")
# ax.set_ylabel(r"Lift Coefficient ($C_L$)")
# ax.legend()
# ax.grid()
# niceplots.adjust_spines(ax, outward=True)
# plt.savefig(f'/home/jexalto99/code/MDO_lab_env/ThesisCode/validation/wingprop_validation/00_results/figures/wingprop_validation_J07962.png')