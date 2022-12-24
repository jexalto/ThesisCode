import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import niceplots

from wingprop_analysis import wingprop
from wing_analysis import wing

file = '/home/jexalto99/code/MDO_lab_env/HELIX_verification/TUD_data/prowim/SecIIIC_ModelII/SecIIIC_Fig24_CLCD_ModelII_conventional_ReD640k.txt'

data = pd.read_csv(file, delimiter=',', skiprows=22)

n=0
index1 = n*19
index2 = (n+1)*19
aoa = data['AoA'][index1:index2]
CL_Jinf = data['CL'][index1:index2]
CD_Jinf = data['CD'][index1:index2]
J_inf = data['J'][index1+1]

# alphas, CL_num_Jinf, CD_num_Jinf = wing()

n=1
index1 = n*19
index2 = (n+1)*19
aoa = data['AoA'][index1:index2]
CL_J1 = data['CL'][index1:index2]
CD_J1 = data['CD'][index1:index2]
J_1 = data['J'][index1+1]

# alphas, CL_num_J1, CD_num_J1 = wingprop(J=1.)

n=4
index1 = n*19
index2 = (n+1)*19
aoa = data['AoA'][index1:index2]
CL_J0796 = data['CL'][index1:index2]
CD_J0796 = data['CD'][index1:index2]
J_0796 = data['J'][index1+1]

# alphas, CL_num_J0796, CD_num_J0796 = wingprop(J=0.7962)

niceplots.setRCParams()

_, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18, 16))
# ax1.plot(alphas, CL_num_Jinf, label='Numerical, prop-off')
# ax1.plot(alphas, CL_num_J1, label=f'Numerical, J=1.0')
# ax1.plot(alphas, CL_num_J0796, label=f'Numerical, J=0.7962')
ax1.scatter(aoa, CL_Jinf, label=f'Experimental, prop-off')
ax1.set_xlabel("Angle of Attack (deg)")
ax1.set_ylabel(r"Lift Coefficient ($C_L$)")
ax1.legend()
# ax1.grid()
niceplots.adjust_spines(ax1, outward=True)
# plt.savefig(f'00_results/figures/numericaldata.png')

# _, ax = plt.subplots(figsize=(10, 7))
# ax2.plot(alphas, CL_num_J1, label=f'Numerical, J=1.0')
ax2.scatter(aoa, CL_J1, label=f'Experimental, J=1.0')
ax2.set_xlabel("Angle of Attack (deg)")
ax2.set_ylabel(r"Lift Coefficient ($C_L$)")
ax2.legend()
# ax2.grid()
niceplots.adjust_spines(ax2, outward=True)
# plt.savefig(f'00_results/figures/wingprop_validation_J1.0.png')

# _, ax = plt.subplots(figsize=(10, 7))
# ax3.plot(alphas, CL_num_J0796, label=f'Numerical, J=0.7962')
ax3.scatter(aoa, CL_J0796, label=f'Experimental, J=0.7962')
ax3.set_xlabel("Angle of Attack (deg)")
ax3.set_ylabel(r"Lift Coefficient ($C_L$)")
ax3.legend()
# ax3.grid()
niceplots.adjust_spines(ax3, outward=True)
# plt.savefig(f'00_results/figures/wingprop_validation_J07962.png')
plt.savefig('00_results/stacked.png')