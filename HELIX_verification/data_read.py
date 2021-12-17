import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import niceplots

data = pd.read_csv('/Users/joaquinexalto/Documents/TU_Delft/code/HELIX_verification/tud_data.csv', dtype=np.float16, sep=',', skiprows=(2), header=None).values

# Polar,  Run,    AoA,    Vinf,   rhoInf,     Tinf,   Pinf,       n,      J=Vinf/nD,  CT=T/rho*n2*D4, CP=P/rho*n3*D5, eta=J*CT/CP

J = data[:, 8]
CT = data[:, 9]
CP = data[:, 10]
eta = data[:, 11]

_, (ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(25, 7))
ax1.scatter(J, CT)

# Plot HELIX Result
ax1.set_xlim([0.6,  1])
ax1.set_ylim([0,    0.16])
ax1.set_xlabel("Advance Ratio [-]")
ax1.set_ylabel(r"Thrust coefficient $[-]$")
niceplots.adjust_spines(ax1, outward=True)

ax2.scatter(J, CP)

# Plot HELIX Result
ax2.set_xlim([0.6,  1])
ax2.set_ylim([0,    0.16])
ax2.set_xlabel("Advance Ratio [-]")
ax2.set_ylabel(r"Thrust coefficient $[-]$")
niceplots.adjust_spines(ax2, outward=True)

ax3.scatter(J, eta)

# Plot HELIX Result
ax3.set_xlim([0.6,  1])
ax3.set_ylim([0,    0.16])
ax3.set_xlabel("Advance Ratio [-]")
ax3.set_ylabel(r"Thrust coefficient $[-]$")
niceplots.adjust_spines(ax3, outward=True)

plt.show()