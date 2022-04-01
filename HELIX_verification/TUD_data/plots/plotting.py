from asyncore import write
import numpy as np
import niceplots
import matplotlib.pyplot as plt
import pandas as pd
import json as js

path = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/TUD_data/prowim/SecIIIA_StingMountedPropeller/SecIIIA_Fig8a_PerfoData_StingProp_AoA0_ReD620k.txt'

file = pd.read_csv(path, header=22, usecols=['Polar', 'Run', 'AoA', 'J=Vinf/nD', 'CT=T/rho*n2*D4', 'CP=P/rho*n3*D5', 'eta=J*CT/CP'], sep=',')
niceColors = niceplots.get_niceColors()

J       = file['J=Vinf/nD']
CT      = file['CT=T/rho*n2*D4']
CP      = file['CP=P/rho*n3*D5']
eta     = file['eta=J*CT/CP']

savedir = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/TUD_data/plots/figure/prowim.png'


niceplots.setRCParams(dark_mode=False, set_dark_background=False)

fig, ax = plt.subplots()
ax.scatter(J, CT, label='CT', marker='o', s=80)
ax.scatter(J, CP, label='CP', marker='x', s=80)
ax.scatter(J, eta, label='eta', marker='*', s=80)
ax.set_xlabel("$J$")
ax.set_ylabel("$C_T/C_P/\eta$", rotation="vertical")
ax.grid()
# ax.set_xticks(np.linspace(0, 2, 5) * np.pi)
fig.legend()
fig.savefig(savedir, dpi=400)
plt.gcf()
plt.close(fig)

def writeJSON(data, fileName):
    """
    This function generates and writes a JSON version of the rotor parameters
    data dictionary to a .json file.
    Parameters
    ----------
    data : dictionary
        Dictionary of rotor sectional and spanwise data
    fileName : string
        Name of file to be output as JSON geometry file
    """
    
    with open(fileName, "w") as file:
        js.dump(data, file, indent=4)

data = {}
jsonsave = '/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/data/prowim_data.json'
data['J']       = J.tolist()
data['CT']      = CT.tolist()
data['CP']      = CP.tolist()
data['eta']     = eta.tolist()

for i in range(len(data['J'])):
    print(data['J'][i], data['CT'][i])
    print('-----------')
writeJSON(data, jsonsave)