import numpy as np
import matplotlib.pyplot as plt

ny=71
y = np.linspace(-5, 5, ny)
y_ = np.zeros((ny-1))

for i in range(ny-1):
    y_[i] = (y[i+1]+y[i])/2

height = 0.1
jet_loc = 0.
jet_radius = 1.0
dir = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures'

plt.scatter(y_, np.ones((ny-1))*height, marker='x', label='VLM mesh')
plt.scatter([jet_loc-jet_radius, jet_loc, jet_loc+jet_radius], np.zeros((3)), marker='|', label='propeller')
plt.legend()
plt.grid()
plt.ylim((-0.025, 0.125))
plt.xlim((-1.1, -0.9))
plt.savefig(dir+'/oversetmesh.png')