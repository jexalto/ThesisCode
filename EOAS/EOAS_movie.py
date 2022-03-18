import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_
import imageio

steps = 15
span_max = 40
r_min = 0.5
Vjet = 150
Vinf = 100
jet_loc = 1.
jet_radius = 0.5

jet_radii = np.linspace(r_min, 2., steps)
jet_locs = np.linspace(-3.6, 3.8, steps)

filename_loc = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/loc/'
filename_radii = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/radii/'

filelist_radii = []
filelist_locs = []

# --- Radius gif ---
if True:
    for index, element in enumerate(jet_radii):
        print(f'jetradius: {element}')
        filelist_radii.append(EOAS_system_(element, jet_loc, Vinf, Vjet, r_min, span_max, filename_radii))

    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/radii.gif', mode='I') as writer:
        for filename in filelist_radii:
            for i in range(2):      # slow down the gif a bit
                image = imageio.imread(filename)
                writer.append_data(image)

# --- Loc gif ---
if False:
    for index, element in enumerate(jet_locs):
        print(f'jetloc: {element}')
        filelist_locs.append(EOAS_system_(jet_radius, element, Vinf, Vjet, r_min, span_max, filename_loc))

    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/locs.gif', mode='I') as writer:
        for filename in filelist_locs:
            for i in range(2):      # slow down the gif a bit
                image = imageio.imread(filename)
                writer.append_data(image)