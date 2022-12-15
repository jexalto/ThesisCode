import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_
import imageio
import os

steps = 15
span_max = 40
r_min = 0.5
Vjet = 150
Vinf = 100
jet_loc = 0.
jet_radius = 1.0

jet_radii = np.linspace(r_min, 2., steps)
jet_locs = np.linspace(0.1, 3.8, steps)
prop_discr = np.arange(5, 25, 1)
nx_ = [3, 5, 7]

filename_loc = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/loc/'
filename_radii = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/radii/'
filename_nx = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/nx/'
filename_propdiscr = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/prop_discr/'

filelist_radii = []
filelist_locs = []
filelist_nx = []
filelist_prop_discr = []

if False:
    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/prop_discr.gif', mode='I') as writer:
        for filename in os.listdir(filename_propdiscr):
            f = os.path.join(filename_propdiscr, filename)
            for i in range(2):      # slow down the gif a bit
                image = imageio.imread(f)
                writer.append_data(image)

# --- Radius gif ---
if False:
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
        filelist_locs.append(EOAS_system_(jet_radius, element, Vinf, Vjet, r_min, span_max, filename_loc))

    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/locs.gif', mode='I') as writer:
        for filename in filelist_locs:
            for i in range(2):      # slow down the gif a bit
                image = imageio.imread(filename)
                writer.append_data(image)

# --- nx gif ---
if False:
    for index, element in enumerate(nx_):
        print(f'nx: {element}')
        filelist_nx.append(EOAS_system_(jet_radius, jet_loc, Vinf, Vjet, r_min, span_max, filename_nx, element))

    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/nx.gif', mode='I') as writer:
        for filename in filelist_nx:
            for i in range(2):      # slow down the gif a bit
                image = imageio.imread(filename)
                writer.append_data(image)

# --- prop_discr gif ---
if True:
    for index, element in enumerate(prop_discr):
        filelist_prop_discr.append(EOAS_system_(jet_radius, jet_loc, Vinf, Vjet, r_min, span_max, filename_propdiscr, 3, element))

    with imageio.get_writer('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/movie/gif_lib/prop_discr.gif', mode='I') as writer:
        for filename in filelist_prop_discr:
            for i in range(3):      # slow down the gif a bit
                image = imageio.imread(filename)
                writer.append_data(image)