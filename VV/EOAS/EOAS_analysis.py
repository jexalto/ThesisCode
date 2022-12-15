import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_

# span = 10m
steps = 15
span_max = 5
r_min = 0.1
Vjet = 150
Vinf = 100
jet_loc_list = np.array([0.], order='F')
jet_radius = 0.237/2
span=0.748

filename = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/'

nx_input = 2
prop_discr = 5
nr_props = len(jet_loc_list)

_ = EOAS_system_(jet_radius, jet_loc_list, Vinf, r_min, span_max, filename, nx_input, prop_discr, span, nr_props)