import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_

# span = 10m
steps = 15
span_max = 5
r_min = 0.1
Vjet = 150
Vinf = 100
jet_loc = 0.1
jet_radius = 0.237/2
span=0.748

filename = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/'

nx = 2

_ = EOAS_system_(jet_radius, jet_loc, Vinf, Vjet, r_min, span_max, filename, nx, 5, span)