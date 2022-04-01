import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_

# span = 10m
steps = 15
span_max = 40
r_min = 0.4
Vjet = 150
Vinf = 100
jet_loc = -1.0
jet_radius = 1.0

filename = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/'

nx = 2

_ = EOAS_system_(jet_radius, jet_loc, Vinf, Vjet, r_min, span_max, filename, nx, 10)