import numpy as np
import matplotlib.pyplot as plt
from EOAS_system import EOAS_system_
import imageio

# span = 10m
steps = 15
span_max = 40
r_min = 0.5
Vjet = 150
Vinf = 100
jet_loc = 1.3
jet_radius = 0.52

filename = '/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/'

_ = EOAS_system_(jet_radius, jet_loc, Vinf, Vjet, r_min, span_max, filename)