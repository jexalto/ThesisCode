import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt

from RethorstCorrection_pyf90.mod_multiprop import multiprop


def correction_(panels_span_VLM, panels_chord_VLM,  panels_overset_wing, panels_jet, span, span_max, r_min, vinf,
                vel_distr_input, radii_input, prop_discr, jet_loc_list, nr_props, nr_radii_input):

    total_correction = np.zeros((panels_chord_VLM*panels_span_VLM, panels_chord_VLM*panels_span_VLM), order='F')
    y_VLM = np.zeros((panels_span_VLM+1), order='F')
    vel_vec = np.zeros((panels_span_VLM), order='F')

    multiprop(span, nr_props, jet_loc_list, vel_distr_input, radii_input, nr_radii_input, prop_discr, vinf, panels_jet, \
                panels_overset_wing, panels_chord_VLM, panels_span_VLM, span_max, r_min, vel_vec, total_correction)
        

    return total_correction, y_VLM, vel_vec


def correction_location(panels_VLM, span, jet_radius, jet_loc, y):
    y_ = np.zeros((panels_VLM))
    
    for index in range(panels_VLM):
        y_[index] = (y[index+1]+y[index])/2
    
    y_ = np.where(y_>(jet_loc+jet_radius), 0, y_)
    y_ = np.where(y_<(jet_loc-jet_radius), 0, y_)
    y_ = np.where(y_==0, 0, 1)

    return y_