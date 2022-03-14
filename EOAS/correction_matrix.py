import numpy as np
from RethorstCorrection_pyf90.mod_overset_interpolation import overset_interpolation


def correction_(panels_VLM, nx, span, jet_loc, jet_radius, Vinf, Vjet, span_max, r_min):
    correction_matrix = np.zeros((panels_VLM, panels_VLM), order='F')
    mesh = np.zeros((panels_VLM+1, 3, nx), order='F')
    panels_jet = 91
    panels_overset_wing = 3101

    overset_interpolation(span, 1., jet_loc, jet_radius, Vinf, Vjet, panels_jet,\
        panels_overset_wing, nx, panels_VLM, span_max, r_min, correction_matrix, mesh)

    return correction_matrix, mesh


def correction_location(panels_VLM, span, jet_radius, jet_loc, y):
    y_ = np.zeros((panels_VLM))
    
    for index in range(panels_VLM):
        y_[index] = (y[index+1]+y[index])/2
    
    y_ = np.where(y_>(jet_loc+jet_radius), 0, y_)
    y_ = np.where(y_<(jet_loc-jet_radius), 0, y_)
    y_ = np.where(y_==0, 0, 1)

    return y_