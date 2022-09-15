import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt

from RethorstCorrection_pyf90.mod_overset_nointerpolation import overset_nointerpolation


def correction_(panels_VLM, panels_chord_VLM, span, jet_loc, jet_radius, Vinf, Vjet, span_max, r_min, panels_overset_wing, panels_jet):
    correction_matrix = np.zeros((panels_chord_VLM*panels_VLM, panels_chord_VLM*panels_VLM), order='F')
    mesh = np.zeros((panels_VLM+1, 3, 2), order='F')
    # ! The next values HAVE to be odd, also the wing/jet panels ratio needs to be larger than 100 !
    nx = panels_chord_VLM+1
    G_ = np.zeros(( (nx-1)*(panels_jet+2*panels_overset_wing), (nx-1)*(panels_jet+2*panels_overset_wing) ), order='F')
    discr_interp = np.zeros((panels_overset_wing*2+panels_jet), order='F')
    y_VLM = np.zeros((panels_VLM+1), order='F')
    
    Vjet = 149.90862405 # this corresponds to the 2 degree angle of attack
    overset_nointerpolation(span, 1., jet_loc, jet_radius, Vinf, Vjet, panels_jet,\
        panels_overset_wing, nx-1, panels_VLM, span_max, r_min, correction_matrix, mesh, G_, discr_interp, y_VLM)

#%% Spacing correction
    spacing_VLM_wing = span/panels_VLM
    spacing_overset_wing = abs(discr_interp[2]-discr_interp[1])
    spacing_overset_jet = abs(discr_interp[panels_overset_wing+4]-discr_interp[panels_overset_wing+3])
    delta_wing = spacing_VLM_wing/spacing_overset_wing
    delta_jet = spacing_VLM_wing/spacing_overset_jet

    index1 = np.argmin(abs(y_VLM-(jet_loc-jet_radius)))-1
    index2 = len(y_VLM)-np.argmin(abs(y_VLM[::-1]-(jet_loc+jet_radius)))-1
    
    # --- This part is necessary if a single moving node is being used ---
    spacing_VLM_wing_left = abs(y_VLM[index1+1]-y_VLM[index1])
    spacing_VLM_wing_right = abs(y_VLM[index2+1]-y_VLM[index2])

    delta_wing_left = spacing_VLM_wing_left/spacing_overset_wing
    delta_wing_right = spacing_VLM_wing_right/spacing_overset_wing

    spacing_VLM_prop_left = abs(y_VLM[index1+2]-y_VLM[index1+1])
    spacing_VLM_prop_right = abs(y_VLM[index2-1]-y_VLM[index2])

    delta_prop_left = spacing_VLM_prop_left/spacing_overset_jet
    delta_prop_right = spacing_VLM_prop_right/spacing_overset_jet

    correction_matrix[:, 0:index1] = correction_matrix[:, 0:index1]*delta_wing
    correction_matrix[:, index1] = correction_matrix[:, index1]*delta_wing_left
    correction_matrix[:, index1+1] = correction_matrix[:, index1+1]*delta_prop_left

    correction_matrix[:, index1+2:index2-1] = correction_matrix[:, index1+2:index2-1]*delta_jet

    correction_matrix[:, index2-1] = correction_matrix[:, index2-1]*delta_prop_right
    correction_matrix[:, index2] = correction_matrix[:, index2]*delta_wing_right
    correction_matrix[:, index2+1:] = correction_matrix[:, index2+1:]*delta_wing

    # --- This part is necessary if all nodes are slightly adjusted ---
    # correction_matrix[:, 0:index1+1] = correction_matrix[:, 0:index1+1]*delta_wing
    # correction_matrix[:, index1+1:index2+1] = correction_matrix[:, index1+1:index2+1]*delta_jet
    # correction_matrix[:, index2+1:] = correction_matrix[:, index2+1:]*delta_wing

    discr_interp_ = discr_interp
    interpolation = sp.interpolate.interp2d(discr_interp_, discr_interp_, G_, kind='linear')

    y_ = np.zeros((panels_VLM))
    for i in range(panels_VLM):
        y_[i] = (y_VLM[i+1]+y_VLM[i])/2

    correction_scipy = np.zeros((len(y_), len(y_)))
    for index_i, element_i in enumerate(y_):
        for index_j in range(panels_VLM-1):
            if y_[index_j]>=jet_loc:
                break
            correction_scipy[index_i, index_j] = interpolation(y_[index_j], element_i)

    # --- artificial peak ---
    peak = np.zeros((len(discr_interp)))
    
    for index, _ in enumerate(peak):
        if G_[index, panels_overset_wing+round(panels_jet/2)]>1.1*G_[index, panels_overset_wing+round(panels_jet/2)-1]:
            peak[index] = 1

    peak_interpolated = np.zeros((len(y_)))

    for index, _ in enumerate(peak_interpolated):
        if correction_scipy[index, index_j-1]<1.1*correction_scipy[index, index_j-2]:
            correction_scipy[index, index_j-1] = correction_scipy[index, index_j-1]*2
            peak_interpolated[index] = 1
        
    panels_left = panels_VLM-index_j
    addition = correction_scipy[:, index_j-panels_left-1:index_j-1]
    correction_scipy[:, index_j:] = addition[:, ::-1]

    point_orig = panels_overset_wing*2
    plt.clf()
    plt.figure(1)
    plt.plot(np.arange(0, len(G_[point_orig,:]), 1), G_[point_orig, :], label='original data')
    # plt.plot(np.arange(0, len(correction_scipy[0,:]), 1), correction_scipy[0,:], label='interpolated')
    plt.legend()
    plt.grid()
    plt.title(f'Scaled original data - point {point_orig}')
    plt.xlabel('Spanwise point')
    plt.ylabel('Correction value')
    plt.xlim((400, 700))
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/peaks/original_data.png')

    plt.clf()
    plt.figure(1)
    # plt.plot(np.arange(0, len(G_[0,:]), 1), G_[15, :], label='original data')
    plt.plot(np.arange(0, len(correction_scipy[0,:]), 1), correction_scipy[0,:], label='interpolated')
    plt.legend()
    plt.grid()
    plt.title(f'VLm panels: {panels_VLM} Overset Wing panels: {panels_overset_wing}, jet panels: {panels_jet}')
    plt.xlabel('Spanwise point')
    plt.ylabel('Correction value')
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/jexalto/code/MDO_lab_env/ThesisCode/EOAS/figures/peaks/correctioncheck.png')    

    return correction_scipy, panels_jet, panels_overset_wing, G_, discr_interp, y_VLM, delta_wing, delta_jet


def correction_location(panels_VLM, span, jet_radius, jet_loc, y):
    y_ = np.zeros((panels_VLM))
    
    for index in range(panels_VLM):
        y_[index] = (y[index+1]+y[index])/2
    
    y_ = np.where(y_>(jet_loc+jet_radius), 0, y_)
    y_ = np.where(y_<(jet_loc-jet_radius), 0, y_)
    y_ = np.where(y_==0, 0, 1)

    return y_