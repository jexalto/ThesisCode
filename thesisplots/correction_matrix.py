import numpy as np
import scipy.interpolate as sp
import matplotlib.pyplot as plt

import niceplots

from RethorstCorrection_pyf90.mod_velocity_distribution_nonsymmetry import velocity_distribution_nosym


def correction_(    span, jet_loc_input, vel_distr_input, radii_input, prop_discr, vinf, panels_jet,
                    panels_overset_wing, panels_chord_vlm, panels_span_vlm, span_max, r_min):

    total_correction = np.zeros((panels_chord_vlm*panels_span_vlm, panels_chord_vlm*panels_span_vlm), order='F')
    vel_vec = np.zeros((panels_span_vlm), order='F')

    velocity_distribution_nosym(    span, jet_loc_input, vel_distr_input, radii_input, prop_discr, vinf, panels_jet,
                                    panels_overset_wing, panels_chord_vlm, panels_span_vlm, span_max, r_min, vel_vec, total_correction)
        

    return total_correction, vel_vec

panels_overset_wing = 2501
panels_jet = 151

steps = 10
step = np.arange(0, steps, 1)
jet_radius=1.0
def x2(steps):
    return 150-(0.001*steps-3)**2
span=10.
jet_loc_input=0.0
vel_distr_input = np.array(x2(step), order='F')
radii_input = np.array(np.linspace(0.001, jet_radius, steps), order='F')
prop_discr=20
vinf=100.
panels_chord_vlm=1
panels_span_vlm_low=100
panels_span_vlm=200
span_max=20.
r_min=0.1
total_correction_high, vel_vec = correction_(    span, jet_loc_input, vel_distr_input, radii_input, prop_discr, vinf, panels_jet,
                                                panels_overset_wing, panels_chord_vlm, panels_span_vlm, span_max, r_min)
x_high = np.arange(0, len(total_correction_high[0,:]), 1)

total_correction_low, vel_vec = correction_(    span, jet_loc_input, vel_distr_input, radii_input, prop_discr, vinf, panels_jet,
                                                panels_overset_wing, panels_chord_vlm, panels_span_vlm_low, span_max, r_min)
x_low = np.arange(0, len(total_correction_low[0,:]), 1)

dir = '/home/jexalto/code/MDO_lab_env/ThesisCode/thesisplots/figures/'
index_low = 0#round(panels_span_vlm_low/4)
index_high = 0#round(panels_span_vlm/4)

niceplots.setRCParams()
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

fig, ax = plt.subplots(nrows=2, figsize=(16, 12))

line = ax[0].plot(x_low, total_correction_low[index_low,:], clip_on=False, color=colors[0])
# ax.vlines(tp, -3, 1.0 - np.cos(omega * tp), linestyle="--", color="gray", zorder=0)
ax[0].set_ylabel("Correction factor", ha="right", rotation="vertical")
ax[0].set_title(f'VLM panels={panels_span_vlm_low}')
ax[0].set_ylim(bottom=1.2*min(total_correction_low[index_low,:]), top=1.2*max(total_correction_low[index_low,:]))
ax[0].grid()

line = ax[1].plot(x_high, total_correction_high[index_high,:], clip_on=False, color=colors[0])
# ax.vlines(tp, -3, 1.0 - np.cos(omega * tp), linestyle="--", color="gray", zorder=0)
ax[1].set_ylabel("Correction factor", ha="right", rotation="vertical")
ax[1].set_xlabel("Spanwise panels", ha="right", rotation="horizontal")
ax[1].set_title(f'VLM panels={panels_span_vlm}')
ax[1].set_ylim(bottom=1.2*min(total_correction_high[index_high,:]), top=1.2*max(total_correction_high[index_high,:]))
ax[1].grid()
niceplots.adjust_spines(ax[0], outward=True)
niceplots.adjust_spines(ax[1], outward=True)

plt.savefig(dir+"correction.png", dpi=400)