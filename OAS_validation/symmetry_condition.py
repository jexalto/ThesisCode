import numpy as np
import scipy.special as sp
from matplotlib import pyplot as plt

from RethorstCorrection_pyf90.mod_gvalues_odd import goo_odd, gjo_odd, goj_odd, gjj_odd

from RethorstCorrection_pyf90.mod_integrals import wjo_int, woo_int, kv_integral

from RethorstCorrection_pyf90.mod_modifiedbessel import secondkind_bessel, firstkind_bessel, firstkind_bessel_der, secondkind_bessel_der

# %% Parameters
span = 10
jet_loc = 4.0
jet_radius = 0.75

span = 10.
crd = 2.

Vinf = 100
Vjet = 105

ny = 12
nx = 2

y_ = np.linspace(0, span, ny)
y = (y_[1:] + y_[0:-1])/2

chord_ = np.linspace(0, crd, nx)
chord = (chord_[1:] + chord_[0:-1])/2

#%% WHAT SIDE OF THE WING DO WE NEED TO CALCULATE?
right = span - jet_loc
left = span - right

if right>left: # How many points do we need to calculate the correction for?
    right = True
    N_ = np.where(y>jet_loc, 1, 0)
    N = sum(N_)
else:
    right = False
    N_ = np.where(y<jet_loc, 1, 0)
    N = sum(N_)

n_ = len(y) - N # how many points need the symmetry condition
n = n_ - 1 # index of point that needs correction

# %% ARE POINTS INSIDE/OUTSIDE JET?
if right:
    y_jet = y[n_:] - (jet_loc + jet_radius)
    y_jet_bool = np.where(y_jet<0, 1, 0)
else:
    y_jet = y[:N] - (jet_loc - jet_radius)
    y_jet_bool = np.where(y_jet>0, 1, 0)

# --- Make sure we get the normalised y value (normalised over jet radius) ---
y_n_ = abs(y_-jet_loc)/jet_radius # if y_n<1 --> inside jet
if right:
    y_n_ = y_n_[n_:]
    c = y_n_[0:-1]
    d = y_n_[1:]
    eta = abs(y[n_:]-jet_loc)/jet_radius # z/r
else:
    y_n_ = y_n_[:N]
    c = np.flip(y_n_[1:])
    d = np.flip(y_n_[0:-1])
    eta = abs(y[:N]-jet_loc)/jet_radius # z/r

# %% DEFINE INPUTS FOR RETHORST CORRECTION PACKAGE
mu = Vinf/Vjet
ksi = chord
radius = jet_radius

# --- Make entire correction matrix ---
# First calculate corrections
correction = np.zeros((N, N))
goo_out = np.array(0.0)
goj_out = np.array(0.0)
gjo_out = np.array(0.0)
gjj_out = np.array(0.0)

for element in range(N):
    if y_jet_bool[element]==1:
        # J-
        for element2 in range(N):
            if y_jet_bool[element2]==1:
                # JJ
                c_ = c[element2]
                d_ = d[element2]
                eta_ = eta[element2]
                ksi = chord[int( (element-0.1)/(ny-1) )]
                # gjj_odd(mu, ksi, eta, c, d, radius, gjj_out)
                correction[element, element2] = element+element2 # 0.5 # gjj_out
            else:
                # JO
                c_ = c[element2]
                d_ = d[element2]
                eta_ = eta[element2]
                ksi = chord[int( (element-0.1)/(ny-1) )]
                # gjo_odd(mu, ksi, eta, c, d, radius, gjo_out)
                correction[element, element2] = element+element2  # 0.51 # gjo_out
    else:
        # O-
        for element2 in range(N):
            if y_jet_bool[element2]==1:
                # OJ
                c_ = c[element2]
                d_ = d[element2]
                eta_ = eta[element2]
                ksi = chord[int( (element-0.1)/(ny-1) )]
                # goj_odd(mu, ksi, eta, c, d, radius, goj_out)
                correction[element, element2] = element+element2 #  1.5 # goj_out
            else:
                # OO
                c_ = c[element2]
                d_ = d[element2]
                eta_ = eta[element2]
                ksi = chord[int( (element-0.1)/(ny-1) )]
                # goo_odd(mu, ksi, eta, c, d, radius, goo_out)
                correction[element, element2] = element+element2  # goo_out

correction_ = np.zeros((ny-1, ny-1))
if right:
    correction_[0:n_, 0:n_+1] = np.flip(correction[1:n_+1, 0:n_+1])/100
    correction_[n_:, n_:] = correction[:, :]/100

G_ = np.zeros(( (ny-1)*(nx-1), (ny-1)*(nx-1) ))

right = True
if right:
    for i in range(nx):
        G_[(ny)*i - i : ny*(i+1)-(i+1), n_ + (ny-1)*i : ny*(i+1)-(i+1)] = 1
else:
    for i in range(nx):
        G_[(ny)*i - i : ny*(i+1)-(i+1), (ny-1)*i : n_ + (ny-1)*i] = 1
a = None
