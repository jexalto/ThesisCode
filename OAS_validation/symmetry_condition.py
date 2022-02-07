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

Vinf = 100
Vjet = 105

ny = 12
nx = 2

y = np.linspace(0, 10, ny)
y = (y[1:] + y[0:-1])/2

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
y_n = abs(y-jet_loc)/jet_radius # if y_n<1 --> inside jet
if right:
    y_n = y_n[n_:]
    c = y_n[0:-1]
    d = y_n[1:]
else:
    y_n = y_n[:N]
    c = np.flip(y_n[1:])
    d = np.flip(y_n[0:-1])

# %% DEFINE INPUTS FOR RETHORST CORRECTION PACKAGE
mu = Vinf/Vjet
ksi = np.ones( (len(y_n)) ) * 1.1
eta = y_n # z/r
radius = 1.0

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
                # gjj_odd(mu, ksi, eta, c, d, radius, gjj_out)
                correction[element, element2] = 0.5 # gjj_out
            else:
                # JO
                # gjo_odd(mu, ksi, eta, c, d, radius, gjo_out)
                correction[element, element2] = 0.51 # gjo_out
    else:
        # O-
        for element2 in range(N):
            if y_jet_bool[element2]==1:
                # OJ
                # goj_odd(mu, ksi, eta, c, d, radius, goj_out)
                correction[element, element2] = 1.5 # goj_out
            else:
                # OO
                # goo_odd(mu, ksi, eta, c, d, radius, goo_out)
                correction[element, element2] = 1.51 # goo_out

G_ = np.zeros(( (ny-1)*(nx-1), (ny-1)*(nx-1) ))
G_1 = np.zeros(( (ny-1)*(nx-1), (ny-1)*(nx-1) ))
right = True
if right:
    for i in range(nx):
        G_[(ny)*i - i : ny*(i+1)-(i+1), n_ + (ny-1)*i : ny*(i+1)-(i+1)] = 1
if right:
    for i in range(nx):
        G_1[(ny)*i - i : ny*(i+1)-(i+1), (ny-1)*i : n_ + (ny-1)*i] = 1

print(G_ + G_1)


# goo_odd(mu, ksi, eta, c, d, radius, goo_out)
# goj_odd(mu, ksi, eta, c, d, radius, goj_out)
# gjo_odd(mu, ksi, eta, c, d, radius, gjo_out)
# gjj_odd(mu, ksi, eta, c, d, radius, gjj_out)
