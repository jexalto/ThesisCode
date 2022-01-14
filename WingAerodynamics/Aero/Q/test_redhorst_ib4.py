import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import sys
import pandas as pd

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from vortices import v_induced_by_horseshoe_vortex
from jet_correction_ib2 import jet_correction_ib

#%% INPUTS
# code for two jets at the wingtips

Vinf = 100              #velocity [m/s]
alpha = np.deg2rad(2)   #aoa [rad]

r0 = 1                  #jet radius [m]
mu = 2/3                #Vinf/Vjet [-]
crd = 1                 #wing chord [m]
A = 10                  #wing aspect ratio [m]
y_jet = 2
tc = 1.06               #wing thickness correction [-]

nr = 8                   #number of panels in jet
nl = 4                  #number of panels outside jet
ni = 24

Vj = Vinf/mu            #Vjet [m/s]
b = A*crd               #wingspan [m]
x = -crd*1/2            #control point x location [m]
dzeta = -crd*1/2/r0     #non-dim control point x [-]

#%% Calculate geometry

#left jet of
yl_ = np.linspace(0, y_jet-r0, nl)
y_l = (yl_[1:]+yl_[:-1])/2
s_l = (yl_[1:]-yl_[:-1])/2

#inside jet
yp_ = np.linspace(y_jet-r0, y_jet+r0, ni)
y_p = (yp_[1:]+yp_[:-1])/2
s_p = (yp_[1:]-yp_[:-1])/2

#right from jet
yr_ = np.linspace(y_jet+r0, b/2, nr)
y_r = (yr_[1:]+yr_[:-1])/2
s_r = (yr_[1:]-yr_[:-1])/2

#full wing
s = np.concatenate((s_r, s_p, s_r))
y = np.concatenate((y_l, y_p, y_r))
y_ = np.concatenate((yl_, yp_[1:], yr_[1:]))

s = np.concatenate((-s[::-1],s))
y = np.concatenate((-y[::-1],y))
y_ = np.concatenate((-y_[::-1],y_[1:]))

#number of panels
N = len(y)

#%% Velocity

V = np.ones(y.size)*Vinf
# V = np.where(np.abs(y-y_jet)<r0, Vj, Vinf)
# V = np.where(np.abs((np.abs(y)-y_jet)/r0)<1, Vj, Vinf)
jet =1 # number of jets

for panel in range(N):  # if panel in jet
    if jet == 1:
        if y_jet - r0 < y[panel] < y_jet+r0:
            xi = np.abs(y_jet - np.abs(y[panel]))
            V[panel] = Vinf - (Vj - Vinf) * 6 * (xi - 1) * xi

    if jet==2:
        if np.abs((np.abs(y[panel])-y_jet)/r0)<1:
            xi = np.abs(y_jet-np.abs(y[panel]))
            V[panel] = Vinf-(Vj-Vinf)*6*(xi-1)*xi

#V = np.where(np.abs((np.abs(y)-y_jet)/r0)<1, Vj, Vinf)
#V = np.where(np.abs(np.abs(y)-y_jet)/r0<1, Vj, Vinf)

#%% Horseshoe induced velocity

F = np.zeros((N,N))

for cp in range(N):
    for panel in range(N):
        P = np.array([-x, y[cp]])       #control point
        A = np.array([0, y_[panel]])    #point 1 horseshoe
        B = np.array([0, y_[panel+1]])  #point 2 horseshoe
        v = v_induced_by_horseshoe_vortex(P,A,B)
        v = -v[0]*(4*np.pi)
        F[cp, panel] = v
           
y_sym = y_[y_>=0]
V_sym = np.interp(y_sym, y, V)

#%% Jet correction

G_ = 0

for i in range(y_sym.size-1):
    mu = V_sym[i-1]/V_sym[i]
    loc = y_sym[i]
    if mu != 1 and loc<y_jet:
        Gi = jet_correction_ib(y_, mu, crd, loc, y_jet)
        G_ += Gi

#%% Apply jet correction to wing and solve

G = np.zeros((N, N))

#mirror the influence matrix G to describe the correction on both sides
#of the jet
N_ = int(N/2)
#G[:N_, :N_] = G_
#G[N_:, N_:] = G_[::-1, ::-1]
G = G_
if jet == 2:
    G = G_ + G[::-1,::-1]

#calculate circulation with jet
RHS = V*alpha
AIC = 1/(4*np.pi)*(F+G)
gamma = np.linalg.solve(AIC, RHS)
cl = 2*gamma*V/(Vinf**2*crd)

#calcualte with uniform jet correction
y_t = np.array([-5., -4.64705882, -4.29411765, -3.94117647, -3.58823529,
              -3.23529412, -2.88235294, -2.52941176, -2.17647059, -1.82352941,
              -1.47058824, -1.11764706, -0.76470588, -0.41176471, -0.05882353,
              0.29411765, 0.64705882, 1., 1.22222222, 1.44444444,
              1.66666667, 1.88888889, 2.11111111, 2.33333333, 2.55555556,
              2.77777778, 3., 3.33333333, 3.66666667, 4.,
              4.33333333, 4.66666667, 5.])
G_u = jet_correction_ib(y_, Vinf/Vj, crd, y_jet-r0, y_jet)
Gu = np.zeros((N, N))
Gu = G_u
if jet ==2:
    Gu = Gu + Gu[::-1, ::-1]
# Gu[:N_, :N_] = G_u[::-1, ::-1]
# Gu[N_:, N_:] = G_u
AIC_u = 1/(4*np.pi)*(F+Gu)
gamma_u = np.linalg.solve(AIC_u, RHS)
cl_u = 2*gamma_u*V/(Vinf**2*crd)

#calculate circulation with jet without rethorst correction
AIC_n = 1/(4*np.pi)*(F)
gamma_n = np.linalg.solve(AIC_n, RHS)
cl_n = 2*gamma_n*V/(Vinf**2*crd)

#calculate circulation without jet
AIC_v = 1/(4*np.pi)*(F)
RHS_v = Vinf*alpha*np.ones(RHS.shape)
gamma_v = np.linalg.solve(AIC_v, RHS_v)
cl_v = 2*gamma_v/(Vinf*crd)

#calculate circulation with uniform jet velocity and correction
#Vu = np.where(np.abs(y-y_jet)<r0, Vj, Vinf)
if jet ==1:
    Vu = np.where(np.abs(y-y_jet)/r0<1, Vj, Vinf)
elif jet ==2:
    Vu = np.where(np.abs(np.abs(y) - y_jet) / r0 < 1, Vj, Vinf)

RHS_uu = Vu*alpha
AIC_uu = 1/(4*np.pi)*(F+Gu)
gamma_uu = np.linalg.solve(AIC_uu, RHS_uu)
cl_uu = 2*gamma_uu*Vu/(Vinf**2*crd)



   
#create figure
plt.figure()
plt.plot(y, cl, y, cl_u, y, cl_n, y, cl_v, y, cl_uu)
plt.grid()
plt.xlabel('$y$ [m]')
plt.ylabel('$C_l$')
plt.xlim([-b/2, b/2])
plt.legend(['LL - with jet correction', 'LL - with  equivalent uniform jet correction', 'LL - no corr', 'LL - no jet', 'LL- uniform jet' ])
plt.savefig('Cl_a=%.1f_deg_N=%d.png' % (np.rad2deg(alpha), N), dpi=300)

plt.figure()
plt.plot(y, gamma, y, gamma_u, y, gamma_n, y, gamma_v)
plt.grid()
plt.xlabel('$y$ [m]')
plt.ylabel('$\Gamma$ [m2/s]')
plt.xlim([-b/2, b/2])
plt.legend(['LL - with jet correction', 'LL - with  equivalent uniform jet correction', 'LL - no corr', 'LL - no jet'])
plt.savefig('gamma_a=%.1f_deg_N=%d.png' % (np.rad2deg(alpha), N), dpi=300)

plt.figure()
plt.plot(y, V)
plt.grid()
plt.xlabel('$y$ [m]')
plt.ylabel('V [m/s]')
plt.xlim([-b/2, b/2])
plt.savefig('V_a=%.1f_deg_N=%d.png' % (np.rad2deg(alpha), N), dpi=300)


