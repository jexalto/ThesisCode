import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import sys
import pandas as pd

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from vortices import v_induced_by_horseshoe_vortex

#%% INPUTS
# code for two jets at the wingtips

Vinf = 100              #velocity [m/s]
alpha = np.deg2rad(2)   #aoa [rad]

r0 = 1                  #jet radius [m]
mu = 2/3                #Vinf/Vjet [-]
crd = 1                 #wing chord [m]
A = 10                  #wing aspect ratio [m]
tc = 1.06               #wing thickness correction [-]

m = 10                   #number of panels in jet
n = 10                  #number of panels outside jet

Vj = Vinf/mu            #Vjet [m/s]
b = A*crd               #wingspan [m]
x = -crd*1/2            #control point x location [m]
dzeta = -crd*1/2/r0     #non-dim control point x [-]

#%% Calculate geometry

#left jet
yl_ = np.linspace(b/2-r0, b/2, m+2)
y_l = (yl_[1:]+yl_[:-1])/2
s_l = (yl_[1:]-yl_[:-1])/2

#inside jet
yp_ = np.linspace(-b/2+r0, b/2-r0, 2*n+1)
y_p = (yp_[1:]+yp_[:-1])/2
s_p = (yp_[1:]-yp_[:-1])/2

#right from jet
yr_ = np.linspace(-b/2, -b/2+r0, m+2)
y_r = (yr_[1:]+yr_[:-1])/2
s_r = (yr_[1:]-yr_[:-1])/2

#full wing
s = np.concatenate((s_r, s_p, s_l))
y = np.concatenate((y_r, y_p, y_l))
y_ = np.concatenate((yr_, yp_[1:-1], yl_))

#number of panels
N = 2*n+2*m+2

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
           
#%% Even vortex jet correction

#calculations are made in a system where the jet is on y=0
y2 = y+b/2-s[0]
y2_ = y_+b/2-s[0]
r0_ = r0-s[0]
#control point y coordinates
eta = y2/r0_
eta = eta.reshape(eta.size, 1)
#horseshoe vortex midpoitn y coordinates
beta = y2/r0_
#horseshoe vortex points
c = y2_[:-1]/r0_
d = y2_[1:]/r0_

#calculate correction
Gjje = 1/r0_*(1-mu**2)/(1+mu**2)*(1/(1/d-eta)-1/(1/c-eta)+1/(1/d+eta)-1/(1/c+eta))
Goje = -1/r0_*(1-mu)**2/(1+mu**2)*(1/(eta-c)-1/(eta-d)+1/(eta+d)-1/(eta+c))
Gjoe = Goje
Gooe = -Gjje

#determine if control point and horseshoe vortex are in jet
beta_in = np.where(beta<1, 1, 0)
eta_in = np.where(eta<1, 0.5, 0)
G_in = beta_in+eta_in

#apply right correction
Ge = np.where(G_in==0, Gooe, 0)
Ge = np.where(G_in==1, Goje, Ge)
Ge = np.where(G_in==0.5, Gjoe, Ge)
Ge = np.where(G_in==1.5, Gjje, Ge)
Ge[:,0] *= 0.5 #At beta=0 there's only one hs vortex

#clear memory
Gjje, Goje, Gjoe, Gooe = None, None, None, None

#%% Odd vortex jet correction

#change location of the horseshoe vortex on the symmetry plane
c[0] = 0
#to avoid singularities
eta = np.where(np.abs(eta)<1e-15, 1e-6, eta)

#create integration variables p, lbda, lb
p = np.arange(0, 15, 1)
p = p.reshape(p.size, 1, 1)

dlbda = 0.1
lbda = np.arange(0, 25, dlbda)
lbda = lbda.reshape(lbda.size, 1, 1, 1)

dlb = 0.002
lb = np.arange(dlb, 1, dlb)
lb = lb.reshape(lb.size, 1, 1, 1, 1)
lb = c*lbda+lb*(d-c)*lbda
dlb = (d-c)*dlb*lbda

#calculate integrals
int_I = sp.iv(2*p+1, lb)/lb*dlb
int_I = np.sum(int_I, axis=0)
int_I = np.where(np.isnan(int_I), 0, int_I)

int_K = sp.kv(2*p+1, lb)/lb*dlb
int_K = np.sum(int_K, axis=0)
int_K = np.where(np.isnan(int_K), 0, int_K)
lb = None

#calculate frequently used terms
Kp = sp.kvp(2*p+1, lbda)
K = sp.kv(2*p+1, lbda)
K_ = sp.kv(2*p+1, lbda*eta)
Ip = sp.ivp(2*p+1, lbda)
I = sp.iv(2*p+1, lbda)
I_ = sp.iv(2*p+1, lbda*eta)
sin = np.sin(dzeta*lbda)
mu1 = (1/((1/mu**2)-1))
mu2 = ((1/mu)-mu)

#calculate 4 corrections
int_jj = K*Kp*I_*sin/(mu1-lbda*I*Kp)*int_I*dlbda
int_jj = np.where(np.isnan(int_jj), 0, int_jj)
int_jj = np.sum(int_jj, axis=0)
sum_jj = (2*p+1)**2*int_jj
sum_jj = np.sum(sum_jj, axis=0)
Gjjo = 8/(r0_*np.pi*eta)*sum_jj
sum_jj = None

int_oj = (1/(mu-lbda*mu2*I*Kp)-1)*K_*sin/lbda*int_I*dlbda
int_oj = np.where(np.isnan(int_oj), 0, int_oj)
int_oj = np.sum(int_oj, axis=0)
sum_oj = (2*p+1)**2*int_oj
sum_oj = np.sum(sum_oj, axis=0)
Gojo = 8/(r0_*np.pi*eta)*sum_oj
sum_oj = None

int_jo = (1/(mu-lbda*mu2*I*Kp)-1)*I_*sin/lbda*int_K*dlbda
int_jo = np.where(np.isnan(int_jo), 0, int_jo)
int_jo = np.sum(int_jo, axis=0)
sum_jo = (2*p+1)**2*int_jo
sum_jo = np.sum(sum_jo, axis=0)
Gjoo = 8/(r0_*np.pi*eta)*sum_jo
sum_jo = None

int_oo = I*Ip*K_*sin/(mu1-lbda*I*Kp)*int_K*dlbda
int_oo = np.where(np.isnan(int_oo), 0, int_oo)
int_oo = np.sum(int_oo, axis=0)
sum_oo = (2*p+1)**2*int_oo
sum_oo = np.sum(sum_oo, axis=0)
Gooo = 8/(r0_*np.pi*eta)*sum_oo
sum_oo = None

#apply right correction
Go = np.where(G_in==0, Gooo, 0)
Go = np.where(G_in==1, Gojo, Go)
Go = np.where(G_in==0.5, Gjoo, Go)
Go = np.where(G_in==1.5, Gjjo, Go)

#Influence of single horseshoe vortex
G_ = Go+Ge

#%% Apply jet correction to wing and solve

G = np.zeros((N, N))

#mirror the influence matrix G to describe the correction on both sides
#of the jet
#G = G_+G_[::-1, ::-1]
G = G_[::-1, ::-1]

#calculate circulation with jet
#V = np.where(np.abs(y)>b/2-r0, Vj, Vinf)
V = np.where(y>b/2-r0, Vj, Vinf)
RHS = V*alpha
AIC = 1/(4*np.pi)*(F+G)
gamma = np.linalg.solve(AIC, RHS)
cl = 2*gamma*V/(Vinf**2*crd)

#calculate circulation with jet without rethorst correction
AIC_n = 1/(4*np.pi)*(F)
gamma_n = np.linalg.solve(AIC_n, RHS)
cl_n = 2*gamma_n*V/(Vinf**2*crd)

#calculate circulation without jet
AIC_v = 1/(4*np.pi)*(F)
RHS_v = Vinf*alpha*np.ones(RHS.shape)
gamma_v = np.linalg.solve(AIC_v, RHS_v)
cl_v = 2*gamma_v/(Vinf*crd)

#read cfd data
data_cfd = pd.read_csv('cfd3.txt')
y_cfd = data_cfd['z[m]'].to_numpy()
y_cfd *= -1
beta_cfd = y_cfd/r0
cl_cfd = data_cfd[' Cl'].to_numpy()
    
#create figure
plt.figure()
plt.plot(y, cl*tc, y, cl_n*tc, y, cl_v*tc, y_cfd, cl_cfd)
plt.grid()
plt.xlabel('$y$ [m]')
plt.ylabel('$C_l$')
plt.xlim([-b/2, b/2])
plt.legend(['LL - Rethorst', 'LL - no corr', 'LL - no jet', 'CFD - with jet'])
plt.savefig('Cl_a=%.1f_deg_N=%d.png' % (np.rad2deg(alpha), N), dpi=300)

plt.figure()
plt.plot(y, gamma, y, gamma_n, y, gamma_v)
plt.grid()
plt.xlabel('$y$ [m]')
plt.ylabel('$\Gamma$ [m2/s]')
plt.xlim([-b/2, b/2])
plt.legend(['LL - Rethorst', 'LL - no corr', 'LL - no jet'])
plt.savefig('Cl_a=%.1f_deg_N=%d_gamma.png' % (np.rad2deg(alpha), N), dpi=300)

