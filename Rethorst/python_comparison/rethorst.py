from jet_correction_ib import jet_correction_ib
import numpy as np

ny = 20
nx = 5

y = np.linspace(0, 20, ny)
mu = 1.2
crd = 4.0
loc = 9.5
ib_loc = 10.0

G = jet_correction_ib(y, mu, crd, loc, ib_loc, dlbda=0.1, dlb=0.005, p_max=10, lbda_max=5)

print(G)

# --- Durham's method ---

Ct = 0.4
a = (np.sqrt(1 + 8/np.pi * Ct) - 1)/2

def Vx(a, V, Rp, x):
    Vx = a * V * (1 + x/(Rp*Rp + x*x))
    return Vx

def Vy(a, V, Rp, x, y):
    Vy = 0.5 * a * V * (1 + (Rp*Rp*y)/(Rp*Rp + x*x)*np.sqrt(Rp*Rp + x*x))
    return Vy

def Vz(a, V, Rp, x, z):
    Vz = 0.5 * a * V * (1 + (Rp*Rp*z)/(Rp*Rp + x*x)*np.sqrt(Rp*Rp + x*x))
    return Vz

def R1(Rp, a, x):
    R1_ = Rp * np.sqrt( (1+a)/( 1+a*( 1+x/np.sqrt(Rp*Rp + x*x) ) ) )
    return R1_

for i in range(10):
    x = i/10
    Rp = 1
    V = 40
    print("Vx %f" % Vx(a, V, Rp, x))
    print("Vy %f" % Vy(a, V, Rp, 1, x))
    print("Vz %f" % Vz(a, V, Rp, 1, x))