import numpy as np
import scipy.special as sp
import time

start  = time.time()

mu = 0.9
dzeta = 0.5
eta = 0.3
r0 = 1.0
c = 1.2
d = 1.3

p = np.arange(0, 4, 1)
p = p.reshape(p.size, 1, 1)

dlbda = 0.1
lbda = np.arange(dlbda, 5, dlbda)
lbda = lbda.reshape(lbda.size, 1, 1, 1)

dlb = 0.005
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

int_oj = (1/(mu-lbda*mu2*I*Kp)-1)*K_*sin/lbda*int_I*dlbda
int_oj = np.where(np.isnan(int_oj), 0, int_oj)
int_oj = np.sum(int_oj, axis=0)
print(int_oj[0:-1, 0, 0])
sum_oj = (2*p+1)**2*int_oj
sum_oj = np.sum(sum_oj, axis=0)
Gojo = 8/(r0*np.pi*eta)*sum_oj

end = time.time()

print('--- PYTHON ---')
print('Gjj value: %f' %(Gojo))
print('Time needed to converge: %f seconds' %(end - start))