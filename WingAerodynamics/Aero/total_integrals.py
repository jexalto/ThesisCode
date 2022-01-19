import numpy as np
import scipy.special as sp
import time

start  = time.time()
c = 1.2
d = 1.3
mu = 0.9
dzeta = 0.5
eta = 0.1
r0 = 1.0

p = np.arange(0, 10, 1)
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

input  = 0.5
order = 2
# print('Iv with input %f and order %f is %f' %(input, order, sp.iv(order, input)))

#calculate 4 corrections
int_jj = K*Kp*I_*sin/(mu1-lbda*I*Kp)*int_I*dlbda
int_jj = np.where(np.isnan(int_jj), 0, int_jj)
int_jj = np.sum(int_jj, axis=0)

sum_jj = (2*p+1)**2*int_jj

sum_jj = np.sum(sum_jj, axis=0)
Gjjo = 8/(r0*np.pi*eta)*sum_jj
end = time.time()

print('--- PYTHON ---')
print('Gjj value: %f' %(Gjjo))
print('Time needed to converge: %f seconds' %(end - start))