import numpy as np
import scipy.special as sp
import time

start = time.time()
#create integration variables p, lbda, lb
mu = 0.9
dzeta = 0.5
eta = 0.3
r0_ = 1.0
c = 1.2
d = 1.3

p = np.arange(0, 10, 1)
p = p.reshape(p.size, 1, 1)

dlbda = 0.1
lbda = np.arange(0, 5, dlbda)
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

end = time.time()
print('--- PYTHON ---')
print('Gjj value: %f' %(Gjjo))
print('Goj value: %f' %(Gojo))
print('Gjo value: %f' %(Gjoo))
print('Gjo value: %f' %(Gooo))
print('Time needed to converge: %f seconds' %(end - start))