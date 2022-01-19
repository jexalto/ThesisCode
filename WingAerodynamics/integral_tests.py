import numpy as np
from numpy.core.function_base import linspace
import scipy.special as sp
import time
from matplotlib import pyplot as plt

order = 1
dlambda = 0.1
lbda = np.arange(dlambda, 5-dlambda, dlambda)

# start = time.time()
# int_I = sp.iv(order, lbda)/lbda*dlambda
# end = time.time()

# startK = time.time()
# int_K = sp.kv(order, lbda)/lbda*dlambda
# endK = time.time()

c = 1.2
d = 1.3

start = time.time()
p = 0

dlbda = 0.1
lbda = np.arange(0, 5, dlbda)

dlb = 0.005
lb = np.arange(dlb, 1, len(lbda))

lb = c*lbda+lb*(d-c)*lbda
dlb = (d-c)*dlb*lbda

int_I = sp.iv(2*p+1, lb[5])/lb[5]*dlb[5]
int_I = np.sum(int_I, axis=0)
int_I = np.where(np.isnan(int_I), 0, int_I)

start = time.time()
it = 100
clambda = c * 0.6
dlambda = d * 0.6
d_lambda = (dlambda-clambda)/it
int_K = 0
for i in range(it):
    lbb = (clambda + (dlambda - clambda)*i/it)
    int_K += sp.kv(2*p+1, lbb)/lbb * d_lambda
end = time.time()

# int_K = np.sum(int_K, axis=0)
# int_K = np.where(np.isnan(int_K), 0, int_K)

# print("--- Python ---")
print("Kv integral: ", int_K)
print('Kv integral: time = %f seconds' % (end-start))