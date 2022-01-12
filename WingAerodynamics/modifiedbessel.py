import numpy as np
import scipy.special as sp
import time

lst =  np.arange(0, 5, 0.1)

def modifiedbessel(order, x):
    bessel = 0

    for m in range(100):
        addition = (x/2)**(2*m + order)/(factorial(m) * factorial(m + order + 1 - 1))
        bessel += addition

    return bessel

def factorial(m):
    output = 1
    for i in range(1, m+1):
        output = output * i
    return float(output)


start = time.time()
mybessel = modifiedbessel(1, lst)
print('--- %s seconds passed: own ---' % round((time.time()-start), 6))

startscipy = time.time()
scipybessel = sp.iv(1, lst)
print('--- %s seconds passed: scipy ---' % round((time.time()-startscipy), 6))

print(' --- Convergence is: %s ---' % round(sum(abs(scipybessel - mybessel)), 15))

i=10
print(factorial(i))
print(sp.factorial(i))

print(scipybessel)
print(mybessel)