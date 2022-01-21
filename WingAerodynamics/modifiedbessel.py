import matplotlib
import numpy as np
import scipy.special as sp
import time
from matplotlib import pyplot as plt

<<<<<<< HEAD
lst =  np.arange(1, 10, 0.1)
order = 4
=======
lst =  np.arange(0, 5.1, 0.1) # np.array([0.1])
>>>>>>> 66d153c2425ddf4f7062e8084c690c655e326a7e

def modifiedbessel(order, x):
    bessel = 0
    pi = 3.141592653589793
    t = 0
    steps=1000
    dt = 2*pi/steps

<<<<<<< HEAD
    for m in range(5):
=======
    for m in range(100):
        # Rewrite factorial into smarter call
>>>>>>> 66d153c2425ddf4f7062e8084c690c655e326a7e
        addition = (x/2)**(2*m + order)/(factorial(m) * factorial(m + order))
        bessel += addition
        print(m, m+order, factorial(m), factorial(m + order))
    # for m in range(steps):
    #     addition = 1/pi * np.exp(x*np.cos(t)) * np.cos(order*t) * dt
    #     bessel += addition
    #     t = t + dt

    return bessel

def factorial(m):
    output = 1
    for i in range(1, m+1):
        output = output * i
    
    return float(output)

def gamma(z):
    discretisation = 100
    dx = 10/discretisation
    x = 0.001
    gamma = 0

    for i in range(discretisation):
        gamma += x**(z-1)*np.exp(-x)*dx
        x += dx

    return gamma

def secondorderbessel(order, x):
    sk_bessel = 0
    t = 0.001
    dt = 0.01
    for i in range(1000):
        sk_bessel += gamma(order+0.5)*(2*x)**order/np.sqrt(np.pi) * np.cos((t))/(t**2 + x**2)**(order+0.5)*dt
        t += dt
    
    return sk_bessel

<<<<<<< HEAD
def secondorderbessel_cosh(order, x):
    sk_bessel = 0
    t = 0.0001
    dt = 0.01
    for i in range(1000):
        sk_bessel += np.exp(-x*np.cosh(t)) * np.cosh(order*t) * dt
        t += dt
    
    return sk_bessel

plt.plot(lst, (modifiedbessel(order, lst)-sp.iv(order, lst))/sp.iv(order, lst), label='hometeam')
# plt.plot(lst, sp.kv(4, lst), label='competition')
plt.legend()
plt.grid()
print('-----')
# print(sp.iv(5, 4.8))
print(modifiedbessel(5, 4.8))
# plt.show()
# alpha = 0.5
# bessel_der = 0.5 * (modifiedbessel(1-1, alpha) + modifiedbessel(1+1, alpha))
# print('Second kind bessel function: ', secondorderbessel_cosh(3, 4.8))
# print('Scipy bessel function: ', sp.ivp(2, alpha))
=======
print(scipybessel)
print(sum(mybessel))
>>>>>>> 66d153c2425ddf4f7062e8084c690c655e326a7e
