import matplotlib
import numpy as np
import scipy.special as sp
import time
from matplotlib import pyplot as plt

lst =  np.arange(1, 10, 0.1)
order = 4

def modifiedbessel(order, x):
    bessel = 0
    pi = 3.141592653589793
    t = 0
    steps=1000
    dt = 2*pi/steps

    for m in range(10):
        addition = (x/2)**(2*m + order)/(factorial(m) * factorial(m + order))
        bessel += addition
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
plt.show()
# alpha = 0.5
# bessel_der = 0.5 * (modifiedbessel(1-1, alpha) + modifiedbessel(1+1, alpha))
print('Second kind bessel function: ', secondorderbessel_cosh(3, 4.8))
# print('Scipy bessel function: ', sp.ivp(2, alpha))