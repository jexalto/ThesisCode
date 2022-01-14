import numpy as np
import scipy.special as sp


def sears_function(sigma, M):
    """Computes the value(s) of the Sears' function S, including a
    compressibility correction.
    Original function written by Tomas Sinnige in MATLAB
    
    Args:
        sigma (float/numpy.array): reduced frequency (can be a vector)
        M (float) Mach number

    Returns:
        S_lf (numpy.array): containing values of compressible low-frequency
            Sears' function
    """
    
    #Compute Prandtl-Glauert compressibility factor
    beta = np.where(M<0.7, (1-M**2)**0.5, (1-0.7**2)**0.5)
    #Compute value incompressible Sears' function
    #compute incompressible result using subfunction SearsFunctionM0
    x = sigma/beta**2
    S_M0 = sears_function_M0(x)
    #correct for erroneous values in case x=0
    S_M0 = np.where(x==0, 1, S_M0) #Sears' function should be 1 at x=0
    #Apply compressibility correction
    #compute value compressibility correction
    f = (1-beta)*np.log(M)+beta*np.log(1+beta)-np.log(2)
    Mcorrection = (sp.jv(0, M**2*x) + 1j*sp.jv(1, M**2*x))
    #compute compressible Sears' function 
    S_lf = S_M0*Mcorrection*np.exp(-1j*x*f)
    
    return S_M0, S_lf

def sears_function_M0(x):
    """Computes the value(s) of the incompressible Sears' function S(x) for the 
    given input value(s) x.
    Original function written by Tomas Sinnige in MATLAB
    
    Args:
        x (numpy.array): input vector containing values for which the
            incompressible Sears' function is computed
            
    Returns:
        S_M0 (numpy.array): containing values of incompressible Sears' function
            for all value(s) of input vector x
    """
    
    #compute incompressible Sears' function
    S_M0 = (1j*x*(sp.kv(0, 1j*x)+sp.kv(1, 1j*x)))**(-1)
    #correct for erroneous values in case x=0
    S_M0 = np.where(x==0, 1, S_M0) #Sears' function should be 1 at x=0

    return S_M0
