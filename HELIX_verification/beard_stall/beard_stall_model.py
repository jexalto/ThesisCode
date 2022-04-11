from helix_pyf90.mod_stall_beard import beard_stall
import numpy as np

def beard(cl, alpha, alpha_0, M):
    cl_stall = 0.
    for iAlpha, Alpha in alpha:
        CL
        beard_stall(cl_stall, cl, Alpha, alpha_0, M)
        print(cl_stall)
    return cl_stall

CL_alpha = 2*np.pi
alpha_0=7.
alpha=np.linspace(-8, 12, 10)
cl=alpha*CL_alpha/360
M=50.

_ = beard(0., 9, alpha_0, M)