import numpy as np
import numpy as np


def actuator_disk(fax, ftgt, rho, Vinf, r, dr, R, omega):
    """ Computes axial and tangential induction factors for an incompressible
    acutator disk based on input axial and tangential loading.
    Original function written by Tomas Sinnige in MATLAB
    
    Args:
        fax (numpy.array): radial distr. of isolated thrust per annulus [N]
            vector of size: 1xlength(r)
        ftgt (numpy.array): radial distr. of isolated torque force per annulus [N]
            vector of size: 1xlength(r)
        rho (float): freestream air density [kg/m3]
        Vinf (float): freestream velocity [m/s]
        r (numpy.array): radial coordinates (with r varying along 2nd dim)
        dr (numpy.array): radial extent of blade segments 
        R_prop (float): propeller radius [m]
        omega (float): propeller rotational speed [rad/s]

    Returns:
        a (numpy.array): axial induction factor (a=dVa/Va) 
            vector of size: 1 x length(r)
        b (numpy.array): tangential induction factor (b=dVt/(omega*r)
            vector of size: 1 x length(r)
    """
    
    #compute axial induction factor - computation is done for both positive
    #and negative loading regimes and then results are merged
    apos = -0.5+0.5*(1+fax/(0.5*rho*Vinf**2*2*np.pi*r*dr))**0.5 #for positive loading
    aneg = -0.5+0.5*(1-fax/(0.5*rho*Vinf**2*2*np.pi*r*dr))**0.5 #for negative loading
    #set induction to positive or negative loading value according to fax
    a = np.where(fax<0, -aneg, apos) 
    #compute tangential induction factor
    b = ftgt/(4*np.pi*rho*Vinf*R*omega*r*dr*(1+a))
    
    return a,b
