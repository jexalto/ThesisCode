import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so
import scipy.special as sp


class SlipstreamConway(object):
    """This class uses an analytical solution of a slipstreamtube model
    to quickly calculate the propeller slipstream induced velocities.
    The method consists of the superposition of solutions of even polynomials,
    given by Conway. This solution does not gives only axial and radial
    velocities.
    
    J.T. Conway (1995) Analytical solutions for the actuator disk with variable
    radial distribution of load
    """
    
    def __init__(self):
        
        self.Vz0 = None         #scaling factor used by Conway
        self.ds = 0.1           #integral discretization setting
        self.s_max = 500        #integral upper limit
        
    def calculate_Vz0(self, Vax, rR, plot=False):
        """Matches a solution of eight even polynomials to the axial velocity
        using a least squares approach.
        
        Args:
            Vax (numpy.array): axial velocities at the propeller plane
            rR (numpy.array): r/R of radial stations 
            plot (bool): shows plot of the polynomial fit
        """
        
        #exponents of polynomials
        mu = np.arange(1, 9, 1)
        mu = mu.reshape(mu.size, 1)
        
        #least squares
        X = (1-rR**2)**mu
        X = np.transpose(X)
        
        Y = Vax.reshape(Vax.size, 1)
        
        X_X1 = np.linalg.inv(np.dot(np.transpose(X), X))
        Vz0 = np.dot(np.dot(X_X1, np.transpose(X)), Y)
        self.Vz0 = Vz0.flatten()
        
        #plot
        if plot:
            Vz0 = self.Vz0
            Vz0 = Vz0.reshape(Vz0.size, 1)
            mu = np.arange(1, Vz0.size+1, 1)
            mu = mu.reshape(mu.size, 1)
            
            V0 = Vz0*(1-(rR)**2)**mu
            V0 = np.sum(V0, axis=0)
            
            plt.figure()
            plt.plot(rR, Vax, rR, V0)
            plt.grid()
            plt.xlabel('r/R [-]')
            plt.ylabel('$V_{ax}$ [m/s]')
            plt.legend(['Input', 'Conway polynomial'])
            plt.savefig('conway_comparison.png',
                        dpi=300)
            
    def obj(self, x, Vax, rR):
        """Calcultes the axial velocity distribution by superposition of even
        polynomials and the error with the real distribution.
        
        Args:
            x (numpy.array): design vector, contains values for Vz0
            Vax (numpy.array): axial velocities at the propeller plane
            rR (numpy.array): r/R of radial stations 
                
        Returns:
            err (float): error between polynomial solution and the actual values
        """
        
        #resize
        x = x.reshape(x.size, 1)
        mu = np.arange(1, x.size+1, 1)
        mu = mu.reshape(mu.size, 1)
        
        #calcualte Vax based on polynomials
        Vax2 = x*(1-rR**2)**mu
        
        #error calculation
        err = np.abs(Vax-np.sum(Vax2, axis=0))
        err = np.sum(err)
        #normalize
        err = err/np.abs(np.mean(Vax))
        
        return err
    
    def v_axial(self, x, rR, R):
        """Calculates the axial velocity in the propeller slipstream (x>0).
        
        Args:
            x (numpy.array): x locations, shape (N,)
            rR (numpy.array): r/R locations, shape (M,)
            R (float): propeller radius
        
        Returns:
            Vax (numpy.array): axial velocities, shape (N,M)
        """
        
        #reshape
        rR = rR.reshape(rR.size, 1)
        
        Vz0 = self.Vz0
        Vz0 = Vz0.reshape(Vz0.size, 1, 1)
    
        mu = np.arange(1, Vz0.size+1, 1)
        mu = mu.reshape(mu.size, 1, 1)
        
        s = np.arange(self.ds, self.s_max, self.ds)
        s = s.reshape(s.size, 1 ,1, 1)
        
        #calculate coefficients
        A1 = np.exp(-s*np.abs(x))*sp.jv(mu+1, s*R)*sp.jv(0, s*rR*R)
        A2 = A1/s**mu
        A3 = -2**mu*sp.gamma(mu+1)/R**(mu-1)
        A4 = (1-(rR)**2)**mu
        A5 = 2*A4+A3*np.trapz(A2, s, axis=0)
        
        #axial velocity
        Vax = np.sum(Vz0*A5, axis=0)
        
        return Vax

    def v_radial(self, x, rR, R):
        """Calculates the radial velocity in the propeller slipstream (x>0).
        
        Args:
            x (numpy.array): x locations, shape (N,)
            rR (numpy.array): r/R locations, shape (M,)
            R (float): propeller radius
        
        Returns:
            Vr (numpy.array): radial velocities, shape (N,M)
        """
        
        #reshape
        rR = rR.reshape(rR.size, 1)
        
        Vz0 = self.Vz0
        Vz0 = Vz0.reshape(Vz0.size, 1, 1)
    
        mu = np.arange(1, Vz0.size+1, 1)
        mu = mu.reshape(mu.size, 1 ,1)
        
        s = np.arange(self.ds, self.s_max, self.ds)
        s = s.reshape(s.size, 1 ,1, 1)
        
        #calculate coefficients
        A1 = np.exp(-s*np.abs(x))*sp.jv(mu+1, s*R)*sp.jv(1, s*rR*R)
        A2 = A1/s**mu
        A3 = -2**mu*sp.gamma(mu+1)/R**(mu-1)
        A4 = A3*np.trapz(A2, s, axis=0)
        
        #radial velocity
        Vr = np.sum(Vz0*A4, axis=0)
        
        return Vr
    