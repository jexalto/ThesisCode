from math import cos, sin, atan
import numpy as np
import os
import scipy.interpolate as si

from Xfoil import Xfoil
from profile_drag import ProfileDrag


class NACA4(object):
    """
    4-digit NACA airfoils generator. Calculates the camber line,
    its gradient, (half) thickness distribution and the position
    of both lower and upper surfaces. Default values are set for
    the NACA 2412.

    Reference:
    [1] http://airfoiltools.com/airfoil/naca4digit
    [2] NACA Report 824, pp 262

    Modified from: https://github.com/aqreed/PyVLM
    """

    def __init__(self, name):
        """Calculate the NACA4 airfoil coefficients
        
        Args:
            name (str): String with four numbers
        """
        
        self.M = int(name[0])/100
        self.P = int(name[1])/10
        self.T = int(name[2:])/100
        self.name = 'NACA'+name

        if self.M < 0 or self.M > 9.5/100:
            msg = 'Max camber M should be within [0, 9.5]'
            raise ValueError(msg)

        if self.P < 0 or self.P > 9/10:
            msg = 'Max camber position P should be within [0, 9]'
            raise ValueError(msg)

        if self.T < 0 or self.T > 40/100:
            msg = 'Thickness T should be within [0, 40]'
            raise ValueError(msg)
            
        self.polar_inviscid = None
        self.polar_laminar = None
        self.polar_turbulent = None

    def camber_line(self, x):
        """Returns the (local) camber line of the airfoil
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            z (float): camber at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        M = self.M
        P = self.P

        if x < P:
            z = (M/P**2) * x * (2*P - x)
        else:
            z = (M / (1 - P)**2) * (1 - 2*P + x * (2*P - x))

        return z

    def camber_gradient(self, x):
        """Returns the (local) camber gradient of the airfoil
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            dz (float): camber gradient at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        M = self.M
        P = self.P
        
        if type(x)==float:
            if x < P:
                dz = (2*M/P**2) * (P - x)
            else:
                dz = (2*M / (1 - P)**2) * (P - x)
        else:
            if P!=0:
                dz1 = (2*M/P**2) * (P - x)
            else:
                dz1 = np.zeros(x.shape)
            dz2 = (2*M / (1 - P)**2) * (P - x)
            dz = np.where(x<P, dz1, dz2)

        return dz

    def thickness(self, x):
        """Returns the (local) half-thickness distribution.
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            t (float): thickness at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        T = self.T

        # Values for a t=20% airfoil
        a0, a1, a2, a3, a4 = 0.2969, -0.126, -0.3516, 0.2843, -0.1015

        # Half thickness (corrected from t=20% values)
        t = (T/0.2) * (a0*x**0.5 + x*(a1 + x*(a2 + x*(a3 + x*a4))))

        return t

    def upper_surface(self, x):
        """Returns the position of the upper surface
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            xu, yu (float): corresponding x and y coordiante
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        theta = atan(self.camber_gradient(x))

        xu = x - self.thickness(x) * sin(theta)
        yu = self.camber_line(x) + self.thickness(x) * cos(theta)

        return xu, yu

    def lower_surface(self, x):
        """Returns the position of the lower surface.
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            xl, yl (float): corresponding x and y coordiante
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        theta = atan(self.camber_gradient(x))

        xl = x + self.thickness(x) * sin(theta)
        yl = self.camber_line(x) - self.thickness(x) * cos(theta)

        return xl, yl
    
    def thickness_correction(self):
        """Returns a thickness correction based on thin airfoil theory and
        linear-vorticity stream function panel method by XFOIL
        
        Retruns:
            corr (float): correction value
        """
        
        #create xfoil instance
        xf = Xfoil()
        xf.NACA = self.name[4:]
        xf.visc = False
        xf.run_dir = 'temp/'
        xf.xfoil_dir = 'data_xfoil/'
        
        #run xfoil
        cl1,__,__,__ = xf.run_xfoil(-2.5)
        cl2,__,__,__ = xf.run_xfoil(2.5)
        
        #find linear fit
        p = np.polyfit([-2.5, 2.5], [cl1, cl2], 1)
        cla = np.rad2deg(p[0])
        
        #base correction on slope
        corr = cla/(2*np.pi)
        
        return corr
    
    def init_polars(self, polar_dir):
        """Creates polar instances based on precalculated Xfoil data.
        
        Args:
            polar_dir (str): main directory where Xfoil polars are stored,
                the directory must contain folders called 'inviscid', 'laminar'
                and 'turbulent'
        """
        
#        self.polar_inviscid = ProfileDrag(self.name,
#                                          polar_dir+'/inviscid/',
#                                          False)
        #self.polar_laminar = ProfileDrag(self.name,
        #                                 polar_dir+'/laminar/')            # TURNED OFF
        self.polar_turbulent = ProfileDrag(self.name,
                                           polar_dir+'/turbulent/')
        
    def cl_thin_airfoil(self, alpha):
        """Calculates the lift coefficient using thin airfoil theory
        
        Args:
            alpha (float): angle of attack in degrees
        
        Returns:
            cl (float): lift coefficient
        """
        
        theta = np.linspace(0, np.pi, 100)
        dzdx = self.camber_gradient((1-np.cos(theta))/2)
        a0 = -1/np.pi*np.trapz(dzdx*(np.cos(theta)-1), theta)
        cl = 2*np.pi*(np.deg2rad(alpha)-a0)
        
        return cl
        

class Airfoil(object):
    """Contains all airfoil data. Airfoil coordinates can be given to an
    instance of the class by using the input_selig function. Based on this
    airfoil properties can be calculated.
    """
                    
    def input_selig(self, x, y, name):
        """Converts the airfoil coordinates in Selig format (with chordwise
        coordiantes from 1-->0-->1, with corresponding coordinates of first
        the upper surface followed by the lower surface). The number of x and y
        coordinates must be odd, with the leading edge coordinate being in the
        middle of the arrays.
        
        Args:
            x (numpy.array): chordwise values (1-->0-->1)
            y (numpy.array): upper surface values followed by lower surface
                values corresponding to x
            name (str): airfoil name
        """
        
        #find the array location of the airfoil leading edge
        N = np.size(x)
        n = np.where(x==0)[0][0] +1#int(N/2)+1
        
        #save inputs
        self.name = name
        self.x_s = x
        self.y_s = y
        
        #put upper surface coordinates in order of increasing x
        self.x = np.sort(x[:n])
        self.upper = y[:n][np.argsort(x[:n])]
        #interpolate lower surface coordinates, so x match for both upper and
        #lower surfaces
        lw = si.pchip(x[n-1:], y[n-1:])
        self.lower = lw(self.x)
        
        #camber line
        self.z = (self.upper+self.lower)/2
        #thickness
        self.t = self.upper-self.lower
        
        #create interpolation functions
        self.f_z =  si.pchip(self.x, self.z)
        self.f_dz = self.f_z.derivative(1)
        self.f_t =  si.pchip(self.x, self.t)
        
        self.polar_inviscid = None
        self.polar_laminar = None
        self.polar_turbulent = None
        
    def camber_line(self, x):
        """Returns the (local) camber line of the airfoil
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            z (float): camber at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)
        
        z = self.f_z(x)

        return z

    def camber_gradient(self, x):
        """Returns the (local) camber gradient of the airfoil
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            dz (float): camber gradient at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        dz = self.f_dz(x)

        return dz

    def thickness(self, x):
        """Returns the (local) half-thickness distribution.
        
        Args:
            x (float): chordwise position between 0 and 1
            
        Returns:
            t (float): thickness at x
        """

        if np.any(x < 0) or np.any(x > 1):
            msg = 'Argument x should be within [0, 1]'
            raise ValueError(msg)

        t = self.f_t(x)

        return t   
    
    def thickness_correction(self):
        """Returns a thickness correction based on thin airfoil theory and
        linear-vorticity stream function panel method by XFOIL
        
        Retruns:
            corr (float): correction value
        """
        
        #directory for analysis
        if not os.path.exists('temp/'):
            os.mkdir('temp/')
            
        #write coordinate file
        f = open('temp/%s.dat' % self.name, 'w')
        f.write('%s\n' % self.name)
        for i in range(self.x_s.size):
            f.write('%.6f %.6f\n' % (self.x_s[i], self.y_s[i]))
        f.close()
        
        #create xfoil instance
        xf = Xfoil()
        xf.airfoil = self.name
        xf.visc = False
        xf.run_dir = 'temp/'
        xf.xfoil_dir = 'data_xfoil/'
        
        #run xfoil
        cl1,__,__,__ = xf.run_xfoil(-2.5)
        cl2,__,__,__ = xf.run_xfoil(2.5)
        
        #find linear fit
        p = np.polyfit([-2.5, 2.5], [cl1, cl2], 1)
        cla = np.rad2deg(p[0])
        
        #base correction on slope
        corr = cla/(2*np.pi)
        
        return corr
    
    def init_polars(self, polar_dir):
        """Creates polar instances based on precalculated Xfoil data.
        
        Args:
            polar_dir (str): main directory where Xfoil polars are stored,
                the directory must contain folders called 'inviscid', 'laminar'
                and 'turbulent'
        """
        
#        self.polar_inviscid = ProfileDrag(self.name,
#                                          polar_dir+'/inviscid/',
#                                          False)
#        self.polar_laminar = ProfileDrag(self.name,
#                                         polar_dir+'/laminar/')                    # TURNED OFF
        self.polar_turbulent = ProfileDrag(self.name,
                                           polar_dir+'/turbulent/')
    
    def cl_thin_airfoil(self, alpha, M):
        """Calculates the lift coefficient using thin airfoil theory
        
        Args:
            alpha (float): angle of attack in degrees
            M (float): MAch number used to compressibility correction
        
        Returns:
            cl (float): lift coefficient
        """
        
        theta = np.linspace(0, np.pi, 100)
        dzdx = self.f_dz((1-np.cos(theta))/2)
        a0 = -1/np.pi*np.trapz(dzdx*(np.cos(theta)-1), theta)
        cl = 2*np.pi*(np.deg2rad(alpha)-a0)
        cl = cl*1/np.sqrt(1-M**2)               #compressibiliy correction
        
        return cl
    