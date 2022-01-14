import numpy as np
import scipy.interpolate as si


class Inflow(object):
    """Contains the inflow flow field for the propeller, using the
    interpolation function, the flowfield can be evaluated at any point
    """
    
    def __init__(self):
        #non-dimensional velocities in x,y,z direction
        self.u = None
        self.v = None
        self.w = None
        
        #corresponding y and z coordinates
        self.y = None
        self.z = None
        
        #resolution of the flowfield when one axis is not defined
        self.res_flowfield = 0.01
        
    def reset(self):
        self.u, self.v, self.w = None, None, None
        self.y, self.z = None, None
        
    def line_input(self, y, z, u, v, w):
        """If the velocity field is only defined on one line in y or z 
        direction, this function will convert the input to a flattened grid
        which can be used by the interpolation function.
        
        y or z is given as a (N,) array and u,v,w are also (N,) arrays
        
        Args:
            y (numpy.array): None or y coordinates of line
            z (numpy.array): None or z coordinates of line
            u (numpy.array): velocity in x direction along line
            v (numpy.array): velocity in y direction along line
            w (numpy.array): velocity in z direction along line
        """
        
        self.reset()
        
        #check which axis is missing and set missing axis
        if y is None:
            y = np.arange(-1, 1+self.res_flowfield, self.res_flowfield)
            case = 'y'
        if z is None:
            z = np.arange(-1, 1+self.res_flowfield, self.res_flowfield)
            case = 'z'
        
        #dimensions of the grid
        m = y.size
        n = z.size
        
        #convert to grid
        if case=='y':
            u = np.ones((n, m))*u.reshape((n,1))
            v = np.ones((n, m))*v.reshape((n,1))
            w = np.ones((n, m))*w.reshape((n,1))
        elif case=='z':
            u = np.ones((n, m))*u
            v = np.ones((n, m))*v
            w = np.ones((n, m))*w
            
        y = np.ones((n, m))*y
        z = np.ones((n, m))*z.reshape((n,1))
        
        #flatten grids
        self.u = u.flatten()
        self.v = -v.flatten() #convert to left-handed coordinate system
        self.w = w.flatten()
        self.y = -y.flatten() #convert to left-handed coordinate system
        self.z = z.flatten()
        
        
    def grid_input(self, y, z, u, v, w):
        """Converts an input given as a grid to a flattened grid which
        can be used by the interpolation function.
        
        z is given as a (N,) array and y as a (M,) array. u,v,w are matrices
        with shape (N,M)
        
        Args:
            y (numpy.array): y coordinates of grid
            z (numpy.array): z coordinates of grid
            u (numpy.array): velocity in x direction
            v (numpy.array): velocity in y direction
            w (numpy.array): velocity in z direction
        """
        
        self.reset()
        
        #dimensions of grid
        m = y.size
        n = z.size
        
        #convert to grid
        y = np.ones((n, m))*y
        z = np.ones((n, m))*z.reshape((n,1))
        
        #flatten grids
        self.u = u.flatten()
        self.v = -v.flatten() #convert to left-handed coordinate system
        self.w = w.flatten()
        self.y = -y.flatten() #convert to left-handed coordinate system
        self.z = z.flatten()
        
    def interpolation(self, y, z):
        """Interpolates the defined flowfield to x and y coordinates
        
        Args:
            y, z (numpy.array): y and z coordinates with the same shape
            
        Returns:
            u_int, v_int, w_int (numpy.array): interpolated velocities
                corresponding to the y and z coorindates
        """
        
        u_int = si.griddata((self.y, self.z), self.u, (y, z),
                            method='linear',
                            fill_value=0)
        v_int = si.griddata((self.y, self.z), self.v, (y, z),
                            method='linear',
                            fill_value=0)
        w_int = si.griddata((self.y, self.z), self.w, (y, z),
                            method='linear',
                            fill_value=0)
        
        return u_int, v_int, w_int
        