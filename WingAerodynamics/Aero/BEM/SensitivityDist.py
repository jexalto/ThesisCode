import numpy as np
import scipy.interpolate as si


class SensitivityDist(object):
    """A class that contains a propeller performance map (CT and CQ for varying
    J). To calculate the CT and CQ for any J, the interpolate function is used.
    """
    
    def __init__(self):
        
        self.J = None   #advance ratio
        self.rR = None  #r/R
        self.dCT = None #local thrust coefficient on grid defined by J and rR
        self.dCQ = None #local torque coefficient on grid defined by J and rR
        
    def interpolate(self, Ji, rRi):
        """Interpolates the propeller sensitivity distribution.
        
        Args:
            Ji (float/numpy.array): advance ratio(s)
                If it is a float, it is used for all radial stations
                If it is an array, the shape must correspond to rRi
            rRi (numpy.array): radial stations
            
        Returns:
            dCTi, dCQi (numpy.array): interpolated thrust and torque values
        """
        
        #find out if Ji is a float or array
        if type(Ji)==float:
            form = '2D'
        else:
            if Ji.shape[1]!=rRi.size:
                raise Exception('Input error')
            else:
                form = '3D'
         
        Jout = False
        #check if the value(s)  of J are within the range, otherwise replace
        #with the minimum or maximum
        if form=='3D':
            if np.max(Ji)>np.max(self.J) or np.min(Ji)<np.min(self.J):
                Jout = True
                Ji = np.where(Ji>np.max(self.J),
                              np.max(self.J),
                              Ji)
                
                Ji = np.where(Ji<np.min(self.J),
                              np.min(self.J),
                              Ji)
        else:
            if Ji>np.max(self.J):
                Jout = True
                Ji = np.max(self.J)
            elif Ji<np.min(self.J):
                Jout = True
                Ji = np.min(self.J)
        
        #print warning if value(s) in Ji are out of range
        if Jout:
            s = 'Advance ratio out of interpolation zone, J is limited!'
            print(s)
        
        #matrix size
        m = self.J.size
        n = rRi.size
        
        if form=='3D':
            #empty grids
            dCT_ = np.zeros((m,n))
            dCQ_ = np.zeros((m,n))
            
            #interpolate on the radial values
            for i in range(m):
                dCT_[i, :] = si.pchip(self.rR, self.dCT[i, :])(rRi)
                dCQ_[i, :] = si.pchip(self.rR, self.dCQ[i, :])(rRi)
            
            #grid to save results
            dCTi = np.zeros(Ji.shape)
            dCQi = np.zeros(Ji.shape)
            
            #interpolate to J for each radius
            for i in range(n):
                dCTi[:, i] = si.pchip(self.J, dCT_[:, i])(Ji[:, i])
                dCQi[:, i] = si.pchip(self.J, dCQ_[:, i])(Ji[:, i])
            
        else:
            #empty arrays
            dCT_ = np.zeros(self.rR.size)
            dCQ_ = np.zeros(self.rR.size)
            
            #interpolate to J for earch radius
            for i in range(self.rR.size):
                dCT_[i] = si.pchip(self.J, self.dCT[:, i])(Ji)
                dCQ_[i] = si.pchip(self.J, self.dCQ[:, i])(Ji)
                
            #interpolate on the radial values
            dCTi = si.pchip(self.rR, dCT_)(rRi)
            dCQi = si.pchip(self.rR, dCQ_)(rRi)
        
        return dCTi, dCQi
        