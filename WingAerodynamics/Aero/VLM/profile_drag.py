import os
import numpy as np
import pandas as pd
import scipy.interpolate as si


class ProfileDrag(object):
    """Class that contains polar data from Xfoil analysis. Functions cl and
    profile drag can be used to find values of lift and drag for a certain
    alpha, M and Re.
    """
    
    def __init__(self, airfoil_name, airfoil_folder, visc=True):
        
        self.airfoil_name = airfoil_name        #corresponds to the file names
        self.airfoil_folder = airfoil_folder    #folder with polar files
        self.visc = visc                        #boolean
        
        self.init_interp_func()
    
    def init_interp_func(self):
        """Reads the data files and creates an interpolation function"""
        
        if self.visc:
            f_cl, M_lst, alpha_lst, Re_lst = self.interp_func('cl')
            f_cd,__,__,__ = self.interp_func('cd')
            self.Re_lst = Re_lst
        else:
            f_cl, M_lst, alpha_lst = self.interp_func_inviscid('cl')
            f_cd,__,__ = self.interp_func_inviscid('cd')
        
        #save
        self.f_cl = f_cl
        self.f_cd = f_cd
        self.M_lst = M_lst
        self.alpha_lst = alpha_lst
        
    def cl(self, M, alpha, Re=0):
        """Returns the lift coefficient
        
        Args:
            M (float): Mach number
            alpha (float): angle of attack [deg]
            Re (float): Reynolds number
            
        Returns:
            cl[0] (float): interpolated lift coefficient at input conditions
        """
        
        M_lst = self.M_lst
        alpha_lst = self.alpha_lst
        
        f_cl = self.f_cl
        
        #check limits of the interpolation function
        if self.visc:
            Re_lst = self.Re_lst
            Re = self.check_limit(Re_lst, Re)
        M = self.check_limit(M_lst, M)
        alpha = self.check_limit(alpha_lst, alpha)
        
        #interpolate
        if self.visc:
            cl = f_cl(np.array([M, alpha, Re]))
        else:
            cl = f_cl(np.array([M, alpha]))
        
        return cl[0]
    
    def profile_drag(self, M, cl, alpha, Re=0):
        """Returns the profile drag coefficient based on lift coefficient. If
        multiple matching drag values are found for one lift coefficient, the
        angle of attack will be used to determine the drag value.
        
        Args:
            M (float): Mach number
            cl (float): lift coefficient
            alpha (float): angle of attack [deg]
            Re (float): Reynolds number
            
        Returns:
            cd[0] (float): interpolated profile drag coefficient at input conditions
        """
        
        M_lst = self.M_lst
        alpha_lst = self.alpha_lst
        
        f_cl = self.f_cl
        f_cd = self.f_cd
        
        #check limits of the interpolation function
        if self.visc:
            Re_lst = self.Re_lst
            Re = self.check_limit(Re_lst, Re)
        M = self.check_limit(M_lst, M)
        alpha = self.check_limit(alpha_lst, alpha)
        
        #get range of lift coefficients at M and Re
        if self.visc:
            pts = np.zeros([len(alpha_lst), 3])
            pts[:,0] = M
            pts[:,1] = alpha_lst
            pts[:,2] = Re
        else:
            pts = np.zeros([len(alpha_lst), 2])
            pts[:,0] = M
            pts[:,1] = alpha_lst
        
        cl_lst = f_cl(pts)
        
        #check if cl is in limits
        cl = self.check_limit(cl_lst, cl)
        
        #find machting lift coefficient values
        match = []
        for i in range(len(cl_lst)-1):
            cl_2 = cl_lst[i:i+2]
            if cl<max(cl_2) and cl>min(cl_2):
                match.append(i)
        
        #determine the which combination of cl and angle of attack is closest
        #to the input cl and angle of attack
        if len(match)>1:
            alpha_match = alpha_lst[match]
            ind = match[np.argmin(np.abs(alpha_match-alpha))]
        else:
            ind = match[0]
        
        #interpolate angle of attack
        da_dcl = (alpha_lst[ind+1]-alpha_lst[ind])/(cl_lst[ind+1]-cl_lst[ind])
        alpha2 = alpha_lst[ind]+da_dcl*(cl-cl_lst[ind])
        
        #calculate cd using interpolated angle of attack
        if self.visc:
            cd = f_cd(np.array([M, alpha2, Re]))
        else:
            cd = f_cd(np.array([M, alpha2]))            
        
        return cd[0]
    
    def check_limit(self, lst, val):
        """Returns a value within a interval. If the value is outside the
        interval, the minimum or maximum will be returned.
        
        Args:
            lst (list/numpy.array): list of possible values
            val (float): value to be checked
        
        Returns:
            val (float): input value within the boundries
        """

        if val>max(lst):
            print('value' + str(val) + ' is too high')
            val = max(lst)-0.001*abs(max(lst))
        elif val<min(lst):
            print('value' + str(val) + ' is too low')
            val = min(lst)+0.001*abs(min(lst))
        
        return val
    
    def interp_func(self, var):  
        """Reads the polars from the folder and constructs a interpolation
        function using the data.
        
        Args:
            var (str): variable for which the interpolation function is constructed
            
        Returns:
            f (scipy.interpolate.RegularGridInterpolator): the interpolation function
            M_lst (numpy.array): list with available Mach numbers
            alpha_lst (numpy.array): list of available angles of attack
            Re_lst (numpy.array): list of available Reynolds numbers
        """          
        
        airfoil_folder = self.airfoil_folder
        airfoil_name = self.airfoil_name
        #print(airfoil_folder, airfoil_name)
        
        #get file names
        dir_lst = os.listdir(airfoil_folder)
        dir_lst = [i for i in dir_lst if i[:len(airfoil_name)]==airfoil_name]
        
        Re_lst = []
        M_lst = []
        
        #determine the available values of Mach and Reynolds numbers
        for d in dir_lst:
            d = d[:-len(d.split('.')[-1])-1]
            spl = d.split('_')
            for s in spl:
                if s[:2] == 'Re':
                    Re_lst.append(float(s[2:]))
                if s[0] == 'M':
                    M_lst.append(float(s[1:]))
            
        Re_lst = np.unique(Re_lst)
        M_lst = np.unique(M_lst)
        
        #read the data
        for i,M in enumerate(M_lst):
            Re_data = None
            for j,Re in enumerate(Re_lst):
                #file name identifier
                fname = '_'.join([airfoil_name, 'Re%d' % Re, 'M%0.2f' % M])
                missing = True
                for d in dir_lst:
                    if fname in d:
                        #read data from file
                        if airfoil_folder is None:
                            df = pd.read_csv(d, header=6)
                        else:
                            df = pd.read_csv(airfoil_folder+'/'+d, header=6)
                        df = df.set_index('alpha')
                        missing = False
                
                if missing:
                    print('File not found')
                    print(fname)
                else:
                    #put data in a pandas DataFrame
                    if Re_data is None:
                        Re_data = df[var].rename(Re)
                    else:
                        Re_data = pd.merge(Re_data, df[var].rename(Re),
                                           how='outer',
                                           left_index=True,
                                           right_index=True)
            Re_data = Re_data.astype('float64')
            #use interpolation to fill missing data
            Re_data = Re_data.interpolate(method='linear',
                                          axis=0,
                                          limit_area='inside')
            Re_data = Re_data.interpolate(method='linear',
                                          axis=1,
                                          limit_direction='both')
            
            #create a grid to store data
            if i==0:
                alpha_lst = np.array(Re_data.index.drop_duplicates())
                grid = np.zeros([len(M_lst),
                                 len(alpha_lst),
                                 len(Re_lst)])
            #Re_data = Re_data.drop_duplicates()
            grid[i] = Re_data.values
         
        #create interpolation function
        f = si.RegularGridInterpolator((M_lst, alpha_lst, Re_lst), grid)
        
        return f, M_lst, alpha_lst, Re_lst
    
    def interp_func_inviscid(self, var):
        """Reads the polars from the folder and constructs a interpolation
        function using the data.
        
        Args:
            var (str): variable for which the interpolation function is constructed
            
        Returns:
            f (scipy.interpolate.RegularGridInterpolator): the interpolation function
            M_lst (numpy.array): list with available Mach numbers
            alpha_lst (numpy.array): list of available angles of attack
        """   
        
        airfoil_folder = self.airfoil_folder
        airfoil_name = self.airfoil_name

        
        #get file names
        dir_lst = os.listdir(airfoil_folder)
        dir_lst = [i for i in dir_lst if i[:len(airfoil_name)]==airfoil_name]
        
        #read the data
        for i, d in enumerate(dir_lst):
            if airfoil_folder is None:
                df = pd.read_csv(d, header=2)
            else:
                df = pd.read_csv(airfoil_folder+'/'+d, header=2)
            df = df.set_index('alpha')
    
            #get Mach number from file name
            d = d[:-len(d.split('.')[-1])-1]
            spl = d.split('_')
            for s in spl:
                if s[0] == 'M':
                    M = float(s[1:])
            
            #store data
            if i==0:
                M_data = df[var].rename(M)
            else:
                M_data = pd.merge(M_data, df[var].rename(M),
                                  how='outer',
                                  left_index=True,
                                  right_index=True)
        
        #get available Mach numbers and angles of attack                
        M_lst = np.sort(M_data.columns)
        M_data = M_data.reindex(columns=M_lst)
        alpha_lst = M_data.index
        
        #use interpolation to fill missing data
        M_data = M_data.interpolate(method='linear',
                                    axis=0,
                                    limit_area='inside')
        M_data = M_data.interpolate(method='linear',
                                    axis=1,
                                    limit_direction='both')   
        
        #create interpolation function
        f = si.RegularGridInterpolator((M_lst, alpha_lst),
                                       np.transpose(M_data.values))
        
        return f, M_lst, alpha_lst
