import numpy as np
import os
import pandas as pd
import scipy.interpolate as si


class BEM(object):
    """A BEM code based on the method used by XROTOR. The induced velocities
    are calculated using the Graded Momentum Formulation (GRAD) of XROTOR.
    """
    
    def __init__(self):
        
        #inputs
        self.R = 3              #radius [m]
        self.B = 6              #number of blades [-]
        self.rR_hub = 0.2       #hub r/R [-]
        self.Vinf = 30          #inflow velocity [m/s]
        self.omega = 800        #rotational speed [rad/s]
        self.rho = 1.225        #air density [kg/m3]
        self.mu = 1.7894e-05    #reference viscosity [Pa*s]
        self.beta = 30          #propeller pitch [deg]
        self.a = 340            #speed of sound [m/s]
              
        self.rR_b = None        #r/R[-]
        self.cR_b = None        #c/R[-]
        self.theta_b = None     #blade twist [deg]
        
        #controil point variables
        self.rR = None          #r/R [-]
        self.cR = None          #c/R [-]
        self.r = None           #radius [m]
        self.c = None           #chord [m]
        self.theta = None       #blade twist [deg]
        
        #settings
        self.airfoil_folder = 'data_aero_prop/' #storage of airfoil polars
        self.tol = 1e-12        #tolerance of the residual
        self.N_iter_max = 500   #maximum number of iterations
        self.N_point = 40       #number of bladewise control points
        
        #lift and drag polar functions
        self.f_cl = None        #lift interpolation function
        self.f_cd = None        #drag interpolation function
        
        #storage variables for results
        self.res = None
        self.res_sum = None
        
    def polar_cl(self, idx):
        """Reads polars from the airfoil folder and constructs an interpolation
        function.
        
        Args:
            idx (int): Radial station number, which corresponds to the file
                names of the polars
        
        Returns:
            f_cl (interp2d instance): Returns the lift coefficient for Reynolds
                number and angle of attack
        """
    
        #Get Reynolds number for section
        airfoil_folder = self.airfoil_folder
        #Get polar files for the right section
        dir_lst = os.listdir(airfoil_folder)
        dir_lst = [i for i in dir_lst if i[2:2+int(np.log10(idx)+1)]==str(idx)]
        
        #Read data for all available Reynolds numbers
        aero_data = {}
        for file in dir_lst:    
            #Read data
            df = pd.read_csv(airfoil_folder+'/'+file,
                             header=6)
            df = df.set_index('alpha')
            #Get the Reynolds number of the data
            f = open(airfoil_folder+'/'+file, 'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Re' in line:
                    line = line.split()
                    Re = float(line[1])
                    aero_data[Re] = df
                    break
        
        #Get available Reynolds numbers
        Re_lst = list(aero_data.keys())
        Re_lst.sort()
    
        #Create one pandas.DataFrame with data for all Reynolds numbers
        for i,Re in enumerate(Re_lst):
            if i==0:
                data = aero_data[Re]['cl'].rename(Re)
            else:
                data = pd.merge(data, aero_data[Re]['cl'].rename(Re),
                                how='outer',
                                left_index=True,
                                right_index=True)

        data = data.astype('float64')

        #Fill NaN values with data using interpolation/extrapolation
        data = data.interpolate(method='linear',
                                axis=0,
                                limit_area='inside')
        data = data.interpolate(method='linear',
                                axis=1,
                                limit_direction='both')
        
        #Interpolate for angle of attack values
        alpha_lst = list(data.index)
        f_cl = si.interp2d(Re_lst, alpha_lst,
                           data.values,
                           kind='linear')
        
        return f_cl
    
    
    def polar_cd(self, idx):
        """Reads polars from the airfoil folder and constructs an interpolation
        function.
        
        Args:
            idx (int): Radial station number, which corresponds to the file
                names of the polars
        
        Returns:
            f_cl (interp2d instance): Returns the drag coefficient for Reynolds
                number and angle of attack
        """
    
        #Get Reynolds number for section
        airfoil_folder = self.airfoil_folder
        #Get polar files for the right section
        dir_lst = os.listdir(airfoil_folder)
        dir_lst = [i for i in dir_lst if i[2:2+int(np.log10(idx)+1)]==str(idx)]
        
        #Read data for all available Reynolds numbers
        aero_data = {}
        for file in dir_lst:    
            #Read data
            df = pd.read_csv(airfoil_folder+'/'+file,
                             header=6)
            df = df.set_index('alpha')
            #Get the Reynolds number of the data
            f = open(airfoil_folder+'/'+file, 'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Re' in line:
                    line = line.split()
                    Re = float(line[1])
                    aero_data[Re] = df
                    break
        
        #Get available Reynolds numbers
        Re_lst = list(aero_data.keys())
        Re_lst.sort()

        #Create one pandas.DataFrame with data for all Reynolds numbers
        for i,Re in enumerate(Re_lst):
            if i==0:
                data = aero_data[Re]['cd'].rename(Re)
            else:
                data = pd.merge(data, aero_data[Re]['cd'].rename(Re),
                                how='outer',
                                left_index=True,
                                right_index=True)

        data = data.astype('float64')
        #Fill NaN values with data using interpolation/extrapolation
        data = data.interpolate(method='linear',
                                axis=0,
                                limit_area='inside')
        data = data.interpolate(method='linear',
                                axis=1,
                                limit_direction='both')
        
        #Interpolate for angle of attack values
        alpha_lst = list(data.index)
        f_cd = si.interp2d(Re_lst, alpha_lst,
                           data.values,
                           kind='linear')
        
        return f_cd
    
    def init_polars(self):
        """Constructs automatically the polar functions by reading polars
        for each radial station from the airfoil folder.
        """
        
        f_cl = {}
        f_cd = {}
        
        #Construct interpolation functions
        for i in range(self.rR_b.size):
            f_cl[i] = self.polar_cl(i+1)
            f_cd[i] = self.polar_cd(i+1)
            
        self.f_cl = f_cl
        self.f_cd = f_cd
            
    def calculate_residual(self, psi):
        """Calculates the residual. Circulation is calculated using the local
        velocity at the blade and the lift coefficient. This is compared to
        the circulation calculated with the tangential induced velocity using
        Helmholtz relation.
        
        Args:
            psi (float): dummy variable used for iteration scheme
        """
        
        #propeller geometry
        r = self.r
        R = self.R
        B = self.B
        beta = np.deg2rad(self.theta+self.beta)
        
        #inflow velocities
        Ua = self.Vinf
        Ut = self.omega*r
        
        #total inflow velocity
        U = (Ua**2+Ut**2)**0.5
        
        #total axial and tangential velocities
        Wa = 0.5*Ua+0.5*U*np.sin(psi)
        Wt = 0.5*Ut+0.5*U*np.cos(psi)
        
        #tangential induced velocity
        vt = Ut-Wt
        
        #angle of attack [rad]
        alpha = beta-np.arctan2(Wa, Wt)
        
        #total velocity
        W = (Wa**2+Wt**2)**0.5
        M = W/self.a
        
        #local Reynolds number
        Re = self.rho*W*self.c/self.mu
        
        #apply modified Prandtl's tip loss factor
        lambda_w = self.rR*Wa/Wt
        f_T = -B/2*(1-r/R)*1/lambda_w
        F_T = 2/np.pi*np.arccos(np.exp(f_T))
        
        #final correction, avoid singularities
        F = F_T
        F_ = np.interp(r, r[~np.isnan(F)], F[~np.isnan(F)])
        F = np.where(np.isnan(F), F_, F)
        F = np.where(F<1e-7, 1e-5, F)
        
        #calculate circulation using Helmholtz' relation
        gamma = vt*(4*np.pi*r)/B*F*(1+(4*lambda_w*R/(np.pi*B*r))**2)**0.5
        
        #calculate lift coefficient
        cl = np.zeros(r.shape)

        r_b = self.rR_b*R
        for i, r_ in enumerate(r):
            idx = np.argsort(np.abs(r_b-r_))
            if idx[0]>idx[1]:
                i1 = idx[1]
                i2 = idx[0]
            else:
                i1 = idx[0]
                i2 = idx[1]
            cl1 = self.f_cl[i1](Re[i], np.rad2deg(alpha[i]))[0]
            cl2 = self.f_cl[i2](Re[i], np.rad2deg(alpha[i]))[0]
            cl[i] = np.interp(r_, [r_b[i1], r_b[i2]], [cl1, cl2])
        
        #apply Prandtl-Glauert correction
        cl = np.where(M<0.7, cl/(1-M**2)**0.5, cl/(1-0.7**2)**0.5)
        
        #residual
        R_ = gamma - 0.5*W*self.c*cl
            
        return R_

    def analysis(self):
        """Runs and iterative scheme that reduces the residual in circulation
        at the control points to zero.
        """
        
        #create polars
        if self.f_cl is None or self.f_cd is None:
            self.init_polars()
        
        conv = False
        N_iter=0
        
        #create cosine spacing
        cos_b = (1-np.cos(np.linspace(0, np.pi, self.N_point+1)))/2
        rR_b = self.rR_hub+cos_b*(1-self.rR_hub)
        
        #calculate geometry at the control points
        self.rR = (rR_b[1:]+rR_b[:-1])/2
        self.cR = si.pchip(self.rR_b, self.cR_b)(self.rR)
        self.theta = si.pchip(self.rR_b, self.theta_b)(self.rR)
        self.r = self.rR*self.R
        self.c = self.cR*self.R
        
        #first guess of dummy variable psi
        psi = np.arctan2(self.Vinf, self.omega*self.r)

        while not conv:
            N_iter=N_iter+1
            
            #calculate residual
            R_ = self.calculate_residual(psi)
            
            #calculate gradient
            dpsi = 0.01
            R2_ = self.calculate_residual(psi+dpsi/2)
            R3_ = self.calculate_residual(psi-dpsi/2)
            dR_ = R2_-R3_
            
            #update psi
            psi -= R_/(dR_/dpsi)
            
            #avoid psi values outside boundaries
            val = (psi<0.5*np.pi)*(psi>-0.125*np.pi)
            if np.any(val):
                psi_ = np.interp(self.r, self.r[val], psi[val])
                psi = np.where(psi>0.5*np.pi, psi_, psi)
                psi = np.where(psi<-0.125*np.pi, psi_, psi)
            
            #check for convergence
            if np.max(np.abs(R_))<self.tol:
                R_ = self.calculate_residual(psi)
                print('Converged after %d iterations' % N_iter)
                print('max residual: %.4e' % np.max(np.abs(R_)))
                conv = True
                conv_f = True
            
            #check if maximum number of iterations is reached
            if N_iter==self.N_iter_max:
                R_ = self.calculate_residual(psi)
                print('Not converged')
                print('max residual: %.4e' % np.max(np.abs(R_)))
                conv = True
                conv_f = False
        
        #save results
        self.save_results(psi, conv_f)
    
    def save_results(self, psi, conv_f):
        """Saves the results of the analysis
        
        Args:
            psi (float): converged value of dummy variable psi
        """
        
        #create dataframe
        cols=['r', 'rR', 'Va', 'Vt', 'alpha', 'cR', 'beta', 'Re', 'Mach', 'Cl',
              'Cd', 'V', 'gamma', 'dL', 'dD', 'dR', 'theta', 'phi', 'dT', 'dQ',
              'dCT', 'dCQ']
        res = pd.DataFrame(columns=cols)
        
        #propeller geometry
        r = self.r
        R = self.R
        B = self.B
        beta = np.deg2rad(self.theta+self.beta)
        
        #inflow velocities
        Ua = self.Vinf
        Ut = self.omega*r
        
        #total inflow velocity
        U = (Ua**2+Ut**2)**0.5
        
        #total axial and tangential velocities
        Wa = 0.5*Ua+0.5*U*np.sin(psi)
        Wt = 0.5*Ut+0.5*U*np.cos(psi)
        
        #axial and tangential induced velocities
        va = Wa-Ua
        vt = Ut-Wt
        
        #angle of attack [rad]
        alpha = beta-np.arctan2(Wa, Wt)
        
        #total velocity
        W = (Wa**2+Wt**2)**0.5
        M = W/self.a
        
        #local Reynolds number
        Re = self.rho*W*self.c/self.mu
        #print(Re)
        
        #apply modified Prandtl's tip loss factor
        lambda_w = self.rR*Wa/Wt
        f_T = -B/2*(1-r/R)*1/lambda_w
        F_T = 2/np.pi*np.arccos(np.exp(f_T))
        
        #final correction, avoid singularities
        F = F_T
        F_ = np.interp(r, r[~np.isnan(F)], F[~np.isnan(F)])
        F = np.where(np.isnan(F), F_, F)
        F = np.where(F<1e-7, 1e-5, F)
        
        #circulation
        gamma = vt*(4*np.pi*r)/B*F*(1+(4*lambda_w*R/(np.pi*B*r))**2)**0.5

        #calculate cl and cd        
        cl = np.zeros(r.shape)
        cd = np.zeros(r.shape)
        r_b = self.rR_b*R
        for i, r_ in enumerate(r):
            idx = np.argsort(np.abs(r_b-r_))
            if idx[0]>idx[1]:
                i1 = idx[1]
                i2 = idx[0]
            else:
                i1 = idx[0]
                i2 = idx[1]
            cl1 = self.f_cl[i1](Re[i], np.rad2deg(alpha[i]))[0]
            cl2 = self.f_cl[i2](Re[i], np.rad2deg(alpha[i]))[0]
            cl[i] = np.interp(r_, [r_b[i1], r_b[i2]], [cl1, cl2])
            cd1 = self.f_cd[i1](Re[i], np.rad2deg(alpha[i]))[0]
            cd2 = self.f_cd[i2](Re[i], np.rad2deg(alpha[i]))[0]
            cd[i] = np.interp(r_, [r_b[i1], r_b[i2]], [cd1, cd2])
        
        #apply Prandtl-Glauert correction
        cl = np.where(M<0.7, cl/(1-M**2)**0.5, cl/(1-0.7**2)**0.5)
        cd = np.where(M<0.7, cd/(1-M**2)**0.5, cd/(1-0.7**2)**0.5)
        
        #save radial station results to dataframe
        res['r'] = self.r
        res['rR'] = self.rR
        
        res['beta'] = np.rad2deg(beta)
        res['theta'] = self.theta
        res['alpha'] = np.rad2deg(alpha)
        res['phi'] = res['beta']-res['alpha']
        
        res['V'] = W
        res['Va'] = va*F
        res['Vt'] = vt*F
        
        res['Va_i'] = va
        res['Vt_i'] = vt
        
        res['cR'] = self.cR
        res['Re'] = Re
        res['Mach'] = M
        
        res['gamma'] = gamma
        
        res['Cl'] = cl
        res['Cd'] = cd

        res['dL'] = cl*0.5*self.rho*res['V']**2*self.c
        res['dD'] = cd*0.5*self.rho*res['V']**2*self.c
        res['dR'] = (res['dL']**2+res['dD']**2)**0.5
        
        phi = np.deg2rad(res['phi'].to_numpy())
        res['dT'] = res['dL']*np.cos(phi)-res['dD']*np.sin(phi)
        res['dQ'] = (res['dL']*np.sin(phi)+res['dD']*np.cos(phi))*res['r']
          
        norm_T = (self.rho*(self.omega/(2*np.pi))**2*(2*self.R)**4)
        norm_Q = (self.rho*(self.omega/(2*np.pi))**2*(2*self.R)**5)
        #thrust and torque coefficients
        res['dCT'] = res['dT']/norm_T
        res['dCQ'] = res['dQ']/norm_Q
        
        #save total propeller results
        res_sum = {}
        r_b = self.rR_b*self.R
        res_sum['T'] = np.trapz(si.pchip(self.r, res['dT'])(r_b), r_b)*self.B
        res_sum['Q'] = np.trapz(si.pchip(self.r, res['dQ'])(r_b), r_b)*self.B
        res_sum['P'] = res_sum['Q']*self.omega
        res_sum['Ct'] = res_sum['T']/norm_T
        res_sum['Cq'] = res_sum['Q']/norm_Q
        res_sum['Cp'] = res_sum['Cq']*2*np.pi
        res_sum['J'] = np.pi*self.Vinf/(self.omega*self.R)
        res_sum['eta'] = res_sum['Ct']/res_sum['Cp']*res_sum['J']
        res_sum['R'] = self.R
        res_sum['r_hub'] = self.R*self.rR_hub
        res_sum['rho'] = self.rho
        res_sum['a'] = self.a
        res_sum['Vinf'] = self.Vinf
        res_sum['N_blades'] = self.B
        res_sum['rpm'] = self.omega/(2*np.pi)*60
        res_sum['converged'] = conv_f
        
        self.res = res
        self.res_sum = res_sum
 