import numpy as np
import os
import pandas as pd
import scipy.interpolate as si

import matplotlib.pyplot as plt


class BEM(object):
    
    def __init__(self):
        self.R = 3
        self.B = 6
        self.rR_hub = 0.2
        self.Vinf = 30
        self.omega = 800
        self.rho = 1.225
        self.mu = 1.46e-5
        self.rR = np.linspace(0.2, 1, 12)
        self.beta = 30
        self.theta = None
        self.cR = None
        self.a = 340
        
        self.rR_b = None
        self.cR_b = None
        self.theta_b = None
        
        self.airfoil_folder = 'data_aero_prop/'
        self.tol = 1e-12
        self.N_iter_max = 500
        
        self.f_cl = None
        self.f_cd = None
        
        self.res = None
        self.res_sum = None
        
    def polar_cl(self, idx):
    
        #Get Reynolds number for section
        airfoil_folder = 'data_aero_prop/'
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
    
        #Get Reynolds number for section
        airfoil_folder = 'data_aero_prop/'
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
        f_cl = {}
        f_cd = {}
        for i in range(self.rR_b.size):
            f_cl[i] = self.polar_cl(i+1)
            f_cd[i] = self.polar_cd(i+1)
            
        self.f_cl = f_cl
        self.f_cd = f_cd
            
    
    def calculate_residual(self, psi):
        
        r = self.r
        R = self.R
        B = self.B
        
        beta = np.deg2rad(self.theta+self.beta)
        
        Ua = self.Vinf
        Ut = self.omega*r
        
        U = (Ua**2+Ut**2)**0.5
        
        Wa = 0.5*Ua+0.5*U*np.sin(psi)
        Wt = 0.5*Ut+0.5*U*np.cos(psi)
        
        vt = Ut-Wt
        
        alpha = beta-np.arctan2(Wa, Wt)
        
        W = (Wa**2+Wt**2)**0.5
        
        Re = self.rho*W*self.c/self.mu
        
        lambda_w = self.rR*Wa/Wt
        
#        f_T = -B/2*(1-r/R)*1/lambda_w
#        f_R = -B/2*(r/R-self.rR_hub)*1/lambda_w
        
        f_T = -B/2*(1-r/R)*(1+1/lambda_w**2)**0.5
        f_R = -B/2*(r/R-self.rR_hub)*(1+1/lambda_w**2)**0.5

        F_T = 2/np.pi*np.arccos(np.exp(f_T))
        F_R = 2/np.pi*np.arccos(np.exp(f_R))
        
        #final correction
        F = F_T#*F_R
        F_ = np.interp(r, r[~np.isnan(F)], F[~np.isnan(F)])
        F = np.where(np.isnan(F), F_, F)
        F = np.where(F<1e-7, 1e-5, F)
        
        gamma = vt*(4*np.pi*r)/B*F*(1+(4*lambda_w*R/np.pi*B*r)**2)**0.5
        
#        Re_ = si.pchip(self.rR, Re, extrapolate=True)(self.rR_b)
#        alpha_ = si.pchip(self.rR, alpha, extrapolate=True)(self.rR_b)
#
#        
#        cl_ = np.zeros(self.rR_b.size)
#        for i in range(cl_.size):
#            cl_[i] = self.f_cl[i](Re_[i], np.rad2deg(alpha_[i]))[0]
#        
#        cl = si.pchip(self.rR_b, cl_)(self.rR)
        
        
        cl = np.zeros(r.shape)
        for i in range(r.size):
            cl1 = self.f_cl[i](Re[i], np.rad2deg(alpha[i]))[0]
            cl2 = self.f_cl[i+1](Re[i], np.rad2deg(alpha[i]))[0]
            cl[i] = (cl1+cl2)/2
                    
        R_ = gamma - 0.5*W*self.c*cl
            
        return R_

    def analysis(self):
        
        if self.f_cl is None or self.f_cd is None:
            self.init_polars()
        
        conv = False
        N_iter=0
        
        self.rR = (self.rR_b[1:]+self.rR_b[:-1])/2
        self.cR = si.pchip(self.rR_b, self.cR_b)(self.rR)
        self.theta = si.pchip(self.rR_b, self.theta_b)(self.rR)
        self.r = self.rR*self.R
        self.c = self.cR*self.R
                
        psi = np.arctan2(self.Vinf, self.omega*self.r)

        while not conv:
            N_iter=N_iter+1
        
            R_ = self.calculate_residual(psi)
            
            dpsi = 0.01
            R2_ = self.calculate_residual(psi+dpsi/2)
            R3_ = self.calculate_residual(psi-dpsi/2)
            
            dR_ = R2_-R3_
            
            psi -= R_/(dR_/dpsi)
            val = (psi<0.5*np.pi)*(psi>-0.125*np.pi)
            if np.any(val):
                psi_ = np.interp(self.r, self.r[val], psi[val])
                psi = np.where(psi>0.5*np.pi, psi_, psi)
                psi = np.where(psi<-0.125*np.pi, psi_, psi)
            
            if np.max(np.abs(R_))<self.tol:
                R_ = self.calculate_residual(psi)
                print('Converged after %d iterations' % N_iter)
                print('max residual: %.4e' % np.max(R_))
                conv = True
            
            if N_iter==self.N_iter_max:
                R_ = self.calculate_residual(psi)
                print('Not converged')
                print('max residual: %.4e' % np.max(R_))
                conv = True
                
        cols=['r', 'rR', 'Va', 'Vt', 'alpha', 'cR', 'beta', 'Re', 'Mach', 'Cl',
              'Cd', 'a', 'ap', 'V', 'gamma', 'dL', 'dD', 'dR', 'theta', 'phi',
              'dT', 'dQ', 'dCT', 'dCQ']
        res = pd.DataFrame(columns=cols)
        
        Ua = self.Vinf
        Ut = self.omega*self.r
        
        U = (Ua**2+Ut**2)**0.5
        
        Wa = 0.5*Ua+0.5*U*np.sin(psi)
        Wt = 0.5*Ut+0.5*U*np.cos(psi)
        
        lambda_w = self.rR*Wa/Wt
        
        f_T = -self.B/2*(1-self.rR)*1/lambda_w
        f_R = -self.B/2*(self.rR-self.rR_hub)*1/lambda_w
        
        F_T = 2/np.pi*np.arccos(np.exp(f_T))
        F_R = 2/np.pi*np.arccos(np.exp(f_R))
        
        #final correction
        F = F_T#*F_R
        F_ = np.interp(self.r, self.r[~np.isnan(F)], F[~np.isnan(F)])
        F = np.where(np.isnan(F), F_, F)
        F = np.where(F<1e-7, 1e-5, F)
        
        res['r'] = self.r
        res['rR'] = self.rR
        
        res['beta'] = self.theta+self.beta
        res['theta'] = self.theta
        res['alpha'] = res['beta']-np.rad2deg(np.arctan2(Wa, Wt))
        res['phi'] = res['beta']-res['alpha']
        
        res['V'] = (Wa**2+Wt**2)**0.5
        res['Va'] = 2*(Wa-Ua)*F
        res['Vt'] = 2*(Ut-Wt)*F
        
        res['a'] = res['Va']/2/self.Vinf
        res['ap'] = res['Vt']/2/(self.omega*self.r)
        
        res['cR'] = self.cR
        res['Re'] = self.rho*res['V']*self.c/self.mu
        res['Mach'] = res['V']/self.a
        
        res['gamma'] = res['Vt']*(2*np.pi*self.r)/self.B*(1+(4*lambda_w*self.R/np.pi*self.B*self.r)**2)**0.5
    
        cl = np.zeros(self.r.shape)
        cd = np.zeros(self.r.shape)
        for i in range(self.r.size):
            cl1 = self.f_cl[i](res['Re'].iloc[i], res['alpha'].iloc[i])[0]
            cl2 = self.f_cl[i+1](res['Re'].iloc[i], res['alpha'].iloc[i])[0]
            cl[i] = (cl1+cl2)/2
            cd1 = self.f_cd[i](res['Re'].iloc[i], res['alpha'].iloc[i])[0]
            cd2 = self.f_cd[i+1](res['Re'].iloc[i], res['alpha'].iloc[i])[0]
            cd[i] = (cd1+cd2)/2
        
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
        
        res_sum = {}
        r_b = self.rR_b*self.R
        res_sum['T'] = np.trapz(si.UnivariateSpline(self.r, res['dT'])(r_b), r_b)*self.B
        res_sum['Q'] = np.trapz(si.UnivariateSpline(self.r, res['dQ'])(r_b), r_b)*self.B
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
        
        
        self.res = res
        self.res_sum = res_sum