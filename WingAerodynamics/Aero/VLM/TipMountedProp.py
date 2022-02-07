import copy
import numpy as np
import pandas as pd
import sys
    
if 'BEM/' not in sys.path:
    sys.path.append('BEM/')
    
if 'VLM/' not in sys.path:
    sys.path.append('VLM/')
    
from create_unsteady_rotor_model import create_unsteady_rotor_model


class PropWing(object):
    """Uses integrates the analysis of a tip-mounted propeller-wing system"""
    
    def __init__(self, vlm, bem, prop_loc, prop_angle=0, rotation_dir=1):
        """Args:
            vlm (VLM class): vlm instance with wing already addad
            bem (BEM class): bem instance with propeller geometry added
            prop_loc (numpy.array): x,y,z vector of propeller location in the
                wing coordinate system
            prop_angle (float): angle of attack of the propeller in deg
            rotation_dir (float): 1 for outborad-up, -1 for inboard-up
        """
        
        #check flow properties
        if vlm.rho != bem.rho:
            bem.rho = vlm.rho
            print('Different values found for density')
        if vlm.Vinf != bem.Vinf:
            bem.Vinf = vlm.Vinf
            print('Different values found for freestream velocity')
        if vlm.mu != bem.mu:
            bem.mu = vlm.mu
            print('Different values found for dynamic viscosity')
        if vlm.a != bem.a:
            bem.a = vlm.a
            print('Different values found for speed of sound')
        
        self.vlm = vlm
        self.bem = bem
        self.urot = None
        
        self.rotation_dir = rotation_dir
        
        self.prop_angle = prop_angle
        self.prop_loc = prop_loc
        
        prop_angle_ = -np.deg2rad(prop_angle)
        self.prop_rot = np.array([[np.cos(prop_angle_), 0, np.sin(prop_angle_)],
                                   [0, 1, 0],
                                   [-np.sin(prop_angle_), 0, np.cos(prop_angle_)]])
    
        self.settings = {'dJ':                  0.5,
                         'N_prop_map':          3,
                         'slipstream_length':   2,
                         'dlbda':               0.1,
                         'dlb':                 0.005,
                         'p_max':               10,
                         'lbda_max':            5,
                         'N_time_step':         30,
                         'N_slipstream_x':      12,
                         'conway_s_max':        500,
                         'conway_ds':           0.1}
        
        self.N_iter_max = 4
        
    def prop_to_wing(self, P):
        """Coordinate transform from propeller to wing coordinate system
        
        Args:
            P (numpy.array): point in propeller coordinate system
        
        Returns:
            p (numpy.array): point in wing coordinate system
        """
        
        p = P+self.prop_loc
        p = np.matmul(self.prop_rot, p)
        
        return p

    def wing_to_prop(self, P):
        """Coordinate transform from wing to propellercoordinate system
        
        Args:
            P (numpy.array): point in wing coordinate system
        
        Returns:
            p1,p2 (numpy.array): points in propeller coordinate system for the 
                left and right propeller
        """
        
        p1 = P-self.prop_loc
        p2 = P-self.prop_loc*np.array([1,-1,1])
        p1 = np.matmul(self.prop_rot.T, p1) #right prop
        p2 = np.matmul(self.prop_rot.T, p2) #left prop
        p2 = p2*np.array([1,-1,1])
        
        return p1, p2
    
    def analysis(self, alpha):
        """Analysis the propeller-wing system. Solution will use iterations
        for to converge the propller-wing induced and wing-propeller induced
        velocities
        
        Args:
            alpha (float): angle of attack in deg
        """
        
        bem = self.bem
        vlm = self.vlm
        
        #pass settings
        for key in vlm.jet_corr_settings.keys():
            vlm.jet_corr_settings[key] = self.settings[key]
        
        #geometry data
        Vinf = vlm.Vinf
        R = bem.R
        alpha_wing = alpha
        alpha_prop = alpha+self.prop_angle
        
        bem.analysis()
        
        if not bem.res_sum['converged']:
            conv_bem = False
        else:
            conv_bem = True
            #create a propeller map
            if self.urot is None:
                dJ = self.settings['dJ']
                N_prop_map = self.settings['N_prop_map']
                urot = create_unsteady_rotor_model(bem, dJ, N_prop_map)
                urot.N_time_step = self.settings['N_time_step']
                urot.sign_rotation = self.rotation_dir
                self.urot = urot
            else:
                urot = self.urot
            
            #slipstream geometry
            sl = self.settings['slipstream_length']
            dsl = bem.Vinf/(bem.omega/(2*np.pi))/self.settings['N_slipstream_x']
            #prevent slipstream with only one x location
            if dsl>=sl:
                dsl = sl/2.001
            x = np.arange(0, sl, dsl)
            z = np.zeros(x.shape)
            #slipstream settings
            urot.slipstream.cnw.ds = self.settings['conway_ds']
            urot.slipstream.cnw.s_max = self.settings['conway_s_max']
                
            #set propeller inflow field
            u_aoa = np.cos(np.deg2rad(alpha_prop))
            w_aoa = np.sin(np.deg2rad(alpha_prop))
            grid_y = np.linspace(-R, R, 50)
            grid_z = np.linspace(-R, R, 50)
            u_in = np.ones((grid_y.size, grid_z.size))*u_aoa
            v_in = np.zeros((grid_y.size, grid_z.size))
            w_in = np.ones((grid_y.size, grid_z.size))*w_aoa
            
            #data frame to store convergence data
            conv = pd.DataFrame(columns=['dCT', 'dCQ', 'dCL', 'dCD'])
            
            ct0, cq0, cl0, cd0 = 10, 10, 10, 10
            
            tmp_dct = {}
            conv_iter = False
            for N in range(self.N_iter_max):
                
                #analyse prop
                urot.inflow.grid_input(grid_y/R, grid_z/R, u_in, v_in, w_in)
                urot.analysis()
                
                urot.calculate_circulation(bem.res['Va_i'].to_numpy(),
                                           bem.res['Vt_i'].to_numpy(),
                                           bem.res['rR'].to_numpy())
                
                urot.init_slipstream(x, z)
        
                #propeller induced velocities on wing
                Vpert = np.zeros((vlm.Vpert.shape[0], 3))
                Vperti = np.zeros((vlm.Vpert.shape[0], 3))
                
                for i in range(vlm.Vpert.shape[0]):
                    P_ = np.array([vlm.Vpert.iloc[i]['x'],
                                   vlm.Vpert.iloc[i]['y'],
                                   0])
                    P1, P2 = self.wing_to_prop(P_)
                    Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
                    Vt2, Vw2 = urot.slipstream.induced_velocity(P2)
                    
                    V1 = np.matmul(self.prop_rot, Vt1)
                    V2 = np.matmul(self.prop_rot, Vt2)
                    
                    V2 = V2*np.array([1,-1,1]) #left prop
                    
                    Vpert[i] = V1+V2
                
                    #velocities for trefftz plane
                    xi = np.max(urot.slipstream.x)/2
                    x0 = P1[0]
                    if urot.slipstream.r_r0 is not None:
                        
                        sl_r_r0 = urot.slipstream.r_r0
                        sl_x = urot.slipstream.x.flatten()
                        sl_z = urot.slipstream.z.flatten()
                        i_xi = np.argmin(np.abs(sl_x-xi))
                        i_x0 = np.argmin(np.abs(sl_x-x0))
                        
                        z0 = np.interp(x0, sl_x, sl_z)
                        zi = np.interp(xi, sl_x, sl_z)
                        r_r0 = sl_r_r0[i_x0].flatten()
                        
                        r0 = (P1[1]**2+P1[2]**2)**0.5
                        
                        if r0>np.max(r_r0)*R:
                            r_r = 1
                        else:
                            i_r = np.argmin(np.abs(r0-r_r0*R))
                            r_r = sl_r_r0[i_xi, 0, i_r]/r_r0[i_r]
                        
                        P1[1] = P1[1]*r_r
                        P1[2] = (P1[2]-z0)*r_r+zi
                        P2[1] = P2[1]*r_r
                        P2[2] = (P2[2]-z0)*r_r+zi
                    
                    P1[0] = xi
                    P2[0] = xi
                    
                    Vt1, Vw1 = urot.slipstream.induced_velocity(P1)
                    Vt2, Vw2 = urot.slipstream.induced_velocity(P2)
                    
                    V1 = np.matmul(self.prop_rot, Vt1)
                    V2 = np.matmul(self.prop_rot, Vt2)
                    
                    V2 = V2*np.array([1,-1,1]) #left prop
                    
                    Vperti[i] = V1+V2
                
                vlm.Vpert['Vx']   = Vpert[:, 0]
                vlm.Vpert['Vz']   = Vpert[:, 2]
                vlm.Vpert['Vzi']   = Vperti[:, 2]
                
    #            b2 = np.max(np.array(vlm.Points)[:, 1])
    #            R = bem.R
    #            vlm.Vpert['turb'] = (b2-np.abs(vlm.Vpert['y'].to_numpy())<R)
                vlm.Vpert['turb'] = True
               
                #analyse wing
                vlm.G = 0
                vlm.alpha = alpha_wing
                vlm.vlm_visc()
                
                #wing induced velocities on propeller
                for iy, y_ in enumerate(grid_y):
                    for iz, z_ in enumerate(grid_z):
                        P_ = np.array([0, y_, z_])
                        P = self.prop_to_wing(P_)
                        
                        V = vlm.induced_velocity(P)
                        V = np.matmul(self.prop_rot.T, V)
                        
                        u_in[iy, iz] = V[0]/Vinf+u_aoa
                        v_in[iy, iz] = V[1]/Vinf
                        w_in[iy, iz] = V[2]/Vinf+w_aoa
                        
                #slipstream deflection
                for i in range(z.size-1):
                    P_ = np.array([x[i], 0, 0])
                    P = self.prop_to_wing(P_)
                    w = vlm.induced_velocity(P)
                    w = np.matmul(self.prop_rot.T, w)
                    
                    dx = x[i+1]-x[i]
                    z[i+1] = z[i]+dx*w[2]/Vinf
                
                #convergence data
                ct1 = urot.prop_us['integral_CT']
                cq1 = urot.prop_us['integral_CQ']
                cl1 = vlm.res['CL']
                cd1 = vlm.res['CD']
                conv.loc[N] = [np.abs(ct0-ct1), np.abs(cq0-cq1), np.abs(cl0-cl1), np.abs(cd0-cd1)]
                
                tmp_dct['urot_%d' % N] = copy.deepcopy(urot)
                tmp_dct['vlm_%d' % N] = copy.deepcopy(vlm)
                
                if conv.loc[N, 'dCT']<1e-4 and conv.loc[N, 'dCQ']<1e-4 and conv.loc[N, 'dCL']<1e-3 and conv.loc[N, 'dCD']<1e-4:
                    print(N)
                    conv_iter = True
                    break
                else:
                    ct0 = ct1
                    cq0 = cq1
                    cl0 = cl1
                    cd0 = cd1            
            
            if not conv_iter:
                N = conv['dCD'].idxmin()
            urot = tmp_dct['urot_%d' % N]
            vlm = tmp_dct['vlm_%d' % N]
            conv = conv.iloc[:N+1]
            
            self.urot = urot
            self.vlm = vlm
            self.conv = conv
            
        return conv_bem
