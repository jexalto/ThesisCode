import numpy as np
import scipy.interpolate as si

from Inflow import Inflow
from SensitivityDist import SensitivityDist
from actuator_disk import actuator_disk
from sears import sears_function
from SlipstreamTubeModel import SlipstreamTubeModel


class UnsteadyRotor(object):
    """Calculate the propeller forces for an non-uniform inflow in a 
    quasi-steady and unsteady way.
    
    Note: Inputs and outputs are given in a right-handed coordinate system,
        with x in the flow direction, y to the left (when looking to the front
        of the propeller) and z upwards. Internally a left-handed coordinate
        system is used with the y axis negative from the right-handed system.
    
    Original code written by Nando van Arnhem in MATLAB
    """
    
    def __init__(self):        
        
        #inputs
        self.Dp = 0.2032        #propeller diameter [m]
        self.Rh = 0.02185       #propeller hub radius [m]
        self.B = 6              #propeller number of blades
        self.n = 115            #propeller rps [s-1]
        self.Vinf = 40          #freestream velocity [m/s]
        self.rho_inf = 1.225    #freestream air density [kg/m3]
        self.a = 340.2941       #speed of sound [m/s]
        self.sign_rotation = -1 #rotation direction, 1=clockwise, -1=counter clockwise
        
        self.rR = None          #propeller r/R distribution
        self.cR = None          #propeller c/R distribution
        
        self.constant_V = SensitivityDist()     #object for propeller sensitivity
        self.constant_RPS = SensitivityDist()   #object for propeller sensitivity
        self.inflow = Inflow()                  #object for inflow velocity
        self.slipstream = SlipstreamTubeModel() #object for slipstreamtube model
        
        #settings
        self.N_time_step = 36   #time step of the analysis [deg]
        
    def set_chord(self, rR_, cR_):
        """Interpolates a c/R distribution for a certain r/R to the r/R used
        in the analysis using interpolation
        
        Args:
            rR_ (np.array): radial stations where the chord is defined
            cR_ (np.array): chord to radius values at r/R
        """
        
        spl_c = si.pchip(rR_, cR_)        
        self.cR = spl_c(self.rR)

    def analysis(self):
        """For the defined inflow, propeller sensitivity maps and propeller
        properties, the propeller will be analyzed.
            Firstly the forces for the isolated propeller will be calculated.
            Secondly using a difference in axial and tangential advance ratio a
                quasi-steady solution will be added. This solution consists of
                a delta thrust and torque which is added to the isolated forces
            Thirdly an unsteady analysis is done where the delta thrust and
                torque from the quasi-steady analysis is recalculated using
                an unsteady correction based on the Sears' function.
        
        The most important results are saved in the following attributes:
            prop_iso: isolated propeller results
            prop_qs: quasi-steady propeller results
            prop_un: unsteady propeller results
        """
        
        #get propeller properties
        rR = self.rR
        Vinf = self.Vinf
        n = self.n
        Dp = self.Dp
        
        #calculate advance ratio
        J_iso = Vinf/(n*Dp)
        self.J_iso = J_iso
        
        #calculate values used to normalize forces and torque
        self.norm_CT = self.rho_inf*n**2*Dp**4
        self.norm_CQ = self.rho_inf*n**2*Dp**5
        
        #get the isolated propeller thrust and torque distribution
        dCT, dCQ = self.constant_V.interpolate(float(J_iso), rR)        
        #calculate the other isolated propeller forces
        self.calculate_prop_iso(dCT, dCQ) #results saved in prop_iso
        
        #define the timesteps as angles
        phi = np.linspace(0, 2*np.pi, self.N_time_step+1)
        phi_ = phi.reshape((phi.size, 1))
        
        #create a grid
        y = -rR*np.sin(phi_)
        z = rR*np.cos(phi_)
        self.phi = phi
        self.y, self.z = y, z
        theta = self.sign_rotation*np.arctan2(y, z)
        
        #determine the velocity field at the propeller axis using the
        #x and y grid
        u,v,w = self.inflow.interpolation(y, z)
        #multiply the non-dimensional inflow field with the freestream velocity
        #to get the absolute velocity field
        u_abs = u*Vinf
        v_abs = v*Vinf
        w_abs = w*Vinf
        self.u, self.v, self.w = u_abs, v_abs, w_abs

        #determine equivalent 'rotational speed' due to the swirling inflow
        Vt = v_abs*np.cos(theta)+w_abs*np.sin(theta)
        delta_n = Vt/(((y**2+z**2)**0.5)*Dp/2)/(2*np.pi)
        
        #change in advance ratio due to the axial and tangential nonuniform inflow
        delta_J_a = u_abs/(n*Dp)-J_iso
        delta_J_t = Vinf/((n+delta_n)*Dp)-J_iso
        
        #interpolate the thrust and torque coefficients for all radial positions
        #and their respective advance ratios
        dCT_a, dCQ_a = self.constant_RPS.interpolate(delta_J_a+J_iso, rR)
        dCT_t, dCQ_t = self.constant_V.interpolate(delta_J_t+J_iso, rR)
        
        #calculate the difference in distributed thrust and torque with
        #respect to the isolated propeller
        delta_dT_a      = dCT_a*self.norm_CT-self.prop_iso['dT']
        delta_dQ_a      = dCQ_a*self.norm_CQ-self.prop_iso['dQ']
        delta_dT_t      = dCT_t*self.norm_CT-self.prop_iso['dT']
        delta_dQ_t      = dCQ_t*self.norm_CQ-self.prop_iso['dQ']
        
        delta_dT = [delta_dT_a, delta_dT_t]
        delta_dQ = [delta_dQ_a, delta_dQ_t]
        
        #calculate the quasi-steady propeller forces
        prop_qs = self.calculate_propeller_forces(delta_dT, delta_dQ)
        #save the advance ratios
        prop_qs['delta_J_a'] = delta_J_a
        prop_qs['delta_J_t'] = delta_J_t
        prop_qs['Vt'] = Vt
        #save result to attribute
        self.prop_qs = prop_qs
        
        #calculate the unsteady propeller forces
        self.prop_us = self.unsteady_correction()
        self.prop_us['Vt'] = Vt

    def calculate_prop_iso(self, dCT, dCQ):
        """Calculates the forces on the isolated propeller based on a
        distributed thrust and torque along the propeller radius. The result is
        stored in the prop_iso attribute.
        
        Args:
            dCT (numpy.array): distributed thrust coefficient along the
                propeller blade
            dCQ (numpy.array): distributed torque coefficient along the
                propeller blade
        """
        
        #get propeller properties
        B = self.B
        n = self.n
        J_iso = self.J_iso
        Dp = self.Dp
        rR = self.rR
        
        #values used to normalize forces and torque
        norm_CT = self.norm_CT
        norm_CQ = self.norm_CQ
        
        #initiate dictionary
        prop_iso = {}
        
        #radius distribution used
        prop_iso['rR'] = rR
        prop_iso['r'] = rR*Dp/2
        
        #distributed properties
        prop_iso['dCT'] = dCT
        prop_iso['dCQ'] = dCQ
        prop_iso['dT'] = prop_iso['dCT']*norm_CT
        prop_iso['dQ'] = prop_iso['dCQ']*norm_CQ
        prop_iso['deta'] = prop_iso['dT']/prop_iso['dQ']*J_iso*Dp/(4*np.pi)
        
        #blade properties
        prop_iso['blade_T'] = np.trapz(prop_iso['dT'], rR*Dp/2)
        prop_iso['blade_Q'] = np.trapz(prop_iso['dQ'], rR*Dp/2)
        prop_iso['blade_CT'] = prop_iso['blade_T']/norm_CT
        prop_iso['blade_CQ'] = prop_iso['blade_Q']/norm_CQ
        prop_iso['blade_CP'] = 2*np.pi*prop_iso['blade_CQ']
        prop_iso['blade_eta'] = prop_iso['blade_CT']/prop_iso['blade_CP']*J_iso
        
        #total properties
        prop_iso['T'] = prop_iso['blade_T']*B
        prop_iso['Q'] = prop_iso['blade_Q']*B
        prop_iso['P'] = prop_iso['blade_Q']*2*np.pi*n
        prop_iso['CT'] = prop_iso['blade_CT']*B
        prop_iso['CQ'] = prop_iso['blade_CQ']*B
        prop_iso['CP'] = prop_iso['blade_CP']*B
        prop_iso['eta'] = prop_iso['blade_eta']
        
        #set class attribute
        self.prop_iso = prop_iso

    def calculate_propeller_forces(self, delta_dT, delta_dQ):
        """Calculates the propeller forces based on a distributed difference in
        thrust and torque based on the isolated propeller. For the quasi-steady
        case two delta_dT and delta_dQ are given, one due to a change in axial
        advance ratio and one due to change in tangential advance ratio. In
        this case the variables are given in a list.
        
        Args:
            delta_dT (list/numpy.array): distributed difference in thrust.
                In the quasi-steady case this is a list: [delta_dT_a, delta_dT_t]
                In the unsteady case this is a array of the values
            delta_dQ (list/numpy.array): distributed difference in torque.
                In the quasi-steady case this is a list: [delta_dQ_a, delta_dQ_t]
                In the unsteady case this is a array of the value
        
        Returns:
            prop (dict): contains all the calculated propeller forces
        """
        
        #check if isolated propeller has been calculated
        if not hasattr(self, 'prop_iso'):
            s = 'Error: Isolated propeller not yet calculated'
            raise Exception(s)
        
        #check the type of input, set case to quasi-steady or unsteady
        if type(delta_dT)==list:
            prop_type = 'qs'
        else:
            prop_type = 'us'
        
        #get propeller properties
        B = self.B
        J_iso = self.J_iso
        Dp = self.Dp
        rR = self.rR
        phi = self.phi
        
        #values used to normalize forces and torque
        norm_CT = self.norm_CT
        norm_CQ = self.norm_CQ
        
        #isolated propeller
        prop_iso = self.prop_iso
        
        r = rR*Dp/2 #actual radius [m]
        #the array phi, but as a column array
        phi_ = phi.reshape((phi.size, 1))
        
        #initiate dictionary
        prop = {}
        
        #set general propeller properties
        prop['phi'] = phi
        prop['y']   = -self.y #convert to right-handed coordinate system
        prop['z']   = self.z
        prop['u']   = self.u
        prop['v']   = self.v
        prop['w']   = self.w
        prop['r']   = r
        
        #determine and set the change in distributed thrust and torque
        if prop_type=='qs':
            prop['delta_dT_a'] = delta_dT[0]
            prop['delta_dQ_a'] = delta_dQ[0]
            prop['delta_dT_t'] = delta_dT[1]
            prop['delta_dQ_t'] = delta_dQ[1]
            
            prop['delta_dT'] = prop['delta_dT_a']+prop['delta_dT_t']
            prop['delta_dQ'] = prop['delta_dQ_a']+prop['delta_dQ_t']
        else:
            prop['delta_dT'] = delta_dT
            prop['delta_dQ'] = delta_dQ
        
        #calculate distributed forces          
        prop['dT']              = prop_iso['dT']+prop['delta_dT']
        prop['dQ']              = prop_iso['dQ']+prop['delta_dQ']
          
        prop['delta_dFz']       = prop['delta_dQ']*np.sin(phi_)/r
        prop['delta_dFy']       = prop['delta_dQ']*np.cos(phi_)/r
        
        if prop_type=='qs':
            prop['delta_dFz_a'] = prop['delta_dQ_a']*np.sin(phi_)/r
            prop['delta_dFz_t'] = prop['delta_dQ_t']*np.sin(phi_)/r      
            prop['delta_dFy_a'] = prop['delta_dQ_a']*np.cos(phi_)/r
            prop['delta_dFy_t'] = prop['delta_dQ_t']*np.cos(phi_)/r
        
        prop['delta_dCT']       = prop['delta_dT']/norm_CT
        prop['dCT']             = prop['dT']/norm_CT
        prop['delta_dCQ']       = prop['delta_dQ']/norm_CQ
        prop['dCQ']             = prop['dQ']/norm_CQ
        
        if prop_type=='qs':
            prop['delta_dCT_a'] = prop['delta_dT_a']/norm_CT
            prop['delta_dCT_t'] = prop['delta_dT_t']/norm_CT
            prop['delta_dCQ_a'] = prop['delta_dQ_a']/norm_CQ
            prop['delta_dCQ_t'] = prop['delta_dQ_t']/norm_CQ

        prop['delta_dCY']       = prop['delta_dFy']/norm_CT
        prop['delta_dCZ']       = prop['delta_dFz']/norm_CT
        
        if prop_type=='qs':
            prop['delta_dCY_a'] = prop['delta_dFy_a']/norm_CT
            prop['delta_dCY_t'] = prop['delta_dFy_t']/norm_CT
            prop['delta_dCZ_a'] = prop['delta_dFz_a']/norm_CT
            prop['delta_dCZ_t'] = prop['delta_dFz_t']/norm_CT
            
        prop['deta']            = prop['dCT']/(prop['dCQ']*2*np.pi)*J_iso
        prop['delta_deta']      = prop['deta']-prop_iso['deta']
        
        #calculate blade forces
        prop['blade_delta_T']   = np.trapz(prop['delta_dT'], r)
        prop['blade_delta_Q']   = np.trapz(prop['delta_dQ'], r)
        prop['blade_delta_Fy']  = np.trapz(prop['delta_dFy'], r)
        prop['blade_delta_Fz']  = np.trapz(prop['delta_dFz'], r)
        
        prop['blade_delta_CT']  = prop['blade_delta_T']/norm_CT
        prop['blade_delta_CQ']  = prop['blade_delta_Q']/norm_CQ
        prop['blade_delta_CY']  = prop['blade_delta_Fy']/norm_CT
        prop['blade_delta_CZ']  = prop['blade_delta_Fz']/norm_CT
        
        if prop_type=='qs':
            prop['blade_delta_T_a']  = np.trapz(prop['delta_dT_a'], r) 
            prop['blade_delta_Q_a']  = np.trapz(prop['delta_dQ_a'], r)
            prop['blade_delta_T_t']  = np.trapz(prop['delta_dT_t'], r)
            prop['blade_delta_Q_t']  = np.trapz(prop['delta_dT_t'], r)
        
            prop['blade_delta_CT_a'] = prop['blade_delta_T_a']/norm_CT
            prop['blade_delta_CQ_a'] = prop['blade_delta_Q_a']/norm_CQ
            prop['blade_delta_CT_t'] = prop['blade_delta_T_t']/norm_CT
            prop['blade_delta_CQ_t'] = prop['blade_delta_Q_t']/norm_CQ

        prop['blade_CT']        = prop['blade_delta_CT']+prop_iso['blade_CT']
        prop['blade_CQ']        = prop['blade_delta_CQ']+prop_iso['blade_CQ']
        prop['blade_eta']       = prop['blade_CT']/(2*np.pi*prop['blade_CQ'])*J_iso
        prop['blade_delta_eta'] = prop['blade_eta']-prop_iso['blade_eta']
        
        #calculate total forces
        prop['integral_delta_T']    = np.mean(prop['blade_delta_T'][:-1])*B
        prop['integral_delta_CT']   = prop['integral_delta_T']/norm_CT
        prop['integral_T']          = prop['integral_delta_T']+prop_iso['T']
        prop['integral_CT']         = prop['integral_delta_CT']+prop_iso['CT']
        
        prop['integral_delta_Q']    = np.mean(prop['blade_delta_Q'][:-1])*B
        prop['integral_delta_CQ']   = prop['integral_delta_Q']/norm_CQ
        prop['integral_Q']          = prop['integral_delta_Q']+prop_iso['Q']
        prop['integral_CQ']         = prop['integral_delta_CQ']+prop_iso['CQ']
        
        prop['integral_CP']         = prop['integral_CQ']*2*np.pi
        
        prop['integral_eta']        = prop['integral_CT']/(prop['integral_CP'])*J_iso
        prop['integral_delta_eta']  = prop['integral_eta']-prop_iso['eta']
        
        prop['integral_delta_Fy']   = np.mean(prop['blade_delta_Fy'][:-1])*B
        prop['integral_delta_Fz']   = np.mean(prop['blade_delta_Fz'][:-1])*B
        prop['integral_delta_CY']   = prop['integral_delta_Fy']/norm_CT
        prop['integral_delta_CZ']   = prop['integral_delta_Fz']/norm_CT
        
        #determine the RMS value for the unsteady case
        if prop_type=='us':
            #RMS is determined for three variables
            rms_dct = {'delta_T': 'blade_delta_T',
                       'delta_Fy': 'blade_delta_Fy',
                       'delta_Fz': 'blade_delta_Fz'}
            
            #divide based on the number of blades, with each 'blade' section
            #consisting of 12 subsections or timesteps
            dphi_rms = 2*np.pi/B/12
            phi_rms = np.arange(0, 2*np.pi, dphi_rms)
            
            #initiaze dictionary
            rms = {}
            rms['phi'] = phi_rms[:12]
            
            for key in rms_dct.keys():
                #interpolate results to phi_rms
                interp = si.pchip(phi, prop[rms_dct[key]])(phi_rms)
                #reshape so the blades are sorted by timestep
                interp = interp.reshape((B, 12))
                #sum the forces per timestep
                forces = np.sum(interp, axis=0)
                rms[key+'_forces'] = forces
                #calculate RMS
                rms[key] = np.std(forces-np.mean(forces))
                
            prop['rms'] = rms
            
        return prop
        
    def unsteady_correction(self):
        """Corrects quasi-steady loading results for unsteady effects using
        Sears' function.
        Original function written by Tomas Sinnige in MATLAB
        
        Returns:
            prop_us (dict): contains all the calculated propeller forces
        """
        
        #get propeller properties
        B = self.B
        n = self.n
        Dp = self.Dp
        Vinf = self.Vinf
        
        #get isolated propeller and quasi steady propeller
        prop_iso = self.prop_iso
        prop_qs = self.prop_qs
        
        #radius and chord distribution in [m]
        r = self.rR*Dp/2
        c = self.cR*Dp/2
        
        #define edges of blade segments
        r_edges = 0.5*(r[1:]+r[:-1])
        r_edges = np.append(self.Rh, r_edges)
        r_edges = np.append(r_edges, Dp/2)
        #get the spanwise extent of the blade sections
        dr = np.diff(r_edges)
        #first section is off due to varying hub radius along blade chord 
        #-> set equal to size of last segment 
        dr[0] = dr[-1]
        
        #compute induction and flow angle for this solution
        faxISO = B*prop_iso['dT']*dr        #estimated thrust force
        ftgtISO = B*prop_iso['dQ']*dr/r     #estimated torque 
        
        #compute induction factors isolated prop
        omega = n*2*np.pi
        aISO, bISO = actuator_disk(faxISO, ftgtISO, 
                                   self.rho_inf, Vinf, r, dr, Dp/2, omega)

        #define number of fourier coefficients
        phi = prop_qs['phi'][:-1]
        Ncoeff = phi.size
              
        #get corresponding wave numbers
        if (Ncoeff%2)==0: #even
            k_vGn1 = np.arange(0,Ncoeff/2+1,1)
            k_vGn2 = np.arange(-Ncoeff/2+1, 0, 1)
            k_vGn = np.append(k_vGn1, k_vGn2)
        else: #odd
            k_vGn1 = np.arange(0,np.floor(Ncoeff/2)+1,1)
            k_vGn2 = np.arange(-np.floor(Ncoeff/2), 0, 1)
            k_vGn = np.append(k_vGn1, k_vGn2)
          
        #compute value of sears function for all harmonics
        k_vGn_ = k_vGn.reshape((k_vGn.size, 1)) #change to column array
        Va = Vinf*(1+aISO)              #axial velocity including induced effects [m/s]
        Vt = omega*r*(1-bISO)           #tangential velocity including induced effects [m/s]
        Vheli = (Va**2+Vt**2)**0.5      #helicoidal velocities [m/s]
        Mheli = Vheli/self.a            #helicoidal Mach numbers
        ktheta = k_vGn_*omega/Vheli     #circumferential wave numbers
        sigma = ktheta*c/2              #reduced circumferential wave numbers
        
        #compute Sears function for input reduced frequencies and helicoidal
        #Mach numbers
        __,S_M = sears_function(sigma, Mheli)

        #Fourier transform the loading time histories (thrust and torque)
        delta_dTqs = prop_qs['delta_dT'][:-1, :] #Quasi steady thrust, skip the last angle
        delta_dQqs = prop_qs['delta_dQ'][:-1, :] #Quasi steady torque, skip the last angle
        
        #Fourier coefficients of thrust and torque deltas
        dTss_Fcoeff = np.fft.fft(delta_dTqs, axis=0)/Ncoeff
        dQss_Fcoeff = np.fft.fft(delta_dQqs/r, axis=0)/Ncoeff
        
        #unsteady results for thrust and torque
        #real part taken to get rid of tiny imaginary part due to numerical
        #error in Fourier transforms
        delta_dTus = np.real(np.fft.ifft(dTss_Fcoeff*S_M, axis=0)*Ncoeff)
        delta_dQus = np.real(np.fft.ifft(dQss_Fcoeff*S_M, axis=0)*Ncoeff)*r

#        #NOT USED        
#        dQss_Fcoeffalt = np.fft.fft(delta_dQqs, axis=0)/Ncoeff
#        delta_dQusalt = np.real(np.fft.ifft(dQss_Fcoeffalt*S_M, axis=0))*Ncoeff
        
        #unsteady thrust and torque
        delta_dT = np.vstack((delta_dTus, delta_dTus[0,:]))
        delta_dQ = np.vstack((delta_dQus, delta_dQus[0,:]))

        #calculate propeller forces
        prop_us = self.calculate_propeller_forces(delta_dT, delta_dQ)
        
        return prop_us

    def calculate_circulation(self, Va_i, Vt_i, rR_i):    
        """Calculates the circulation on the blade for the isolated propeller
        and the quasi-steady and unsteady results. Results are saved in the
        existing prop_iso, prop_qs and prop_us dictionaries.
        
        Args:
            Va_i (numpy.array): induced axial velocity at the propeller [m/s]
            Vt_i (numpy.array): induced tangential velocity at the propeller [m/s]
            rR_i (numpy.array): locations where the induced velocity is specified
        """
        
        #get propeller properties
        Dp = self.Dp
        rR = self.rR
        Vinf = self.Vinf
        rho_inf = self.rho_inf
        n = self.n
        
        #get analysis results
        prop_iso = self.prop_iso
        prop_qs = self.prop_qs
        prop_us = self.prop_us
        
        r = rR*Dp/2             #radii [m]
        
        #interpolate induced velocity values
        Va = si.pchip(rR_i, Va_i)(rR)
        Vt = si.pchip(rR_i, Vt_i)(rR)
        
        #calculate inducation factors
        a = Va/Vinf
        ap = Vt/(n*r*2*np.pi)
        
        #velocities at propeller disk
        Va = Vinf*(1+a)
        Vt = r*n*2*np.pi*(1-ap)
        V = (Va**2+Vt**2)**0.5
        prop_iso['V'] = V
        
        #blade pitch
        phi = np.arctan2(Va, Vt)
        prop_iso['phi_in'] = phi
        
        #lift and drag
        prop_iso['dL'] = prop_iso['dT']*np.cos(phi)+prop_iso['dQ']/r*np.sin(phi)
        prop_iso['dD'] = -prop_iso['dT']*np.sin(phi)+prop_iso['dQ']/r*np.cos(phi)
        
        #circulation
        prop_iso['gamma'] = prop_iso['dL']/(rho_inf*V)
        
        #velocities at propeller disk
        Va = prop_qs['u']*(1+a)
        Vt = (r*n*2*np.pi+prop_qs['Vt'])*(1-ap)
        V = (Va**2+Vt**2)**0.5
        prop_qs['V'] = V
        prop_us['V'] = V
        
        #blade pitch
        phi = np.arctan2(Va, Vt)
        prop_qs['phi_in'] = phi
        prop_us['phi_in'] = phi
        
        #lift and drag
        prop_qs['dL'] = prop_qs['dT']*np.cos(phi)+prop_qs['dQ']/r*np.sin(phi)
        prop_qs['dD'] = -prop_qs['dT']*np.sin(phi)+prop_qs['dQ']/r*np.cos(phi)
        prop_us['dL'] = prop_us['dT']*np.cos(phi)+prop_us['dQ']/r*np.sin(phi)
        prop_us['dD'] = -prop_us['dT']*np.sin(phi)+prop_us['dQ']/r*np.cos(phi)
        
        #circulation
        prop_qs['gamma'] = prop_qs['dL']/(rho_inf*V)
        prop_us['gamma'] = prop_us['dL']/(rho_inf*V)
        
        #save results
        self.prop_iso = prop_iso
        self.prop_qs = prop_qs
        self.prop_us = prop_us
    
    def init_slipstream(self, x, z, contraction=True):
        """Creates a slipstreamtube model for the current propeller results
        
        Args:
            x (numpy.array): x coordinates of the slipstream
            z (numpy.array): z coordinates of the slipstream
            contraction (bool): automatically calculates the contraction if True    
        """
        
        #get inputs
        st = self.slipstream
        st.B = self.B
        st.Vinf = self.Vinf
        st.omega = self.n*2*np.pi
        st.R = self.Dp/2
        st.r_hub = self.Rh
        st.sign_rotation = self.sign_rotation
        
        #create slipstream geometry
        r = self.prop_us['r']
        phi = self.prop_us['phi'][:-1]
        gamma = self.prop_us['gamma'][:-1]
        
        #create slipstreamtube model
        st.slipstream_geom(phi, r, x, z, gamma, contraction)
        self.slipstream = st
