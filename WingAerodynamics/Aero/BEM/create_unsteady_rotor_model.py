import numpy as np
import copy

from UnsteadyRotor import UnsteadyRotor


def create_unsteady_rotor_model(bem_, dJ, N):
    """Creates an unsteady rotor model from a BEM instance. The BEM instance
    is used to create a propeller performance map for varying J with Vinf and
    RPS.
    
    Args:
        bem_ (BEM instance): instance of BEM class, conditions are taking from here
        dJ (float): the range of J in the propeller map [J-dJ, J+dJ]
        N (int): the number of points used to create the propeller map
            N=2 means 2 points for each dJ, so this leads to 8 evaluations plus
            1 evaluation at the center
    
    Returns:
        urot (UnsteadyRotor instance): UnsteadyRotor with the same conditions
            as the BEM instance and the propeller map as input
    """
    
    #create a copy
    bem = copy.deepcopy(bem_)
    
    if bem.res is None:
        bem.analysis()
    
    #variables from original analysis
    Vinf0 = bem.Vinf
    omega0 = bem.omega
    n0 = bem.omega/(2*np.pi)
    J0 = bem.res_sum['J']
    D = 2*bem.R
    
    #calculate the additional advance ratios and corresponding velocity and
    #rotational speeds
    J_range_ = np.linspace(J0-dJ, J0+dJ, 2*N+1)
    J_range = np.delete(np.linspace(J0-dJ, J0+dJ, 2*N+1), N)
    V_range = J_range*n0*D
    omega_range = Vinf0/(J_range*D)*2*np.pi
    
    #arrays to save results
    constant_V_dCT = np.zeros((2*N+1, bem.res['dCT'].size))
    constant_V_dCQ = np.zeros((2*N+1, bem.res['dCT'].size))
    constant_RPS_dCT = np.zeros((2*N+1, bem.res['dCT'].size))
    constant_RPS_dCQ = np.zeros((2*N+1, bem.res['dCT'].size))
    
    #save results for original advace ratio
    constant_V_dCT[N, :] = bem.res['dCT']
    constant_V_dCQ[N, :] = bem.res['dCQ']
    constant_RPS_dCT[N, :] = bem.res['dCT']
    constant_RPS_dCQ[N, :] = bem.res['dCQ']
    
    #calculation for other advance ratios
    for i,J in enumerate(J_range):
        if i<N:
            i_ = i
        else:
            i_ = i+1
        
        for j in range(2): #for constant V and constant RPS
            if j==0:
                bem.omega = omega_range[i]
                bem.Vinf = Vinf0
            else:
                bem.omega = omega0
                bem.Vinf = V_range[i]
                
            #analyse propeller   
            bem.analysis()
            
            #put original advance ratio data if bem does not converge
            if not bem.res_sum['converged']:
                bem.res['dCT'] = constant_V_dCT[N, :]
                bem.res['dCQ'] = constant_V_dCQ[N, :]
            
            #save results
            if j==0:
                constant_V_dCT[i_, :] = bem.res['dCT']
                constant_V_dCQ[i_, :] = bem.res['dCQ']    
            else:
                constant_RPS_dCT[i_, :] = bem.res['dCT']
                constant_RPS_dCQ[i_, :] = bem.res['dCQ']
    
    #create unsteady model with the same input varibles
    urot = UnsteadyRotor()  
    urot.Dp = bem.R*2               #propeller diameter [m]
    urot.Rh = bem.rR_hub*bem.R      #propeller hub radius [m]
    urot.B = bem.B                  #propeller number of blades [-]
    urot.n =  n0                    #propeller rps [1/s]
    urot.Vinf = Vinf0               #freestream velocity [m/s]
    urot.rho_inf = bem.rho          #freestream air density [kg/m3]
    urot.a = bem.a                  #speed of sound [m/s]
    
    #radial distribution to be used
    urot.rR = bem.res['rR'].to_numpy()
    #define chord distribution
    urot.set_chord(bem.res['rR'], bem.res['cR'])
    
    #set propeller maps
    urot.constant_V.J        = J_range_
    urot.constant_V.rR       = bem.res['rR'].to_numpy()
    urot.constant_V.dCT      = constant_V_dCT
    urot.constant_V.dCQ      = constant_V_dCQ
    
    urot.constant_RPS.J      = J_range_
    urot.constant_RPS.rR     = bem.res['rR'].to_numpy()
    urot.constant_RPS.dCT    = constant_RPS_dCT
    urot.constant_RPS.dCQ    = constant_RPS_dCQ
    
    return urot
