import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import copy
import time
    
if 'BEM/' not in sys.path:
    sys.path.append('BEM/')
    
if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

sys.path.insert(1, '/home/jexalto/code/MDO_lab_env/ThesisCode/WingAerodynamics/Aero/VLM')
sys.path.insert(1, '/home/jexalto/code/MDO_lab_env/ThesisCode/WingAerodynamics/Aero/BEM')
    
from vlm_or import PyVLM    # original non-adapted vlm
from BEM import BEM
from TipMountedProp import PropWing

bem_dct = {}
vlm_dct = {}
vlm_dct2 = {}
conv_dct = {}

### ISOLATED PROP
dir_prop_wing = '/home/jexalto/code/MDO_lab_env/ThesisCode/WingAerodynamics/Aero/Validation/SecIIIC_ModelII/'

fig24 = pd.read_csv(dir_prop_wing+'SecIIIC_Fig24_CLCD_ModelII_tipMounted_ReD640k.txt',
                   header=20)
fig24 = fig24.drop(0)
fig24 = fig24.drop(columns=['config'])
fig24 = fig24.astype('float')

fig24_val = pd.DataFrame(columns=['Polar', 'AoA', 'CL', 'CD'], index=np.arange(0, 3*10, 1))
aoa = np.linspace(-8, 16, 5)

c_r = 0.24
c_t = 0.24
b = 0.730*2*0.952
#b = 0.310*2
R = 0.1184
rR_hub = 0.148 #0.0175/R
mu = 1.8121e-05
D = 2*R
B = 4
rR_beta = 0.75
pitch = 23.9-1.2 #deg
sign_rotation = 1 #1 = outboard up, -1 = inboard up
prop_geom = pd.read_csv('/home/jexalto/code/MDO_lab_env/ThesisCode/WingAerodynamics/Aero/Validation/SecIIIC_ModelII/prop_geom.txt')

x = 0.853*R
z = 0.042*R


le_1 = np.array([0,0])
le_3 = np.array([(c_r-c_t)/2,b/2])

le_2 = le_1+np.array([(le_3[0]-le_1[0])/(le_3[1]-le_1[1])*(le_3[1]-R),
                      le_3[1]-R])
c_2 = c_r+(c_t-c_r)/(le_3[1]-le_1[1])*(le_3[1]-R)

le = [le_1, le_2, le_3]
chord = [c_r, c_2, c_t]

prop_loc = np.array([-x, b/2, z]) #prop_loc = np.array([-x-c_t, b/2, z])
prop_ang = 0

n = [1,1] #chordwise
#n = [8,8]
m = [10,10] #spanwise
opt = ['nsin', 'uni']

# le = [le_1, le_3]
# chord = [c_r, c_t]
# n = [8]
# m = [10]
vlm = PyVLM()

airfoil = np.loadtxt('prop_airfoils/NACA64(2)-A015.dat', skiprows=1)
vlm.airfoil.input_selig(airfoil[:, 0], airfoil[:, 1], 'NACA64(2)-A015')

vlm.add_wing(le, chord, n, m, opt)
vlm.polar_dir = 'data_aero_wing/val/'

bem = BEM()
        
bem.R = R
bem.B = B
bem.beta = pitch
        
theta = prop_geom['theta'].to_numpy()
rR = prop_geom['rR'].to_numpy()
theta0 = np.interp(rR_beta, rR, theta)
theta = theta-theta0

bem.theta_b = theta
bem.cR_b = prop_geom['cR'].to_numpy()
bem.rR_b = rR
bem.rR_hub = rR_hub
bem.airfoil_folder = 'data_aero_prop/val/'

ind = 0
#for p in [0, 2, 4]:
#aoa = np.linspace(-12, 25, 20)
#aoa = np.linspace(-10, 34, 20)
#aoa = [5]

for p in [4]:
#for p in [0, 2, 4]:
    
    if p==0:
        for alpha in aoa:
            vlm.Vinf = fig24[fig24['polar']==p+1]['Vinf'].mean()
            vlm.rho = fig24[fig24['polar']==p+1]['rhoInf'].mean()
            vlm.a = (1.4*287*fig24[fig24['polar']==p+1]['Tinf'].mean())**0.5
            vlm.mu = mu
            
            vlm.alpha = alpha
            vlm.Vpert['turb'] = True
            
            vlm.vlm_visc()
            
#            vlm.vlm()
#            vlm.strip_properties()
            
            fig24_val.loc[ind]['CL'] = vlm.res['CL']
            fig24_val.loc[ind]['CD'] = vlm.res['CD']
            fig24_val.loc[ind]['Polar'] = p+1
            fig24_val.loc[ind]['AoA'] = alpha
            
            vlm_dct2[alpha] = copy.deepcopy(vlm)
            
            ind+=1
    else:
        
        vlm.Vinf = fig24[fig24['polar']==p+1]['Vinf'].mean()
        vlm.rho = fig24[fig24['polar']==p+1]['rhoInf'].mean()
        vlm.a = (1.4*287*fig24[fig24['polar']==p+1]['Tinf'].mean())**0.5
        vlm.mu = mu
        
        bem.Vinf = fig24[fig24['polar']==p+1]['Vinf'].mean()
        bem.omega = fig24[fig24['polar']==p+1]['n'].mean()*2*np.pi
        bem.rho = fig24[fig24['polar']==p+1]['rhoInf'].mean()
        bem.a  = (1.4*287*fig24[fig24['polar']==p+1]['Tinf'].mean())**0.5
        bem.mu = mu
    
        propwing = PropWing(vlm, bem, prop_loc)
        propwing.N_iter_max = 4
#        propwing.dJ = 0.2
    
        #for alpha in aoa:
        for alpha in [2]:
            
            # propwing.settings = {'dJ':                  0.3,
            #                      'N_prop_map':          3,
            #                      'slipstream_length':   2,
            #                      'dlbda':               0.1/0.8,
            #                      'dlb':                 0.005,
            #                      'p_max':               int(10*0.8),
            #                      'lbda_max':            5*0.8,
            #                      'N_time_step':         int(30*0.8),
            #                      'N_slipstream_x':      int(12*0.8),
            #                      'conway_s_max':        500*0.8,
            #                      'conway_ds':           0.1/0.8}

            propwing.settings = {'dJ': 0.5,
                             'N_prop_map': 3,
                             'slipstream_length': 2,
                             'dlbda': 0.1,
                             'dlb': 0.005,
                             'p_max': 10,
                             'lbda_max': 5,
                             'N_time_step': 30,
                             'N_slipstream_x': 12,
                             'conway_s_max': 500,
                             'conway_ds': 0.1}
            t0 = time.time()
            propwing.analysis(alpha)
            print(t0-time.time())
            
            prop_us = propwing.urot.prop_us
            vlm = propwing.vlm
            
            L_t = prop_us['integral_T']*np.sin(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.cos(np.deg2rad(alpha))
            D_t = -prop_us['integral_T']*np.cos(np.deg2rad(alpha))+prop_us['integral_delta_Fz']*np.sin(np.deg2rad(alpha))
            
            qS = 0.5*vlm.res['rho']*vlm.res['Vinf']**2*vlm.res['S']
            
            fig24_val.loc[ind]['CL'] = vlm.res['CL']+2*L_t/qS
            fig24_val.loc[ind]['CD'] = vlm.res['CD']+2*D_t/qS
            fig24_val.loc[ind]['Polar'] = p+1
            fig24_val.loc[ind]['AoA'] = alpha
            
            vlm_dct[alpha] = copy.deepcopy(vlm)
            bem_dct[alpha] = copy.deepcopy(prop_us)
            conv_dct[alpha] = copy.deepcopy(propwing.conv)
            
            ind+=1

#data = pd.read_csv('Validation/metal_wing_cl.csv', header=None, names=['a', 'cl'])
#data1 = pd.read_csv('Validation/metal_wing_cd.csv', header=None, names=['cd', 'cl'])
#data1 = data1.sort_values(by=['cl'])
#data2 = data1[data1['cl']>0.9]
#data2 = data2.sort_values(['cd'])
#data1 = data1[data1['cl']<0.9]
#data1 = data1.append(data2)
#
##data3 = pd.read_csv('Validation/yb2=0.21, a=2.csv', header=None, names=['xc', 'cp'])
##
##vlm = vlm_dct2[2]
#    
#
#
#plt.figure()
#plt.plot(data['a'], data['cl'])
#plt.plot(fig24_val1[fig24_val1['Polar']==1]['AoA'], fig24_val1[fig24_val1['Polar']==1]['CL'])
#plt.plot(fig24_val[fig24_val['Polar']==1]['AoA'], fig24_val[fig24_val['Polar']==1]['CL'])
##plt.legend(['data Koomen', 'VLM no thickness corr', 'VLM with thickness corr'])
#plt.legend(['data Koomen', 'VLM no viscous corr', 'VLM with viscous corr'])
##plt.legend(['data Koomen', 'VLM'])
#plt.xlabel(r'$\alpha$ [deg]')
#plt.ylabel('$C_L$ [-]')
#plt.grid()
#plt.savefig('no_nacelle_cl.png', dpi=300)
#
#plt.figure()
#plt.plot(data1['cd'], data1['cl'])
#plt.plot(fig24_val1[fig24_val1['Polar']==1]['CD'], fig24_val1[fig24_val1['Polar']==1]['CL'])
#plt.plot(fig24_val[fig24_val['Polar']==1]['CD'], fig24_val[fig24_val['Polar']==1]['CL'])
##plt.legend(['data Koomen', 'VLM no thickness corr', 'VLM with thickness corr'])
#plt.legend(['data Koomen', 'VLM no viscous corr', 'VLM with viscous corr'])
##plt.legend(['data Koomen', 'VLM'])
#plt.xlabel('$C_D$ [-]')
#plt.ylabel('$C_L$ [-]')
#plt.grid()
#plt.savefig('no_nacelle_cd.png', dpi=300)

#plt.figure()
#plt.plot(data3['xc'], data3['cp'], 'x')
#plt.plot(np.hstack((x_, x_))/c_r, 1-vi**2, 'x')
#plt.gcap.invert_yaxis()
#plt.grid()
#plt.xlabel('$x/c$ [-]')
#plt.ylabel('$C_p$ [-]')


plt.figure()
plt.plot(fig24[fig24['polar']==1]['AoA'], fig24[fig24['polar']==1]['CL'])       
plt.plot(fig24[fig24['polar']==3]['AoA'], fig24[fig24['polar']==3]['CL'])       
plt.plot(fig24[fig24['polar']==5]['AoA'], fig24[fig24['polar']==5]['CL'])       
plt.plot(fig24_val[fig24_val['Polar']==1]['AoA'], fig24_val[fig24_val['Polar']==1]['CL'])    
plt.plot(fig24_val[fig24_val['Polar']==3]['AoA'], fig24_val[fig24_val['Polar']==3]['CL'])    
plt.plot(fig24_val[fig24_val['Polar']==5]['AoA'], fig24_val[fig24_val['Polar']==5]['CL'])
plt.legend(['prop-off Sinnige et al.', 'J=0.9 Sinnige et al.', 'J=0.7 Sinnige et al.', 'prop-off', 'J=0.9', 'J=0.7']) 
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('$C_L$ [-]')
plt.grid()
plt.savefig('fig24_cl_p.png', dpi=300)

plt.figure()
plt.plot(fig24[fig24['polar']==1]['AoA'], fig24[fig24['polar']==1]['CD'])       
plt.plot(fig24[fig24['polar']==5]['AoA'], fig24[fig24['polar']==3]['CD']) 
plt.plot(fig24[fig24['polar']==5]['AoA'], fig24[fig24['polar']==5]['CD'])   
plt.plot(fig24_val[fig24_val['Polar']==1]['AoA'], fig24_val[fig24_val['Polar']==1]['CD'])           
plt.plot(fig24_val[fig24_val['Polar']==3]['AoA'], fig24_val[fig24_val['Polar']==3]['CD'])           
plt.plot(fig24_val[fig24_val['Polar']==5]['AoA'], fig24_val[fig24_val['Polar']==5]['CD']) 
plt.legend(['prop-off Sinnige et al.', 'J=0.9 Sinnige et al.', 'J=0.7 Sinnige et al.', 'prop-off', 'J=0.9', 'J=0.7']) 
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel('$C_D$ [-]')
plt.grid()
plt.savefig('fig24_cd_p.png', dpi=300)

plt.figure()
plt.plot(fig24[fig24['polar']==1]['CD'], fig24[fig24['polar']==1]['CL'])       
plt.plot(fig24[fig24['polar']==3]['CD'], fig24[fig24['polar']==3]['CL'])  
plt.plot(fig24[fig24['polar']==5]['CD'], fig24[fig24['polar']==5]['CL'])   
plt.plot(fig24_val[fig24_val['Polar']==1]['CD'], fig24_val[fig24_val['Polar']==1]['CL'])           
plt.plot(fig24_val[fig24_val['Polar']==3]['CD'], fig24_val[fig24_val['Polar']==3]['CL'])           
plt.plot(fig24_val[fig24_val['Polar']==5]['CD'], fig24_val[fig24_val['Polar']==5]['CL']) 
plt.legend(['prop-off Sinnige et al.', 'J=0.9 Sinnige et al.', 'J=0.7 Sinnige et al.', 'prop-off', 'J=0.9', 'J=0.7']) 
plt.xlabel('$C_D$ [-]')
plt.ylabel('$C_L$ [-]')
plt.grid()
plt.savefig('fig24_clcd_p.png', dpi=300)


plt.figure()
plt.plot(fig24[fig24['polar']==3]['AoA'], fig24[fig24['polar']==1]['CL'].values-fig24[fig24['polar']==3]['CL'])       
plt.plot(fig24[fig24['polar']==5]['AoA'], fig24[fig24['polar']==1]['CL'].values-fig24[fig24['polar']==5]['CL'])       
plt.plot(fig24_val[fig24_val['Polar']==3]['AoA'], fig24_val[fig24_val['Polar']==1]['CL'].values-fig24_val[fig24_val['Polar']==3]['CL'])    
plt.plot(fig24_val[fig24_val['Polar']==5]['AoA'], fig24_val[fig24_val['Polar']==1]['CL'].values-fig24_val[fig24_val['Polar']==5]['CL'])
plt.legend(['J=0.9 Sinnige et al.', 'J=0.7 Sinnige et al.', 'J=0.9', 'J=0.7']) 
plt.xlabel(r'$\alpha$ [deg]')
plt.ylabel(r'$\Delta C_L$ [-]')
plt.grid()
plt.savefig('fig24_cl_diff_p.png', dpi=300)

            
            
#plt.figure()
#plt.plot(fig24[fig24['polar']==1]['AoA'], fig24[fig24['polar']==1]['CL'])       
#plt.plot(fig24_val1[fig24_val1['Polar']==1]['AoA'], fig24_val1[fig24_val1['Polar']==1]['CL'])    
#plt.plot(fig24_val2[fig24_val2['Polar']==1]['AoA'], fig24_val2[fig24_val2['Polar']==1]['CL'])   
#plt.plot(fig24_val3[fig24_val3['Polar']==1]['AoA'], fig24_val3[fig24_val3['Polar']==1]['CL'])    
#plt.plot(fig24_val4[fig24_val4['Polar']==1]['AoA'], fig24_val4[fig24_val4['Polar']==1]['CL'])    
#plt.legend(['prop-off Sinnige et al.', 'VLM b=1 tc=1', 'VLM b=0.952 tc=1', 'VLM b=1 tc=1.11', 'VLM b=0.952 tc=1.11'])  
#plt.xlabel(r'$\alpha$ [deg]')
#plt.ylabel('$C_L$ [-]')
#plt.grid()
#plt.savefig('fig24_cl_p.png', dpi=300)
#
#plt.figure()
#plt.plot(fig24[fig24['polar']==1]['AoA'], fig24[fig24['polar']==1]['CD'])       
#plt.plot(fig24_val1[fig24_val1['Polar']==1]['AoA'], fig24_val1[fig24_val1['Polar']==1]['CD'])   
#plt.plot(fig24_val2[fig24_val2['Polar']==1]['AoA'], fig24_val2[fig24_val2['Polar']==1]['CD'])
#plt.plot(fig24_val3[fig24_val3['Polar']==1]['AoA'], fig24_val3[fig24_val3['Polar']==1]['CD']) 
#plt.plot(fig24_val4[fig24_val4['Polar']==1]['AoA'], fig24_val4[fig24_val4['Polar']==1]['CD']) 
#plt.legend(['prop-off Sinnige et al.', 'VLM b=1 tc=1', 'VLM b=0.952 tc=1', 'VLM b=1 tc=1.11', 'VLM b=0.952 tc=1.11'])  
#plt.xlabel(r'$\alpha$ [deg]')
#plt.ylabel('$C_D$ [-]')
#plt.grid()
#plt.savefig('fig24_cd_p.png', dpi=300)
#
#plt.figure()
#plt.plot(fig24[fig24['polar']==1]['CD'], fig24[fig24['polar']==1]['CL'])       
#plt.plot(fig24_val1[fig24_val1['Polar']==1]['CD'], fig24_val1[fig24_val1['Polar']==1]['CL'])  
#plt.plot(fig24_val2[fig24_val2['Polar']==1]['CD'], fig24_val2[fig24_val2['Polar']==1]['CL'])  
#plt.plot(fig24_val3[fig24_val3['Polar']==1]['CD'], fig24_val3[fig24_val3['Polar']==1]['CL'])           
#plt.plot(fig24_val4[fig24_val4['Polar']==1]['CD'], fig24_val4[fig24_val4['Polar']==1]['CL'])   
#plt.legend(['prop-off Sinnige et al.', 'VLM b=1 tc=1', 'VLM b=0.952 tc=1', 'VLM b=1 tc=1.11', 'VLM b=0.952 tc=1.11'])  
#plt.xlabel('$C_D$ [-]')
#plt.ylabel('$C_L$ [-]')
#plt.grid()
#plt.savefig('fig24_clcd_p.png', dpi=300)

#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['gamma'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['gamma'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('gamma')
#plt.grid()
#
#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['cl'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['cl'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('cl')
#plt.grid()
#
#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['cdi'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['cdi'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('cdi')
#plt.grid()
#
#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['cdp'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['cdp'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('cdp')
#plt.grid()
#
#
#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['alpha_eff'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['alpha_eff'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('alpha_eff [deg]')
#plt.grid()
#
#plt.figure()
#plt.plot(vlm_dct[alpha].Strip.index/b*2, vlm_dct[alpha].Strip['V_ind'], vlm_dct2[alpha].Strip.index/b*2, vlm_dct2[alpha].Strip['V_ind'])
#plt.legend(['with prop', 'without prop'])
#plt.xlabel('y [m]')
#plt.ylabel('V_ind [m/s]')
#plt.grid()