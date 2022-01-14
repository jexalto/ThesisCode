import numpy as np
import sys

if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

from Xfoil import Xfoil

#%% wing


""" Can be used to generate wing polars"""
modes = ['turbulent']
for mode in modes:
    for M in np.arange(0.5, 0.55, 0.05):
        for Re in np.arange(1e7, 1.25e7, 2.5e6):
            xf = Xfoil()
            xf.run_dir = 'prop_airfoils'
            xf.airfoil = 'NACA663-418'
            #xf.NACA = '2412'
            xf.visc = True
            #        xf.N_panel = 160
            xf.Re = Re
            xf.M = M
            if mode is 'turbulent':
               xf.N = 9 #0.01
               xf.vacc = 0
               xf.x_trip_bot = 0.05
               xf.x_trip_up = 0.05

            elif mode is 'laminar':
                xf.N = 9
                xf.vacc = 0
                xf.x_trip_bot = 1
                xf.x_trip_up = 1

            result_u = xf.create_polar(0, 20, 0.1)
            result_l = xf.create_polar(0, -20, 0.1)

            if result_u is None or result_l is None:
                print(Re, M)
            else:
                result_l = result_l.iloc[1:]
                result_l = result_l.iloc[::-1]

                result = result_l.append(result_u, ignore_index=True)

                name = xf.airfoil + '_Re%d_M%0.2f_N%d' % (xf.Re, xf.M, xf.N*100)

                name = name+'_xTRU%0.2f_xTRL%0.2f' % (xf.x_trip_up, xf.x_trip_bot)
                name = name+'.txt'

                f = open('data_aero_wing/ATR72/' + mode + '/'+name, 'w')
                #f = open('new_polars/' + mode + '/' + name, 'w')

                f.write(xf.airfoil + '\n')
                f.write('Re %12.5g\n' % xf.Re)
                f.write('M %12.5g\n' % xf.M)
                f.write('N %12.5g\n' % xf.N)
                f.write('VACC %12.5g VACCflag %d\n' % (xf.vacc, 1))
                f.write('XTR %12.5g %12.5g XTRflag %d\n' % (xf.x_trip_bot, xf.x_trip_up, 1))
                f.write('alpha,cl,cd,cdp,cm\n')
                for i in result.index:
                    row = (result.loc[i, 'alpha'],
                           result.loc[i, 'CL'],
                           result.loc[i, 'CD'],
                           result.loc[i, 'CDp'],
                           result.loc[i, 'CM'])
                    f.write('%12.5g,%12.5g,%12.5g,%12.5g,%12.5g\n' % row)

                f.close()
