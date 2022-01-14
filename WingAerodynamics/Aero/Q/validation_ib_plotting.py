import pickle
import numpy as  np
import pandas as pd
import matplotlib.pyplot as plt



dir_prop_wing = 'Validation/SecIIIC_ModelII/'

fig24 = pd.read_csv(dir_prop_wing+'SecIIIC_Fig24_CLCD_ModelII_conventional_ReD640k.txt',
                   header=20)
fig24 = fig24.drop(0)
fig24 = fig24.drop(columns=['config'])
fig24 = fig24.astype('float')


ps = [1,2,3,4,5]     # polars to plot


def plot_polars(ps, config, folder, save=False, conv=True):
    save_dir = 'Validation/Results/' + folder +'/'

    if config == 'wt':
        name = 'tipMounted'
    elif config =='ib':
        name = 'conventional'
    fig24 = pd.read_csv(dir_prop_wing + 'SecIIIC_Fig24_CLCD_ModelII_' + name + '_ReD640k.txt',
                        header=20)
    fig24 = fig24.drop(0)
    fig24 = fig24.drop(columns=['config'])
    fig24 = fig24.astype('float')

    with open('Validation/Results/'+ folder+ '/fig24_val_data_' + config, 'rb') as fig24_val_file:
        fig24_val = pickle.load(fig24_val_file)
    Js = ['inf', '1.0', '0.9', '0.8', '0.7']
    color = ['b', 'r', 'g', 'c', 'm']
    leg = []
    plt.figure(1)
    for p in ps:
        col = color[p-1]
        J = Js[p-1]
        plt.plot(fig24[fig24['polar']==p]['AoA'], fig24[fig24['polar']==p]['CL'], color=col)
        plt.plot(fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['AoA'], fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['CL'], '--', color=col)
        leg.append('J='+ J + ' Sinnige et al')
        leg.append('J='+ J)

    plt.legend(leg)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_L$ [-]')
    plt.xlim(-8,16)
    plt.ylim(-0.5,1.5)
    plt.grid()
    if save:
        plt.savefig(save_dir + 'fig24_cl_p' + str(ps) + '.png')

    plt.figure(2)
    leg = []
    for p in ps:
        col = color[p-1]
        J = Js[p-1]
        plt.plot(fig24[fig24['polar']==p]['AoA'], fig24[fig24['polar']==p]['CD'], color=col)
        plt.plot(fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['AoA'], fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['CD'], '--', color=col)
        leg.append('J='+ J + ' Sinnige et al')
        leg.append('J='+ J)

    plt.legend(leg)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel('$C_D$ [-]')
    plt.grid()
    if save:
        plt.savefig(save_dir +'fig24_cd_p' + str(ps) + '.png')

    plt.figure(3)
    leg = []
    for p in ps:
        col = color[p - 1]
        J = Js[p-1]
        plt.plot(fig24[fig24['polar']==p]['CD'], fig24[fig24['polar']==p]['CL'], color=col)
        plt.plot(fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['CD'], fig24_val[fig24_val['Polar']==p][fig24_val['conv']==True]['CL'], '--', color=col)
        leg.append('J='+ J + ' Sinnige et al')
        leg.append('J='+ J)

    plt.legend(leg)
    plt.xlabel('$C_D$ [deg]')
    plt.ylabel('$C_L$ [-]')
    plt.grid()
    if save:
        plt.savefig(save_dir +'fig24_cl_cd_p' + str(ps) + '.png')
        plt.close('all')
    return fig24_val


def combine_df():
    config = ''
    with open('Validation/Results/fig24_data_' + config + '1', 'rb') as fig24_val_file:
       fig24_val0 = pickle.load(fig24_val_file)

    with open('Validation/Results/fig24_data_' + '2', 'rb') as fig24_val_file:
        fig24_val1 = pickle.load(fig24_val_file)

    with open('Validation/Results/fig24_data_' + config + '3', 'rb') as fig24_val_file:
        fig24_val2 = pickle.load(fig24_val_file)

    with open('Validation/Results/fig24_data_'+ config + '4', 'rb') as fig24_val_file:
        fig24_val3 = pickle.load(fig24_val_file)

    with open('Validation/Results/fig24_data_'+ config + '5', 'rb') as fig24_val_file:
        fig24_val4 = pickle.load(fig24_val_file)


    df = pd.DataFrame(columns=['Polar', 'AoA', 'CL', 'CD', 'conv'], index=np.arange(0, 100, 1))
    ind = 0
    for Polar in np.arange(1,6):
        if not np.isnan(Polar):
            for fig_val in [ fig24_val0, fig24_val1, fig24_val2, fig24_val3, fig24_val4]: #, fig24_val3, fig24_val4]:

                data = fig_val[fig_val['Polar'] == Polar]
                for aoa in data['AoA']:

                    data_ = data[data['AoA'] == aoa]
                    i = data_['AoA'].keys()[0]
                    df.loc[ind]['CL'] = data_['CL'][i]
                    df.loc[ind]['CD'] = data_['CD'][i]
                    df.loc[ind]['Polar'] = data_['Polar'][i]
                    df.loc[ind]['AoA'] = data_['AoA'][i]
                    df.loc[ind]['conv'] = data_['conv'][i]
                    ind += 1
    with open('Validation/Results/fig24_val_data_ib' + config, 'wb') as fig_data_file:
        pickle.dump(df, fig_data_file)

    return df

import sys
if 'VLM/' not in sys.path:
        sys.path.append('VLM/')


# with open('Validation/Results/vlm_dct_' + '3', 'rb') as fig24_val_file:
#     vlm_dct = pickle.load(fig24_val_file)
#
#
# vlm_data= vlm_dct[2]
# #vlm2 = vlm_dct[2][6]
#
# for key in  [0,2,4,6,8,10]:#vlm_data.keys():
#         vlm = vlm_data[key]
#         plt.plot(vlm.Strip.index, vlm.Strip['cl'])


