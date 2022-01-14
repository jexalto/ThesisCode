import matplotlib.pyplot as plt
import pickle
import sys
if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

folder = ''#'Convergence_study/'

file = folder + 'vlm_convstudy_grid3'
with open(file, 'rb') as vlm_file:
    vlm_conv_grid = pickle.load(vlm_file)

plt.figure(1)
plt.grid()
plt.figure(2)
plt.grid()
for i in vlm_conv_grid:
    plt.figure(1)
    vlm = vlm_conv_grid[i]
    plt.scatter(i, vlm.time, color='k')
    plt.xlabel('Propeller grid size[-]')
    plt.ylabel('Time [min]')
    plt.savefig('Convergence_study/grid_time')

    plt.figure(2)
    plt.scatter(i, vlm.res['CL'], color='k')
    plt.xlabel('Propeller grid size[-]')
    plt.ylabel('CL [-]')
    plt.savefig('Convergence_study/grid_cl')


file = folder + 'vlm_convstudy_N_x3'
with open(file, 'rb') as vlm_file:
    vlm_conv_N_x = pickle.load(vlm_file)

plt.figure(3)
plt.grid()
plt.figure(4)
plt.grid()
for i in vlm_conv_N_x:

    vlm = vlm_conv_N_x[i]
    plt.figure(3)
    plt.scatter(i, vlm.time, color='k')
    plt.xlabel('No. Slipstream x positions[-]')
    plt.ylabel('Time [min]')
    plt.savefig('Convergence_study/st_time')

    plt.figure(4)
    plt.scatter(i, vlm.res['CL'], color='k')
    plt.xlabel('No. Slipstream x positions[-]')
    plt.ylabel('CL [-]')
    plt.savefig('Convergence_study/st_cl')


file = folder + 'vlm_convstudy_strip3'
with open(file, 'rb') as vlm_file:
    vlm_conv_strip = pickle.load(vlm_file)

plt.figure(5)
plt.grid()
plt.figure(6)
plt.grid()
for i in vlm_conv_strip:
    plt.figure(5)
    vlm = vlm_conv_strip[i]
    plt.scatter(i, vlm.time, color='k')
    plt.xlabel('No. Strips as factor of original[-]')
    plt.ylabel('Time [min]')
    plt.savefig('Convergence_study/strip_time')

    plt.figure(6)
    plt.scatter(i, vlm.res['CL'], color='k')
    plt.xlabel('No. Strips as factor of original[-]')
    plt.ylabel('CL [-]')
    plt.savefig('Convergence_study/strip_cl')


file = folder + 'vlm_convstudy_stripx3'
with open(file, 'rb') as vlm_file:
    vlm_conv_stripx = pickle.load(vlm_file)

plt.figure(7)
plt.grid()
plt.figure(8)
plt.grid()
for i in vlm_conv_stripx:
    plt.figure(7)
    vlm = vlm_conv_stripx[i]
    plt.scatter(i, vlm.time, color='k')
    plt.xlabel('No. Strips as factor of original[-]')
    plt.ylabel('Time [min]')
    plt.savefig('Convergence_study/panel_time')

    plt.figure(8)
    plt.scatter(i, vlm.res['CL'], color='k')
    plt.xlabel('No. chordwise panels as factor of original[-]')
    plt.ylabel('CL [-]')
    plt.savefig('Convergence_study/panel_cl')



keys = list(vlm_conv_stripx.keys())
vlm0 = vlm_conv_stripx[keys[-1]]
for i in keys:
    vlm1 = vlm_conv_stripx[i]
    plt.figure(10)
    plt.scatter(i, abs((vlm1.res['CL'] - vlm0.res['CL']) / vlm0.res['CL']), color='k')
    plt.xlabel('No. chordwise panels as factor of original[-]')
    plt.ylabel(r'$\Delta$CL/ CL [-]')
plt.grid()
plt.savefig('Convergence_study/panel_dcl')


keys = list(vlm_conv_strip.keys())
vlm0 = vlm_conv_strip[max(keys)]
for i in keys:

    vlm1 = vlm_conv_strip[i]
    plt.figure(11)
    plt.scatter(i, abs((vlm1.res['CL'] - vlm0.res['CL']) / vlm0.res['CL']), color='k')
    plt.xlabel('No. spanwise  strips as factor of original[-]')
    plt.ylabel(r'$\Delta$CL/ CL [-]')
plt.grid()
plt.savefig('Convergence_study/strip_dcl')


keys = list(vlm_conv_grid.keys())
vlm0 = vlm_conv_grid[max(keys)]
for i in keys:

    vlm1 = vlm_conv_grid[i]
    plt.figure(12)
    plt.scatter(i, abs((vlm1.res['CL'] - vlm0.res['CL']) / vlm0.res['CL']), color='k')
    plt.xlabel('Propeller grid size [-]')
    plt.ylabel(r'$\Delta$CL/ CL [-]')
plt.grid()
plt.savefig('Convergence_study/grid_dcl')


keys = list(vlm_conv_N_x.keys())
vlm0 = vlm_conv_N_x[max(keys)]
plt.figure(13)
for i in keys:
    vlm1 = vlm_conv_N_x[i]
    plt.scatter(i, abs((vlm1.res['CL'] - vlm0.res['CL']) / vlm0.res['CL']), color='k')
    plt.xlabel('No. Slipstream x positions [-]')
    plt.ylabel(r'$\Delta$CL/ CL [-]')
plt.grid()
plt.savefig('Convergence_study/st_dcl')