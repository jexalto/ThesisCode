import matplotlib.pyplot as plt
import pickle
import sys
if 'VLM/' not in sys.path:
    sys.path.append('VLM/')

if 'Convergence_study/' not in sys.path:
    sys.path.append('Convergence_study/')

folder = 'Convergence_study/'
case = 'm1_5'

file = folder + 'convergence_' + case
with open(file, 'rb') as conv_file:
    results = pickle.load(conv_file)

plt.figure(1, figsize=(6.5,2.5))
plt.grid()
plt.figure(2, figsize=(6.5,2.5))
plt.grid()
plt.figure(3, figsize=(6.5,2.5))
plt.grid()
plt.subplots_adjust(bottom=0.18)

for i in results:
    plt.figure(1)
    res = results[i]
    plt.scatter(i, res['time'], color='k')

    plt.figure(2)
    plt.scatter(i, res['vlm_res']['CL'], color='k')

    plt.figure(3)
    endi = max(results.keys())
    cl = results[endi]['vlm_res']['CL']
    dcl = abs(res['vlm_res']['CL']-cl)
    plt.scatter(i, dcl/cl, color='k')

if 'n' in case:
    a = 'chordwise'
elif 'm'  in case :
    a = 'spanwise'

plt.figure(1)
plt.ylabel('Time [min]')
plt.xlabel('Factor of ' + a + ' panels[-]')
plt.savefig('Convergence_study/Figures/' + case + '_time')
plt.subplots_adjust(bottom=0.18)


plt.figure(2)
plt.xlabel('Factor of ' + a + ' panels[-]')
plt.ylabel('CL [-]')
plt.savefig('Convergence_study/Figures/' + case + '_cl')
plt.subplots_adjust(bottom=0.18)

plt.figure(3)
plt.xlabel('Factor of ' + a + ' panels[-]')
plt.ylabel(r'$\Delta$CL/CL [-]')
plt.savefig('Convergence_study/Figures/' + case + '_dcl')
plt.subplots_adjust(bottom=0.18)