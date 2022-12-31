import niceplots
import matplotlib.pyplot as plt

niceplots.setRCParams()
_, ax = plt.subplots(figsize=(10, 7))
ax.plot(span_, chord, label='Optimised')
ax.plot(span_, chord_orig_prop, label='Original')
ax.set_xlabel(r'Spanwise location $y$')
ax.set_ylabel(r'$m$')
ax.legend()
ax.grid()
niceplots.adjust_spines(ax, outward=True)
plt.savefig('00_results/figures/prop_results/chord_opt.png')