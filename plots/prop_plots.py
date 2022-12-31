import niceplots
import numpy
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data/prop_geom.txt", sep=',')

niceplots.setRCParams()

_, (ax1,ax2) = plt.subplots(2,1,figsize=(7, 5))
ax1.plot(df['rR'], df['theta'])
ax1.set_ylabel(r'Twist, $\theta$')

ax2.plot(df['rR'], df['cR'])
ax2.set_xlabel(r'Propeller blade location, $x/b$')
ax2.set_ylabel(r'Chord/Radius, $c/R$')

niceplots.adjust_spines(ax1, outward=True)
niceplots.adjust_spines(ax2, outward=True)
plt.savefig('figures/PROWIM_twist.eps', format='eps')