import matplotlib.pyplot as plt
import niceplots
import numpy as np
import scipy.interpolate as interpolate

conway_14 = [1, 1.012, 1.019, 1.022, 1.0218, 1.021, 1.02, 1.0185, 1.0175, 1.01675, 1.016, 1.015, 1.014, 1.0135, 1.012, 1.011]
conway_57 = [1., 0.997, 0.993, 0.99, 0.988, 0.986, 0.984, 0.9825, 0.981, 0.98, 0.979, 0.978, 0.977, 0.976, 0.9755, 0.975]
conway_83 = [1, 0.985, 0.975, 0.971, 0.968, 0.967, 0.9655, 0.96425, 0.96325, 0.9625, 0.96175, 0.961, 0.9605, 0.96, 0.9595, 0.9595]
veldhuis = [1., 0.995, 0.989, 0.985, 0.982, 0.98, 0.978, 0.976, 0.975, 0.974, 0.973, 0.972, 0.971, 0.97075, 0.9705, 0.97025]

radial = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625, 1.875, 2. ]

radial_linspace = np.linspace(0, 2, 100)

f_conway_14 = interpolate.UnivariateSpline(radial, conway_14, k=5)
f_conway_57 = interpolate.UnivariateSpline(radial, conway_57, k=5)
f_conway_83 = interpolate.UnivariateSpline(radial, conway_83, k=5)
f_veldhuis = interpolate.UnivariateSpline(radial, veldhuis, k=5)

def readdata(filename):
    with open(filename, 'r') as f:
    # Initialize an empty list
        data = []
        radial = []
        # Iterate through each line in the file
        for line in f:
            # Split the line into two columns
            col1, col2 = line.split(sep=',')
            col2 = col2.strip('\n')
            radial.append(float(col1))
            data.append(float(col2))
        # Add the two columns to the data list
    return radial, data

radial_14, conway_14 = readdata('data/conway_14.txt')
radial_57, conway_57 = readdata('data/conway_57.txt')
radial_83, conway_83 = readdata('data/conway_83.txt')
radial_veldhuis, veldhuis = readdata('data/veldhuis.txt')

niceplots.setRCParams()

plt, ax = plt.subplots(figsize=(10, 7))
niceplots.adjust_spines(ax, outward=True)
ax.plot(radial_14, conway_14, label=r'Conway $r/R=0.14$')
ax.plot(radial_57, conway_57, label=r'Conway $r/R=0.57$')
ax.plot(radial_83, conway_83, label=r'Conway $r/R=0.83$')
ax.plot(radial_veldhuis, veldhuis, label=r'Veldhuis')
ax.set_xlim(0, 2.5)
ax.set_xlabel(r'$x/R$')
ax.set_ylabel(r'$r/r_0$')
ax.legend()
plt.savefig('figures/conway.eps', fomrat='eps')
plt.savefig('figures/conway.png')