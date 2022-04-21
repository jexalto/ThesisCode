import pandas as pd
import numpy as np

data = pd.read_csv('/home/jexalto/code/MDO_lab_env/ThesisCode/HELIX_verification/Ilinois_data/da4052_9x6.75_geom.txt', dtype=np.float16, sep='  ', skiprows=(1), header=None).values

rR = data[:, 0]
cR = data[:, 1]
beta = data[:, 2]

print(rR, cR, beta)
