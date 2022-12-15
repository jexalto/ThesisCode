import numpy as np
import pandas as pd

def camber(x_input):
    data = pd.read_csv('/home/jexalto99/code/MDO_lab_env/ThesisCode/validation/wingprop_validation/naca62415.csv')

    x_upper = data['x'][:26]-0.5
    y_upper = data['y'][:26]
    x_lower = data['x'][26:]-0.5
    y_lower = data['y'][26:]

    y_upper_output = np.interp(x_input, x_upper, y_upper)
    y_lower_output = np.interp(x_input, x_lower, y_lower)

    camber = (y_upper_output-y_lower_output)/2

    return camber