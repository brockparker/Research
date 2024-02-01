import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

data = pd.read_csv (r'C:\Users\Brock\oceanview\HR4D30811__9__00009.txt')
df = pd.DataFrame(data, columns = ['Current', 'Voltage'])

def gauss(x, a, b, m, s):
    return a * np.exp( -((x-m)**2) / (2 * s**2) ) + b

    
    
popt, pcov = curve_fit(gauss, curr[:-1], voltage[:-1])
cerr = np.sqrt(np.diag(pcov))
    