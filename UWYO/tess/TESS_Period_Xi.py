import numpy as np
import matplotlib.pyplot as plt
from pylab import *

ptimes, mags, merrs = np.loadtxt('/d/users/brock/TESS/xi_Aql_b.csv', skiprows=1,unpack=True,delimiter=',',usecols=[0,3,4])
period = 1.99472

phase = (np.mod((ptimes-2459000),period))/period
plt.plot(ptimes, mags, 'b.', linestyle='')
plt.show()
