import numpy as np
import matplotlib.pyplot as plt
from pylab import *

ptimes, mags, merrs = np.loadtxt('/d/users/brock/TESS/TOI_1581_ASASSN.csv', skiprows=1,unpack=True,delimiter=',',usecols=[0,3,4])
period = 3.63434

phase = (np.mod((ptimes-2459000),period))/period
plt.plot(phase, mags, 'b.', linestyle='')
plt.show()
