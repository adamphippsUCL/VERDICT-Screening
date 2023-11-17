import numpy as np
import matplotlib.pyplot as plt


reals = np.random.normal(loc = 0, scale = 1, size = 10000)
ims = np.random.normal(loc = 0, scale = 1, size = 10000)

signal = 0
signals = np.sqrt( ims**2 + (signal+reals)**2)


plt.figure()
plt.hist(signals, bins = 20)
plt.show()