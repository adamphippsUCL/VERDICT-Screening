import numpy as np
import matplotlib.pyplot as plt

signal = 1
NoiseMagnitdue = 0.1

reals = np.random.normal(loc = 0, scale = NoiseMagnitdue, size = 10000)
ims = np.random.normal(loc = 0, scale = NoiseMagnitdue, size = 10000)
signals = np.sqrt( ims**2 + (signal+reals)**2)

print(np.mean(signals))
plt.figure()
plt.hist(signals, bins = 50)
plt.show()