#!/usr/bin/env python3


import numpy as np
from matplotlib  import pyplot as plt


data = np.loadtxt("cosmoTimeTable.dat", dtype=float, delimiter=",")

t = data[:,0]
a = data[:,1]

plt.figure()
plt.plot(t, a)
plt.xlabel("t")
plt.ylabel("a")
plt.show()
