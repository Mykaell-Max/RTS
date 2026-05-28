import numpy as np
import matplotlib.pyplot as plt

PI = np.pi
x = np.linspace(0,1,33)
exact_temp = 400.0 + 1400.0*(np.sin(PI*(x/1.0)))**2.0

plt.plot(x, exact_temp)