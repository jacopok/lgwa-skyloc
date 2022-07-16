import numpy as np
import matplotlib.pyplot as plt

f = np.logspace(-4, -2, num=1000)

fn = 1e-3

Q = 200

plt.semilogx(f, np.angle(-f**2 / (fn**2 - f**2 + 1j * fn**2 / Q)))

plt.show()