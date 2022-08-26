import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

sns.set_theme(font_scale = 2)

data = np.loadtxt('exactPlot.txt')

x = data[0, :]
u = data[1, :]

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (8, 8))
axs.plot(x, u)
axs.set_xlabel('x')
axs.set_ylabel('u(x)')
fig.tight_layout()
plt.savefig('exactPlot.pdf')
