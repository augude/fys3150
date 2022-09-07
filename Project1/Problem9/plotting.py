#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

sns.set_theme(font_scale = 2)

dataExact = np.loadtxt('/home/agude/Project1/fys3150/Project1/Problem2/exactPlot.txt')
xExact = dataExact[0, :]
uExact = dataExact[1, :]

data10 = np.loadtxt('tridiagAlgo100.txt')
x10 = data10[0, :]
u10 = data10[1, :]

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
axs.plot(x10, u10)
axs.plot(xExact, uExact, linestyle = 'dashed')
axs.text(x = 0.0, y = 0.6, s = 'a')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$v(x)$')

fig.tight_layout()
plt.savefig('tridiagAlgo.pdf')
plt.show()
# %%
