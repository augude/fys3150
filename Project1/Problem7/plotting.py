#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

sns.set_theme(font_scale = 2)

dataExact = np.loadtxt('/home/agude/Project1/fys3150/Project1/Problem2/exactPlot.txt')
xExact = dataExact[0, :]
uExact = dataExact[1, :]

data10 = np.loadtxt('tridiagAlgo10.txt')
x10 = data10[0, :]
u10 = data10[1, :]

data100 = np.loadtxt('tridiagAlgo100.txt')
x100 = data100[0, :]
u100 = data100[1, :]

data1000 = np.loadtxt('tridiagAlgo1000.txt')
x1000 = data100[0, :]
u1000 = data100[1, :]

data10000 = np.loadtxt('tridiagAlgo10000.txt')
x10000 = data10000[0, :]
u10000 = data10000[1, :]


fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 8))
axs[0, 0].plot(x10, u10)
axs[0, 0].plot(xExact, uExact, linestyle = 'dashed')
axs[0, 0].set_xlabel(r'$x$')
axs[0, 0].set_ylabel(r'$v(x)$')

axs[1, 0].plot(x100, u100)
axs[1, 0].plot(xExact, uExact, linestyle = 'dashed')
axs[1, 0].set_xlabel(r'$x$')
axs[1, 0].set_ylabel(r'$v(x)$')

axs[0, 1].plot(x1000, u1000)
axs[0, 1].plot(xExact, uExact, linestyle = 'dashed')
axs[0, 1].set_xlabel(r'$x$')
axs[0, 1].set_ylabel(r'$v(x)$')

axs[1, 1].plot(x10000, u10000)
axs[1, 1].plot(xExact, uExact, linestyle = 'dashed')
axs[1, 1].set_xlabel(r'$x$')
axs[1, 1].set_ylabel(r'$v(x)$')

fig.tight_layout()
#plt.savefig('tridiagAlgo.pdf')
plt.show()
# %%
