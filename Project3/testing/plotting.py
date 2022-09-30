#%%

import numpy as np; import matplotlib.pyplot as plt 
import seaborn as sns
import tikzplotlib

sns.set_theme(font_scale = 2)

data = np.loadtxt('ElectricField.txt')

N = 100
d = 1e4
Ex = np.zeros(shape = (N, N))
Ey = np.zeros(shape = (N, N))
Ez = np.zeros(shape = (N, N))

X = np.zeros(shape = (N, N))
Y = np.zeros(shape = (N, N))
Z = np.zeros(shape = (N, N))

V = np.zeros(shape = (N, N))

for i in range(N**2):
    ind1 = int(data[i, 0])
    ind2 = int(data[i, 1])
    X[ind1, ind2] = data[i, 2]
    Y[ind1, ind2] = data[i, 3]
    Z[ind1, ind2] = data[i, 4] 

    Ex[ind1, ind2] = data[i, 5]
    Ey[ind1, ind2] = data[i, 6]
    Ez[ind1, ind2] = data[i, 7]  
    V[ind1, ind2] = data[i, 8]       
    
fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (12, 10))
plt.contourf(Y, Z, V, cmap = 'cool')
cbar = plt.colorbar(cmap = 'cool')
cbar.set_label(r'Electric potential [$\frac{u(\mu m)^2}{(\mu s)^2 e}$]')
axs.quiver(Y[::5, ::5], Z[::5, ::5], Ey[::5, ::5], Ez[::5, ::5])
axs.set_xlabel(r'y [$\mu m$]')
axs.set_ylabel(r'z [$\mu m$]')
fig.tight_layout()
plt.savefig('PenningTrapSetup.pdf')
plt.show()

# %%
