#%%
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
#%%

sns.set_theme(font_scale = 2)

files = ['/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo10.txt', 
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo100.txt',
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo1000.txt',
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo10000.txt',
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo100000.txt',
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo1000000.txt',
         '/home/agude/Project1/fys3150/Project1/Problem7/tridiagAlgo10000000.txt'
         ]

def exact(x):
    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    return u

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
for index, filename in enumerate(files):
    data = np.loadtxt(filename) #loads the approximations from txt file
    x = data[0, :][1:-1] #extracts the x-points and discards the boundary conditions
    v = data[1, :][1:-1] #the approximation
    u = exact(x) #calculates the exact solution 
    error = ( np.abs( u - v) ) #error 
    axs.semilogy(x, error, label = f'$n = {10**(index + 1)}$')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$\log_{10}(\Delta)$')
fig.tight_layout()
axs.legend(loc = 'upper center', ncol = len(files)//2, fancybox = True
           , bbox_to_anchor = (0.5, 1.2))
plt.show()
plt.savefig('errors.pdf')
#%%
relErrorMax = np.zeros_like(files)
fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
for index, filename in enumerate(files):
    data = np.loadtxt(filename) #loads the approximations from txt file
    x = data[0, :][1:-1] #extracts the x-points and discards the boundary conditions
    v = data[1, :][1:-1] #the approximation
    u = exact(x) #calculates the exact solution 
    relerror = ( np.abs( (u - v)/u ) ) #reletive error 
    relErrorMax[index] = np.max(relerror)
    axs.semilogy(x, relerror, label = f'$n = {10**(index + 1)}$')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$\log_{10}(\epsilon)$')
axs.legend(loc = 'upper center', ncol = len(files)//2, fancybox = True, 
           bbox_to_anchor = (0.5, 1.25))
plt.show()
plt.savefig('relerrors.pdf')
# %%

print(relErrorMax)
#%%
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 10), sharey = True)
steps = 10**(np.arange(1, 8))
axs[0].scatter(steps, relErrorMax)
axs[1].semilogx(steps, relErrorMax)
axs[0].set_xlabel(r'$n_{steps}$')
axs[1].set_xlabel(r'$\log n_{steps}$')
axs[0].set_ylabel(r'$\epsilon$')
axs[0].set_yticks(np.linspace(10**(-6), 10**(-2), 10 ))
plt.show()
# %%

