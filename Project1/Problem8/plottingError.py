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
lgd = axs.legend(loc = 'upper center', ncol = len(files)//2, fancybox = True
           , bbox_to_anchor = (0.5, 1.2))
plt.show()
fig.savefig('errors.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

relErrorMax = np.zeros(len(files))
steps = np.zeros(len(files))
fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
for index, filename in enumerate(files):
    data = np.loadtxt(filename) #loads the approximations from txt file
    x = data[0, :][1:-1] #extracts the x-points and discards the boundary conditions
    v = data[1, :][1:-1] #the approximation
    u = exact(x) #calculates the exact solution 
    relerror = ( np.abs( (u - v)/u ) ) #relative error 
    relErrorMax[index] = np.max(relerror)
    steps[index] = (10**(index + 1))
    axs.semilogy(x, relerror, label = f'$n = {10**(index + 1)}$')
axs.set_xlabel(r'$x$')
axs.set_ylabel(r'$\log_{10}(\epsilon)$')
lgd = axs.legend(loc = 'upper center', ncol = len(files)//2, fancybox = True, 
           bbox_to_anchor = (0.5, 1.25))
fig.savefig('relerrors.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

#%%
print(steps)
print(relErrorMax)

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))
axs.scatter(steps, relErrorMax, s = 40)
axs.plot(steps, relErrorMax, linestyle = "--", alpha = 0.5)
axs.set_yscale('log')
axs.set_xscale('log')
axs.set_xlabel(r'$n_{steps}$')
axs.set_xlabel(r'$\log_{10} n_{steps}$')
axs.set_ylabel(r'$\log_{10} \epsilon$')
fig.savefig('relerrorsmax.pdf', bbox_inches='tight')
plt.show()
