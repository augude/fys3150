#%%
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
#%%

sns.set_theme(font_scale = 2)

files = ['timings.txt']

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
for index, filename in enumerate(files):
    data = np.loadtxt(filename) #loads the approximations from txt file
    Ns = [10**i for i in range(1, 7)]
    general = data[0, :]
    special = data[1, :]
    axs.plot(Ns, general)
    axs.plot(Ns, special)
    axs.scatter(Ns, general, s = 40)
    axs.scatter(Ns, special, s = 40)
    axs.set_xscale('log')
    axs.set_yscale('log')

axs.set_xlabel(r'$N$')
axs.set_ylabel('Time (log(s))')
lgd = plt.legend(["General", "Special"])
plt.show()
fig.savefig('timings_comparison.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

# %%
