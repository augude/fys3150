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
    general = data
    # special = data[1, :]
    axs.semilogy(general)
    # axs[1].plot(special)

axs.set_xlabel(r'$N$')
axs.set_ylabel('Time (log(s))')
lgd = axs.set_title("General Algorithm")
plt.show()
fig.savefig('timings_comparison.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
