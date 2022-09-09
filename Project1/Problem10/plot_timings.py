#%%
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

sns.set_theme(font_scale = 2)

general = np.loadtxt('timingsGeneral.txt')
special = np.loadtxt('timingsSpecial.txt')
n = len(general[0, :])
generalMean = np.mean(general, axis = 1)
generalSTD = np.std(general, axis = 1, ddof = 1)

specialMean = np.mean(special, axis = 1)
specialSTD = np.std(special, axis = 1, ddof = 1)

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))
Ns = [10**i for i in range(1, 7)]
axs.plot(Ns, generalMean, linestyle = '--', alpha = 0.8, color = '#1f77b4')
axs.plot(Ns, specialMean, linestyle = '--', alpha = 0.8, color = 'orange')

axs.fill_between(Ns, generalMean - generalSTD, generalMean + generalSTD, alpha = 0.2, color = '#1f77b4')
axs.fill_between(Ns, specialMean - specialSTD, specialMean + specialSTD, alpha = 0.2, color = 'orange')

axs.scatter(Ns, generalMean, s = 45, c = '#1f77b4')
axs.scatter(Ns, specialMean, s = 45, c = 'orange')

axs.set_xscale('log')
axs.set_yscale('log')

axs.set_xlabel(r'$n$')
axs.set_ylabel('Time (log(s))')
lgd = plt.legend(["General", "Special"])
plt.show()
fig.savefig('timings_comparison.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

# %%
