#%%
import numpy as np
import matplotlib.pyplot as plt 
import tikzplotlib
import seaborn as sns
from scipy.optimize import curve_fit

sns.set_theme()

def fit(x, a, b, c, d):
    return a*x**b + c*x**(b-1) + d


files = np.array(['countSimilarities.txt', 'countSimilaritiesDense.txt'])
fig, axs = plt.subplots(ncols = 1, nrows = 2, figsize = (10, 10), sharex = True)

for index, file in enumerate(files):
    data = np.loadtxt(file)
    listOfN = data[0, :]
    listOfIter = data[1, :]
    
    params, pcov = curve_fit(fit, listOfN, listOfIter)
    print(params)
    fitted = fit(listOfN, params[0], params[1], params[2], params[3])
    
    axs[index].plot(listOfN, listOfIter, color = 'blue', linewidth = 4.0)
    axs[index].plot(listOfN, fitted, linestyle = '--', color = 'red')
    axs[index].set_ylabel(r'I(N)')
    axs[index].legend(['I(N)', r'$I_{fit}(N) = \mathcal{O}(N^2)$'])
axs[0].set_xlabel(r'N')
axs[0].text(x = 80, y = 20, s = 'Tridiagonal matrix')
axs[1].text(x = 80, y = 20, s = 'Dense matrix')
tikzplotlib.clean_figure()
tikzplotlib.save(
    f"convergence.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
plt.show()

# %%