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
fig, axs = plt.subplots(ncols = 2, nrows = 1, figsize = (10, 10), sharey = True)

for index, file in enumerate(files):
    data = np.loadtxt(file)
    listOfN = data[0, :]
    listOfIter = data[1, :]
    
    params, pcov = curve_fit(fit, listOfN, listOfIter)
    print(params)
    fitted = fit(listOfN, params[0], params[1], params[2], params[3])
    
    axs[index].plot(listOfN, listOfIter)
    axs[index].plot(listOfN, fitted, linestyle = '--')
    axs[index].set_xlabel(r'N')
    axs[index].legend(['I(N)', r'$I_{fit} = \mathcal{O}(N^2)$'])
axs[0].set_ylabel(r'I(N)')
axs[0].text(x = 0, y = 16500, s = 'Tridiagonal matrix')
axs[1].text(x = 0, y = 16500, s = 'Dense matrix')
tikzplotlib.clean_figure()
tikzplotlib.save(
    f"{file[:-4]}.tex",
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