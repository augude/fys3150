#%%
import numpy as np
import matplotlib.pyplot as plt 
import tikzplotlib
import seaborn as sns
from scipy.optimize import curve_fit

sns.set_theme()

files = np.array(['countSimilarities.txt', 'countSimilaritiesDense.txt'])
for file in files:
    data = np.loadtxt(file)
    listOfN = data[0, :]
    listOfIter = data[1, :]
    
    fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (10, 10))
    axs.plot(listOfN, listOfIter)
    axs.set_xlabel(r'N')
    axs.set_ylabel(r'I(N)')
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
def fit(x, a, b, c, d):
    return a*x**b + c*x**(b-1) + d

files = np.array(['countSimilarities.txt'])
for file in files:
    data = np.loadtxt(file)
    listOfN = data[0, :]
    listOfIter = data[1, :]
    
    params, pcov = curve_fit(fit, listOfN, listOfIter)
    print(params)
    fitted = fit(listOfN, params[0], params[1], params[2], params[3])

    fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (10, 10))
    axs.plot(listOfN, listOfIter, label = r'$I_(N)$')
    axs.plot(listOfN, fitted, label = f'O(N^{params[1]:.2f})',linestyle = '--')
    axs.set_xlabel(r'N') 
    axs.set_ylabel(r'I(N)')
    plt.legend()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
            f"fitted.tex",
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
