#%%
import numpy as np
import matplotlib.pyplot as plt 
import tikzplotlib
import seaborn as sns

sns.set_theme()

for N in [9, 99]:
    
    n = N + 1
    
    valFile = f"eigenvalues{N}.txt"
    vecFile = f"eigenvectors{N}.txt"
    val = np.loadtxt(valFile) #the vectors and values are already sorted
    vec = np.loadtxt(vecFile)
        
    x = np.arange(0, 1.0 + 1.0/n, 1.0/n)
    
    fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (10, 10))
    for i in range(3):
        v = np.concatenate((np.zeros(1), vec[i, :], np.zeros(1))) #approximation
        v = -1*v*(v[1] < 0) + 1*v*(v[1] >= 0)
        u = np.concatenate((np.zeros(1), vec[N + i, :], np.zeros(1))) #exact
        axs.plot(x, v, color = 'blue')
        axs.plot(x, u, color = 'red', linestyle = '--')
    axs.set_xlabel(r'x [L]')
    axs.set_ylabel(r'v(x)')
    plt.legend(['Approximation', 'Exact'])
    tikzplotlib.clean_figure()
    tikzplotlib.save(
                f"solve{N + 1}.tex",
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
