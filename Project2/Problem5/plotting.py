#%%
import numpy as np
import matplotlib.pyplot as plt 
import tikzplotlib
import seaborn as sns

sns.set_theme()

data = np.loadtxt('countSimilarities.txt')

listOfN = data[0, :]
listOfIter = data[1, :]


fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (10, 10))
axs.plot(listOfN, listOfIter)
axs.set_xlabel(r'Size of matrix')
axs.set_ylabel(r'Number of iterations')
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
