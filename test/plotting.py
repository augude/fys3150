import numpy as np; import matplotlib.pyplot as plt 
import seaborn as sns 

sns.set_theme(font_scale = 2)

data = np.loadtxt("output.txt")

x = data[:, 0]
y = data[:, 1]

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (8, 8))
axs.scatter(x, y, s = 10, color = "red", label = "Discrete values")
axs.plot(x, y, linestyle = "-", color = "blue", label = "Fitted line")
axs.set_xlabel("x")
axs.set_ylabel(r"$e^x$")
fig.tight_layout()
plt.legend()
plt.savefig("plott1.pdf")

