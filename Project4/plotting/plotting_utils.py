import pyarma as pa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def binToDf(filename):
    """ returns a pandas dataframe  

    Args:
        filename (_string_): bin-file with energy and magnetization for each cycle
    """
    data = pa.mat()
    data.load(filename)
    data = np.array(data)
    L = data[0,0]
    T = data[0,1]
    energy = data[1:,0]
    mag = data[1:,1]
    cycles = np.arange(0, len(energy))
    temperature = np.zeros(len(energy))
    temperature[:] = T
    gridsize = np.zeros(len(energy))
    gridsize[:] = L
    energy1mom = np.cumsum(energy)/(cycles + 1)
    energy2mom = np.cumsum(energy**2)/(cycles + 1)
    heatCapacity = 1/(L**2*T*T)*(energy2mom - energy1mom**2)
    mag1mom = np.cumsum(abs(mag))/(cycles + 1)
    mag2mom = np.cumsum(mag**2)/(cycles + 1)
    susceptibility = 1/(L**2*T)*(mag2mom - mag1mom**2)
    df = pd.DataFrame(
            {
                'energy': energy,
                'magnetization': mag,
                'energy1mom': energy1mom,
                'energy2mom': energy2mom,
                'magnetization1mom': mag1mom,
                'magnetization2mom': mag2mom,
                'heatCapacity': heatCapacity,
                'susceptibility': susceptibility,
                'temperature': temperature,
                'gridsize': gridsize
        }
    )
    
    return df

def extract_data(files):
    summary = pd.DataFrame(columns=['L', 'T', 'E', 'M', 'C', 'X'])
    eng = np.zeros(len(files))
    mag = np.zeros(len(files))
    cv = np.zeros(len(files))
    chi = np.zeros(len(files))
    temp = np.zeros(len(files))
    grid = np.zeros(len(files))
    for index, file in enumerate(files):
        df = binToDf(file)
        temp[index] = df.temperature[0] 
        grid[index] = df.gridsize[0]
        eng[index] = df.energy[int(0.2*1e6):].mean()/df.gridsize[0]**2
        eng2 = np.mean(df.energy[int(0.2*1e6):]**2)/df.gridsize[0]**4
        mag[index] = np.mean(abs(df.magnetization[int(0.2*1e6):]))/df.gridsize[0]**2
        mag2 = np.mean(df.magnetization[int(0.2*1e6):]**2)/df.gridsize[0]**4
        cv[index] = 1/(temp[index]*temp[index])*(eng2 - eng[index]**2)*df.gridsize[0]**2
        chi[index] = 1/(temp[index])*(mag2 - mag[index]**2)*df.gridsize[0]**2
    
    summary['L'] = grid
    summary['T'] = temp
    summary['E'] = eng
    summary['M'] = mag
    summary['C'] = cv
    summary['X'] = chi
    return summary

def scatter_plot(xs, ys, xlabel, ylabel, labels, savefig=False, filename='plot.png'):
    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    for index, x in enumerate(xs):
        axs.plot(x, ys[index], linestyle='--')
        axs.scatter(x, ys[index], label=labels[index])

    axs.set_ylabel(ylabel)
    axs.set_xlabel(xlabel)
    lgd = fig.legend(loc = 'lower center', ncol = len(xs), fancybox = True, 
                bbox_to_anchor = (0.5, -0.05))
    fig.tight_layout()
    plt.show()
    if savefig:
        fig.savefig(filename)