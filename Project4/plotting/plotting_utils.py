import pyarma as pa
import numpy as np
import pandas as pd

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
    heatCapacity = 1/(L**2*T**T)*(energy2mom - energy1mom**2)
    mag1mom = np.cumsum(abs(mag))/(cycles + 1)
    mag2mom = np.cumsum(mag**2)/(cycles + 1)
    susceptibility = 1/(L**2*T)*(mag2mom - mag1mom**2)
    df = pd.DataFrame(
            {
                'energy': energy,
                'magnetization': mag,
                'energy1mom': energy1mom,
                'energy2mom': energy2mom,
                'mag1mom': mag1mom,
                'mag2mom': mag2mom,
                'heatCapacity': heatCapacity,
                'susceptibility': susceptibility,
                'temperature': temperature,
                'gridsize': gridsize
        }
    )
    
    return df
