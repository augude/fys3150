#%%
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns; sns.set_theme(font_scale=2)
import sys


def write_potential(M, v0, xwall_thickness, xwall_center, aperature, slit_seperation, type, save = True):
    ind_center = int((M - 2)*xwall_center)
    ind_thick = int((M - 2)*xwall_thickness)
    ind_aperature = int((M - 2)*aperature)
    ind_slit = int((M - 2)*slit_seperation)
    
    if type == 'no':
        pot = np.zeros(shape = (M - 2, M-2)) 
        
    elif type == 'single':
        pot = np.zeros(shape = (M - 2, M-2)) 
        for i in range(ind_center - ind_thick, ind_center + ind_thick + 1):
            for j in range(M - 2):
                if 0 < j < ((M - 2)//2 - ind_aperature//2):
                    pot[i, j] = v0
        pot = np.flip(pot, axis = 1) + pot
        
    elif type == 'double':        
        pot = np.zeros(shape = (M - 2, M-2)) 
        for i in range(ind_center - ind_thick, ind_center + ind_thick + 1):
            for j in range(M - 2):
                if (0 < j < ((M - 2)//2 - ind_slit//2 - ind_aperature) 
                or ( (M - 1)//2 - ind_slit//2 < j < (M - 1)//2 + ind_slit//2 ) 
                or  ((M - 2)//2 + ind_slit//2 + ind_aperature) < j < M):
                    pot[i, j] = v0
                          
    elif type == 'triple':        
        pot = np.zeros(shape = (M - 2, M-2)) 
        for i in range(ind_center - ind_thick, ind_center + ind_thick + 1):
            for j in range(M - 2):
                if (0 < j < ((M - 2)//2 - ind_slit//2 - ind_aperature - ind_slit) 
                or ( ((M - 2)//2 - ind_aperature//2 - ind_slit) < j < ((M - 2)//2 - ind_aperature//2))):
                #or  ((M - 2)//2 + ind_slit//2 + ind_aperature) < j < M):
                    pot[i, j] = v0
        pot = np.flip(pot, axis = 1) + pot
    
    elif type == 'tunnel':
        pot = np.zeros(shape = (M - 2, M-2)) 
        for i in range(ind_center - ind_thick, ind_center + ind_thick + 1):
            for j in range(M - 2):
                pot[i, j] = v0
    else:
        print('Provide a valid argument. Either "no", "single", "double", "triple" or "tunnel".')
        return None 
    if save:
        np.savetxt(f'{type}.dat', pot)
    return pot 

h = float(sys.argv[2])
M = int(1.0/h + 1)
v0 = float(sys.argv[3])
thick = 0.02
center = 0.5
slit_sep = 0.05
aperature = 0.05 

typ = str(sys.argv[1])
pot = write_potential(M, v0, thick, center, aperature, slit_sep, type = typ)

if pot is not None:
    print(f'A grid of size {np.shape(pot)} with {typ}-slit barrier is successfully created.')
    plt.imshow(pot, origin = 'lower')
    plt.grid()

            
# %%