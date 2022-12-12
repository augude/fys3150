import numpy as np; import matplotlib.pyplot as plt; import seaborn as sns; import pyarma as pa
sns.set_theme(font_scale = 2)
from matplotlib.animation import FuncAnimation

def create_animation(type):
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)

    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    x, y = np.meshgrid(x_points, y_points, sparse=True)

    # Array of time points
    dt = 2.5e-5
    t_points = np.arange(0, 1+dt, dt)

    # Some settings
    fontsize = 12
    t_min = t_points[0]
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.absolute(obj[0, :, :])**2))

    # Plot the first frame
    img = ax.imshow(np.absolute(np.transpose(obj[0, :, :]))**2, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm, origin = 'lower')

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.grid(0)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label(r"Probability", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.absolute(obj[i, :, :])**2))
        img.set_norm(norm)

        # Update z data
        img.set_data(np.transpose(np.absolute(obj[i, :, :])**2))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=50, frames=np.arange(0, len(obj[:, 0, 0]), 2), repeat=False, blit=0)

    # Run the animation!
    #plt.show()

    anim.save(f'../output/{type}.gif', writer="Pillow", bitrate=-1, fps=30)
    
def deviating_probability(type, type1, T):
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)

    N = len(obj[:, 0, 0])
    t = np.linspace(0, T, N)
    
    prob = np.array([np.sum(np.absolute(obj[i, :, :])**2) for i in range(N)])

    obj1 = pa.cx_cube()
    obj1.load(f"../output/{type1}_evolution.bin")
    obj1 = np.array(obj1)

    
    prob1 = np.array([np.sum(np.absolute(obj1[i, :, :])**2) for i in range(N)])

    
    fig, axs = plt.subplots(1, 2, figsize = (15, 7))
    axs[0].plot(t, np.abs(1 - prob))
    axs[0].set_xlabel('t')
    axs[0].set_ylabel(r'$|1 - \sum_{i, j}|u_{i, j}|^2|$')
    axs[0].set_title('Without potential barriers', y = 1.1)
    axs[0].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    
    axs[1].plot(t, np.abs(1 - prob1))
    axs[1].set_xlabel('t')
    axs[1].set_ylabel(r'$\log(|1 - \sum_{i, j}|u_{i, j}|^2|)$')
    axs[1].set_title('With double-slit barriers', y = 1.1)
    axs[1].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    
    fig.tight_layout()    
    plt.savefig(f'../output/{type}_{type1}_deviation.pdf', bbox_inches = 'tight')
    plt.show() 
       
    
def measure(type, y, t, T, title = None):
    from scipy.signal import argrelextrema
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)
    
    M = len(obj[0, 0, :])
    N = len(obj[:, 0, 0]) 
    y_ind = int(y*M)
    t_ind = int(t/T)*N - 1
    x = np.linspace(0, 1, M)
    prob = np.absolute(obj[t_ind, y_ind, :])**2
    
    max_max = argrelextrema(prob, np.greater)
    print(x[max_max])
    fig, axs = plt.subplots(1, 1, figsize = (8, 8))
    axs.scatter(x, prob/np.sum(prob))
    axs.plot(x, prob/np.sum(prob), alpha = 0.5)
    print(np.trapz(prob/np.sum(prob)))
    axs.set_xlabel('y')
    axs.set_ylabel(f'$p(y|x = {y:.2f}; t = {t})$')
    if title is not None:
        fig.suptitle(f'{title}')
    else:
        fig.suptitle(f'Detection probability \n after passing {type}-slit barrier')
    fig.tight_layout()    
    plt.savefig(f'../output/{type}_detection.pdf', bbox_inches = 'tight')
    plt.show()
    
def evolution(type, T, plotstyle):
    func = {'abs': lambda x: np.absolute(x)**2, 'real': lambda x: np.real(x), 'imag': lambda x: np.imag(x)}
    name = {'abs': 'Probability', 'real': 'Real part', 'imag': 'Imaginary part'}
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)

    h = 0.005
    x_points = np.arange(0, 1+h, h)
    y_points = np.arange(0, 1+h, h)
    
    x_min, x_max = x_points[0], x_points[-1]
    y_min, y_max = y_points[0], y_points[-1]

    # Create a colour scale normalization according to the max z value in the first frame
    
    N = len(obj[:, 0, 0]) 
        
    fig, axs = plt.subplots(1, 3, figsize = (25, 7))
    t_inds = [0, N//2, N - 1]
    for index, t_ind in enumerate(t_inds):
        t = t_ind*T/N
        norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(func[plotstyle](obj[t_ind, :, :])))
        img = axs[index].imshow(func[plotstyle](np.transpose(obj[t_ind, :, :])), extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm, origin = 'lower')
        axs[index].text(0.95, 0.95, "t = {:.0e}".format(t), color="white", 
                            horizontalalignment="right", verticalalignment="top")
        axs[index].grid(0)
        axs[index].set_xlabel('x')
        axs[index].set_ylabel('y')
    
        # Add a colourbar
        cbar = fig.colorbar(img, ax=axs[index], fraction=0.046, pad=0.02, format = '%.0e')
        cbar.set_label(name[plotstyle], rotation = -90)
        cbar.ax.tick_params()
    fig.suptitle(r'Evolution of $u(x, y)$ at different points in time')
    fig.tight_layout()
    plt.savefig(f'../output/{type}_{plotstyle}.pdf', bbox_inches = 'tight')
    plt.show()
    