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
    norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.absolute(obj[0, :, :])))

    # Plot the first frame
    img = ax.imshow(np.absolute(obj[0, :, :]), extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

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
        norm = plt.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.absolute(obj[i, :, :])))
        img.set_norm(norm)

        # Update z data
        img.set_data(np.absolute(obj[i, :, :]))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=5, frames=np.arange(0, len(obj[:, 0, 0]), 2), repeat=False, blit=0)

    # Run the animation!
    #plt.show()

    anim.save(f'./{type}.gif', writer="Pillow", bitrate=-1, fps=30)
    
def deviating_probability(type, T, dt):
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)

    N = int(T/dt - 1)
    t = np.linspace(0, T, N)
    prob = np.array([np.sum(np.absolute(obj[i, :, :])**2) for i in range(N)])

    fig, axs = plt.subplots(1, 1, figsize = (8, 8))
    axs.semilogy(t, np.abs(prob - 1))
    axs.set_xlabel('t')
    axs.set_ylabel(r'$\log(|1 - |\mathbf{u}^n|^2|)$')
    fig.suptitle('Deviation from normalized state')
    plt.show()
    
def measure(type, y, t, T):
    obj = pa.cx_cube()
    obj.load(f"../output/{type}_evolution.bin")
    obj = np.array(obj)
    
    M = len(obj[0, 0, :])
    N = len(obj[:, 0, 0]) 
    y_ind = int(y*M)
    t_ind = int(t/T)*N - 1
    x = np.linspace(0, 1, M)
    prob = np.absolute(obj[t_ind, y_ind, :])**2
    
    fig, axs = plt.subplots(1, 1, figsize = (8, 8))
    axs.plot(x, prob/np.sum(prob))
    axs.set_xlabel('x')
    axs.set_ylabel(f'$p(x|y = {1 - y:.2f}; t = {t})$')
    fig.suptitle(f'Detection probability \n after passing {type}-slit barrier')
    plt.show()
    
    

