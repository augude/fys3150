#%%

import numpy as np; import matplotlib.pyplot as plt 
import seaborn as sns
import tikzplotlib

sns.set_theme(font_scale = 2)

data = np.loadtxt('ElectricField.txt')

N = 100
Ex = np.zeros(shape = (N, N))
Ey = np.zeros(shape = (N, N))
Ez = np.zeros(shape = (N, N))

X = np.zeros(shape = (N, N))
Y = np.zeros(shape = (N, N))
Z = np.zeros(shape = (N, N))

V = np.zeros(shape = (N, N))

for i in range(N**2):
    ind1 = int(data[i, 0])
    ind2 = int(data[i, 1])
    X[ind1, ind2] = data[i, 2]
    Y[ind1, ind2] = data[i, 3]
    Z[ind1, ind2] = data[i, 4] 

    Ex[ind1, ind2] = data[i, 5]
    Ey[ind1, ind2] = data[i, 6]
    Ez[ind1, ind2] = data[i, 7]  
    V[ind1, ind2] = data[i, 8]       
    
fig, axs = plt.subplots(ncols = 1, nrows = 1, figsize = (12, 10))
plt.contourf(Y, Z, V, cmap = 'cool')
cbar = plt.colorbar(cmap = 'cool')
cbar.set_label(r'Electric potential [$\frac{u(\mu m)^2}{(\mu s)^2 e}$]')
axs.quiver(Y[::5, ::5], Z[::5, ::5], Ey[::5, ::5], Ez[::5, ::5])
axs.set_xlabel(r'y [$\mu m$]')
axs.set_ylabel(r'z [$\mu m$]')
fig.tight_layout()
#plt.savefig('PenningTrapSetup.pdf')
plt.show()

# %%
files = ['testOneParticleFE.txt', 'testOneParticleRK4.txt']

for filename in files:
    test = np.loadtxt(filename)
    t = test[:, 0]
    pos = test[:, 1:4]
    vel = test[:, 4:]

    q = 1; V0 = 9.65e8; B0 = 9.65e1; m = 40.078; d = 1e4
    wz = np.sqrt(2*q*V0/(m*d**2)) 
    w0 = q*B0/m

    v0 = vel[0, 2]; x0 = pos[0, 0]; z0 = pos[0, 2]

    wp = (w0 + np.sqrt(w0**2 - 2*wz**2))/2
    wm = (w0 - np.sqrt(w0**2 - 2*wz**2))/2

    Ap = (v0 + wm*x0)/(wm - wp)
    Am = -(v0 + wp*x0)/(wm - wp)

    f = Ap*np.exp(-1j*wp*t) + Am*np.exp(-1j*wm*t)
    xAna = np.real(f)
    yAna = np.imag(f)
    zAna = z0*np.cos(wz*t)

    fig, axs = plt.subplots(3, 1, figsize = (10, 10))
    axs[0].plot(t, pos[:, 0])
    axs[0].plot(t, xAna, linestyle = '--')
    axs[0].set_ylabel(r'x [$\mu m$]')
    axs[0].set_xlabel(r't [$\mu s$]')

    axs[1].plot(t, pos[:, 1])
    axs[1].plot(t, yAna, linestyle = '--')
    axs[1].set_ylabel(r'y [$\mu m$]')
    axs[1].set_xlabel(r't [$\mu s$]')

    axs[2].plot(t, pos[:, 2])
    axs[2].plot(t, zAna, linestyle = '--')
    axs[2].set_ylabel(r'z [$\mu m$]')
    axs[2].set_xlabel(r't [$\mu s$]')
    plt.legend(['Numerical', 'Analytical'])
    savename = filename.replace('.txt', '')
    fig.tight_layout()
    #plt.savefig(f'{savename}.pdf')
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}.tex",
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
files = ['testDoubleSetupWithInternal.txt', 'testDoubleSetupWithoutInternal.txt']

for filename in files:
    savename = filename.replace('.txt', '')
    savename = savename.replace('testDoubleSetup', '')
    print(savename)
    data = np.loadtxt(filename)
    t = data[:, 0]
    pos1 = data[:, 1:4]
    vel1 = data[:, 4:7]
    pos2 = data[:, 7:10]
    vel2 = data[:, 10:13]
    
    r = np.sqrt( (pos1[:, 0] - pos2[:, 0])**2 + (pos1[:, 1] - pos2[:, 1])**2 + (pos1[:, 2] - pos2[:, 2])**2) 
    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    axs.plot(t, r)
    axs.set_xlabel(r't [$\mu s$]')
    axs.set_ylabel(r'r [$\mu m$]')
    fig.tight_layout()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}Distance.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()

    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    axs.plot(pos1[:, 0], pos1[:, 1])
    axs.plot(pos2[:, 0], pos2[:, 1], linestyle = '--')
    axs.set_xlabel(r'x [$\mu m$]')
    axs.set_ylabel(r'y [$\mu m$]')
    plt.legend(['Particle 1', 'Particle 2'])
    fig.tight_layout()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}XY.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()
    
    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    axs.plot(pos1[:, 0], vel1[:, 0])
    axs.plot(pos2[:, 0], vel2[:, 0], linestyle = '--')
    axs.set_xlabel(r'x [$\mu m$]')
    axs.set_ylabel(r'$v_x$ [$m/s$]')
    plt.legend(['Particle 1', 'Particle 2'])
    fig.tight_layout()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}PhaseX.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()
    
    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    axs.plot(pos1[:, 1], vel1[:, 1])
    axs.plot(pos2[:, 1], vel2[:, 1], linestyle = '--')
    axs.set_xlabel(r'y [$\mu m$]')
    axs.set_ylabel(r'$v_y$ [$m/s$]')
    plt.legend(['Particle 1', 'Particle 2'])
    fig.tight_layout()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}PhaseY.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()
    
    fig, axs = plt.subplots(1, 1, figsize = (10, 10))
    axs.plot(pos1[:, 2], vel1[:, 2])
    axs.plot(pos2[:, 2], vel2[:, 2], linestyle = '--')
    axs.set_xlabel(r'z [$\mu m$]')
    axs.set_ylabel(r'$v_z$ [$m/s$]')
    plt.legend(['Particle 1', 'Particle 2'])
    fig.tight_layout()
    tikzplotlib.clean_figure()
    tikzplotlib.save(
    f"{savename}PhaseZ.tex",
    extra_axis_parameters=[
        "title style={align=center}",
        "xmajorticks=true",
        "ymajorticks=true",
        "mark options={mark size=2.5pt, line width=1.5pt}",
        ],
        strict=True,
    )
    plt.show()
    
    fig = plt.figure(figsize = (10, 10))
    axs = fig.add_subplot(1, 1, 1, projection = '3d')
    axs.plot3D(pos1[:, 0], pos1[:, 1], pos1[:, 2])
    axs.plot3D(pos2[:, 0], pos2[:, 1], pos2[:, 2], linestyle = '--')
    axs.set_xlabel(r'x [$\mu m$]')
    axs.set_ylabel(r'y [$\mu m$]')
    axs.set_zlabel(r'z [$\mu m$]')
    axs.set_xticks([-8000, -4000, 4000, 8000])
    axs.set_yticks([-8000, -4000, 4000, 8000])
    axs.set_zticks([-6000, -3000, 3000, 6000])
    plt.legend(['Particle 1', 'Particle 2'])
    fig.tight_layout()
    plt.savefig(f"{savename}3D.pdf")
    plt.show()
    

# %%
