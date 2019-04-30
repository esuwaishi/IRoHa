# %%
import numpy as np
import matplotlib.pyplot as plt
import longitudinal as lde
#%matplotlib notebook

# %%
turns=500
sin_phi_s=0.5
#phi_s=0
energy_list=[250e9,20e9]
mass=0.938e9
harm=360
voltage=5e6
mcf=0.0018
update_eta=True

#print(height)

# Vicinity of the pi-phi_s
fig,ax=plt.subplots(1,2,figsize=(10,4))

npar=5
result_g1=[]
for i in range(len(energy_list)):
    beam_energy=energy_list[i]
    gamma=beam_energy/mass
    eta=mcf-1/gamma/gamma
    beta=np.sqrt(1-1/gamma/gamma)
    phi_s=np.arcsin(sin_phi_s)
    yf=np.sqrt(np.cos(phi_s)-(np.pi-2*phi_s)*np.sin(phi_s)/2)
    height=2*np.sqrt(voltage/2/np.pi/beta/beta/beam_energy/harm/np.abs(eta))*yf
    
    initial_phi=np.ones(npar)*(np.pi-phi_s)
    initial_delta=np.linspace(height/npar, height*0.99, npar)
    temp=lde.longitudinal_evolve_delta(turns, 
                                       initial_phi,
                                       initial_delta,
                                       sin_phi_s=sin_phi_s, alphac=mcf, E0_ini=beam_energy,
                                       mass=mass, e_volt=voltage, harm=harm,
                                       update_eta=update_eta
                                      )
    result_g1.append(temp)

# Vicinity of the phi_s

result_g2=[]
for i in range(len(energy_list)):
    beam_energy=energy_list[i]
    gamma=beam_energy/mass
    eta=mcf-1/gamma/gamma
    beta=np.sqrt(1-1/gamma/gamma)
    phi_s=np.arcsin(sin_phi_s)
    yf=np.sqrt(np.cos(phi_s)-(np.pi-2*phi_s)*np.sin(phi_s)/2)
    height=2*np.sqrt(voltage/2/np.pi/beta/beta/beam_energy/harm/np.abs(eta))*yf
    ax[i].set_ylim([-2*height,2*height])
    initial_phi=np.ones(npar)*(phi_s)
    initial_delta=np.linspace(0, height*0.99, npar)
    temp=lde.longitudinal_evolve_delta(turns, 
                                       initial_phi,
                                       initial_delta,
                                       sin_phi_s=sin_phi_s, alphac=mcf, E0_ini=beam_energy,
                                       mass=mass, e_volt=voltage, harm=harm, 
                                       update_eta= update_eta
                                      )
    result_g2.append(temp)


for i in range(2):
    ax[i].set_aspect('auto')
    ax[i].set_title("{} GeV".format(energy_list[i]/1e9))
    ax[i].set_xlabel('Phase [rad]')
    
    ax[i].set_xlim([-np.pi,2*np.pi])
    
    #ax[i].set_prop_cycle(plt.cycler('color', plt.cm.Oranges(np.linspace(0.3, 0.8, npar))))
    ax[i].set_prop_cycle(plt.cycler('color', plt.cm.Oranges(np.linspace(0.3, 0.8, npar))))
    ax[i].plot(result_g1[i][0], result_g1[i][1], linestyle=None)
    ax[i].set_prop_cycle(plt.cycler('color', plt.cm.Blues(np.linspace(0.3, 0.8, npar))))
    ax[i].plot(result_g2[i][0], result_g2[i][1], linestyle=None)
    ax[i].set_xticks([-np.pi, -0.5*np.pi, 0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
    ax[i].set_xticklabels([r"$-\pi$", r"$-\frac{1}{2}\pi$", "$0$", r"$\frac{1}{2}\pi$",
                     r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
ax[0].set_ylabel(r'$\delta=dp/p_0$')
plt.tight_layout()

# %%
fig,ax=plt.subplots()
ax.set_xlabel('Phase [rad]')
ax.set_ylabel('Potential [Arb. unit]')
phi=np.linspace(-np.pi,4*np.pi,100)
phi_s=np.pi/6
ax.plot(phi, (np.cos(phi)-np.cos(phi_s)+np.sin(phi_s)*(phi-phi_s)), label=r'$\eta<0$; $\phi_{stable}=5\pi/6$')
phi_s=0
ax.plot(phi, (np.cos(phi)-np.cos(phi_s)+np.sin(phi_s)*(phi-phi_s)), label=r'$\eta<0$; $\phi_{stable}=\pi$')
ax.axhline(y=0, ls='--', c='g')
ax.set_xticks([-np.pi, -0.5*np.pi, 0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
ax.set_xticklabels([r"$-\pi$", r"$-\frac{1}{2}\pi$", "$0$", r"$\frac{1}{2}\pi$",
                     r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
ax.legend(loc='best')

# %% [markdown]
# ## シンクロトロン振動

# $$
#  \omega_s = \sqrt{-\frac{\omega_{RF}^2\eta e V \cos(\phi_s)}{2\pi h \beta^2 E}}
# $$

# %%
