# %%
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import longitudinal as lde
#%matplotlib notebook

# %%
delta = 0.025
delta = 0.1
x = np.arange(-2*np.pi, 2*np.pi, delta)
y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = Y**2/2
Z2 = 2*np.sin(X/2)**2
Z = Z1 + Z2

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, 12, colors='w')
#ax.grid()
#ax.clabel(CS, inline=1, fontsize=10)
#ax.set_title('Simplest default with labels')

plt.show()

# %%
fig, ax = plt.subplots()

E0_e = 0.5109989461e6/1e3
E0_p = 938e6/1e3
Et = np.linspace(0,1e12, 100000000)/1e3
gamma_e = Et/E0_e
gamma_p = Et/E0_p

y1 = np.sqrt(1-1/gamma_e**2)
y2 = np.sqrt(1-1/gamma_p**2)

ax.semilogx(Et, y1)
ax.semilogx(Et, y2)
ax.grid()

fig.savefig("test.pdf")

plt.show()

# %%
alpha_p_LER = 3.20e-4
alpha_p_HER = 4.55e-4
gamma_t_LER = 1/np.sqrt(alpha_p_LER)
gamma_t_HER = 1/np.sqrt(alpha_p_HER)
E0 = 0.5109989461e6
Et_LER = gamma_t_LER * E0
Et_HER = gamma_t_HER * E0

print(gamma_t_LER, gamma_t_HER)
print(Et_LER*1e-6, Et_HER*1e-6)

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
