# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

x = np.linspace(0, 20, 100)
plt.plot(x, np.sin(x))
plt.show()

# %% [markdown]
# ## シンクロトロン振動

# $$
#  \omega_s = \sqrt{-\frac{\omega_{RF}^2\eta e V \cos(\phi_s)}{2\pi h \beta^2 E}}
# $$