import matplotlib.pyplot as plt 
from matplotlib.patches import FancyArrowPatch
import numpy as np 
import numdifftools as nd

feature_x = np.arange(-50, 50, 2)
feature_y = np.arange(-50, 50, 2)

x, y = np.meshgrid(feature_x, feature_y)
z = 0.5*(y-x)**2 + 0.5*(1-x)**2
u = 2*x - y - 1
v = y - x

# Normalize all gradients to focus on the direction not the magnitude
norm = np.linalg.norm(np.array((u, v)), axis=0)
u = u / norm
v = v / norm

fig, ax = plt.subplots(1, 1)
ax.set_aspect(1)
ax.plot(feature_x, feature_y, c='k')
ax.quiver(x, y, u, v, units='xy', scale=0.5, color='gray')
ax.contour(x, y, z, 10, cmap='jet', lw=2)

arrow = FancyArrowPatch((35, 35), (35+34*0.2, 35+0), arrowstyle='simple',
                        color='r', mutation_scale=10)  
ax.add_patch(arrow)  # NOTE: this gradient is scaled to make it better visible

plt.show() 