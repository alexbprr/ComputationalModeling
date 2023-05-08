import numdifftools as nd
  
g = lambda x:(x**4)+x + 1
grad1 = nd.Gradient(g)([1])
print("Gradient of x ^ 4 + x+1 at x = 1 is ", grad1)
  
def rosen(x): 
    return (1-x[0])**2 +(x[1]-x[0]**2)**2
  
grad2 = nd.Gradient(rosen)([1, 2])
print("Gradient of (1-x ^ 2)+(y-x ^ 2)^2 at (1, 2) is ", grad2)

import matplotlib.pyplot as plt 
from matplotlib.patches import FancyArrowPatch
import numpy as np 

feature_x = np.arange(-50, 50, 1) 
feature_y = np.arange(-50, 50, 1) 

# Creating 2-D grid of features 
[X, Y] = np.meshgrid(feature_x, feature_y) 


fig = plt.figure(figsize=(9,7))# plt.subplots(1, 1) 
ax = fig.gca()

def z(arg):
    X = arg[0]
    Y = arg[1]
    #return 0.5*np.array((Y-X)*(Y-X) + 0.5*(1-X)*(1-X))
    return X*Y

Z = z([X,Y])
# plots contour lines 
ax.contour(X, Y, Z, 10, cmap = 'jet') 
ax.grid(True)
ax.axis('scaled')
#ax.clabel(cp, inline=1, fontsize=10)  
ax.set_title('Contour Plot') 
ax.set_xlabel('feature_x') 
ax.set_ylabel('feature_y') 

gz = nd.Gradient(z)

x = 10
y = 10
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='k', mutation_scale=10)
ax.add_patch(arrow)

x = -40
y = 40
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='b', mutation_scale=10)
ax.add_patch(arrow)

x = -10
y = -10 
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='k', mutation_scale=10)
ax.add_patch(arrow)

x = 10
y = -10
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='k', mutation_scale=10)
ax.add_patch(arrow)

x = -10
y = 10
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='k', mutation_scale=10)
ax.add_patch(arrow)

x = 40
y = 40
dz1_dx = gz([x,y])
arrow = FancyArrowPatch((x, y), (x+dz1_dx[0], y+dz1_dx[1]), arrowstyle='simple', color='r', mutation_scale=10)
ax.add_patch(arrow)

plt.show() 