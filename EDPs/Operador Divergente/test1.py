from sympy.vector import CoordSys3D
from sympy.vector import divergence
from sympy.vector import gradient

R = CoordSys3D('R') #Sistema de Coordenadas 3D definido em um frame 'R' 
vecField1 = 3*R.i + 4*R.j + 5*R.k
print(vecField1)
div1 = divergence(vecField1)
print(div1)
vecField2 = R.x**2*R.i + 2* R.x*R.z *R.k  #(x²,0,2.x.z)
print(vecField2) #2x + 0 + 2x = 4x  
div2 = divergence(vecField2)
print(div2)

f = R.x**2 + R.y**2*R.x 
vecField = f*R.i + f*R.j + f*R.k
div = divergence(vecField)
print(div)

grad = gradient(f)
print(grad)

print('gradient at (0,2): ' + str(grad.evalf(subs={R.x:0, R.y:2})))

import numdifftools as nd
import numpy as np
import matplotlib.pyplot as plt
# Divergence of vector-valued function f evaluated at x
#def div(f,x):
#    jac = nd.Jacobian(f)(x)
#    return jac[0,0] + jac[1,1] + jac[2,2]

dim= 20
xarray= np.arange(-dim,dim)
yarray= np.arange(-dim,dim)
x,y = np.meshgrid(xarray,yarray)
#vx = (x-L)/((x-L)**2+y**2) - (x+L)/((x+L)**2 +y**2)
#vy = y/((x-L)**2+y**2) - y/((x+L)**2 +y**2)

threshold = 0.33
Mx = np.abs(vx) > threshold
My = np.abs(vy) > threshold
vx = np.ma.masked_array(vx, mask=Mx)
vy = np.ma.masked_array(vy, mask=My)

plt.figure()
plt.quiver(x,y, vx, vy, pivot='mid')
plt.xlabel("$x$-axis")
plt.ylabel("$y$-axis")
plt.axis('equal')
plt.show()