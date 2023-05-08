from sympy.vector import CoordSys3D
from sympy.vector import divergence
from sympy.vector import gradient
import numpy as np
import matplotlib.pyplot as plt

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
print('divergent at (0,2): ' + str(div.evalf(subs={R.x:0, R.y:2})))

grad = gradient(f)
print(grad)
print('gradient at (0,2): ' + str(grad.evalf(subs={R.x:0, R.y:2})))