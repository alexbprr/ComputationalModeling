import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pandas as pd

dt = 1e-03
tfinal = 50
T = tfinal*(1.0/dt)
t = numpy.arange(0,tfinal+dt,dt)
numIteracoesGravadas = 10
interval = T/numIteracoesGravadas
print(t)

dx = 0.1
xdim = 10
xsize = int(xdim/dx)
xarray = np.arange(0,xdim,dx)
x = np.arange(0,xsize,1)

print(xsize)
print(xarray)
print(x)

maxDensity = 1.0 
keqDensity = 0.5

r_b = 0.
Db = 2
Dn = 0.05
Xn = 0.

c = 0.

"""
wbf = 1.4
wff = 1.
wnf = 1.3
wbb = 1.
wfb = 1.4
wnb = 1.
wbn = 1.4
wfn = 4.
wnn = 1.
"""

def verifyDensity(value):
    if (value < 0.): 
        return 0.
    else: 
        return value

def verifyMaxDensity(value):
    if (value > maxDensity):
        value = maxDensity
    return value

def h(v):
    return 0.1*v/(1 + v)

def diffusionX(u, x, t):
    if (x == 0):    
        resX = (u[t][x+1] - u[t][x])/(dx*dx);    
    elif (x == (xsize -1)):    
        resX = (u[t][x-1] - u[t][x])/(dx*dx);    
    else:    
        resX = (u[t][x+1] - 2*u[t][x] + u[t][x-1])/(dx*dx); #dentro do dominio mas fora da extremidade    
    return resX;  

alpha = 1    
bfd = 0
nfd = 0
gradbfd = 0
bfe = 0
nfe = 0
gradbfe = 0
beta = 10 #Influencia de n na difusao de B 

def fronteiraDireita(b, n, x, t):
    global bfd, nfd, gradbfd
    bfd = (b[t][x+1] + b[t][x])/2
    nfd = (n[t][x+1] + n[t][x])/2
    gradbfd = (b[t][x+1] - b[t][x])

def fronteiraEsquerda(b, n, x, t):
    global bfe, nfe, gradbfe
    bfe = (b[t][x] + b[t][x-1])/2
    nfe = (n[t][x] + n[t][x-1])/2
    gradbfe = (b[t][x] - b[t][x-1])
      
def difNaoLinear(b, n, x, t):
    global bfd, nfd, gradbfd, bfe, nfe, gradbfe
    if(x == 0):  
        fronteiraDireita(b, n, x, t)
        bfe = (b[t][x])/2
        nfe = (n[t][x])/2
        gradbfe = 0                
        return alpha*(bfd/(1 + beta*nfd)*gradbfd)/(dx*dx) 
    elif(x == xsize-1):         
        fronteiraEsquerda(b, n, x, t)
        bfd = (b[t][x])/2
        nfd = (n[t][x])/2
        gradbfd = 0
        return alpha*(- bfe/(1 + beta*nfe)*gradbfe)/(dx*dx ) 
    else:
        fronteiraDireita(b, n, x, t)
        fronteiraEsquerda(b, n, x, t)
        return alpha*(bfd/(1 + beta*nfd)*gradbfd - bfe/(1 + beta*nfe)*gradbfe) /(dx*dx ) 
    
limits = [30,35,80,90]
concB = [0.7,0.2]
concN = [0.]
concF = [0.]

bstate = []
fstate = []
nstate = []
    
tAntigo = 0
tAtual = 1

b = np.zeros((2,xsize))
f = np.zeros((2,xsize))
n = np.zeros((2,xsize))

for i in range(0,xsize):    
    if (i >= limits[0] and i < limits[1]): 
        b[tAntigo][i] = concB[0]
        
    f[tAntigo][i] =  concF[0]
    n[tAntigo][i] =  concN[0]    
 
    #updateW(b,f,n,wb,wf,wn,gwb,gwf,gwn,i,tAntigo)

bstate.append(b[0])
fstate.append(f[0])
nstate.append(n[0])

vmax_color=10
cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=vmax_color)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
color_index = 0
colorVal = scalarMap.to_rgba(color_index)

fig1 = plt.figure(figsize=(9,7))
ax1 = fig1.gca()
plt.axis([0, x[xsize-1]+0.1, 0, maxDensity])
plt.xlabel('x(in mm)',fontsize=15)
plt.ylabel('concentration',fontsize=15)
ax1.plot(x,b[tAntigo],label="{:.0f}".format(t[0]), color=colorVal) 

fig2 = plt.figure(figsize=(9,7))
ax2 = fig2.gca()
plt.axis([0, x[xsize-1]+0.1, 0, maxDensity])
plt.xlabel('x(in mm)',fontsize=15)
plt.ylabel('concentration',fontsize=15)
ax2.plot(x,n[tAntigo],label="{:.0f}".format(t[0]), color=colorVal) 
    
for step in range(1,int(T+1)):  

    if(step % 2 == 0):
        tAtual = 0
        tAntigo = 1
    else: 	
        tAtual = 1
        tAntigo = 0

    for i in range(0,xsize):        
        b[tAtual][i] = b[tAntigo][i] + ( r_b*b[tAntigo][i] + Db*difNaoLinear(b,n,i,tAntigo) )*dt 
        b[tAtual][i] = verifyMaxDensity(b[tAtual][i])
        
        f[tAtual][i] = 1; 
        f[tAtual][i] = verifyMaxDensity(f[tAtual][i])

        if(i == 50):# or i == 70
            s = 10
        else:
            s = 0
            
        n[tAtual][i] = n[tAntigo][i] + ( s - c*b[tAntigo][i]*n[tAntigo][i] + Dn*diffusionX(n,i,tAntigo) )*dt
        n[tAtual][i] = verifyMaxDensity(n[tAtual][i])
        
        #updateW(b,f,n,wb,wf,wn,gwb,gwf,gwn,x,tAtual)
    bstate.append(b[tAtual])
    fstate.append(f[tAtual])
    nstate.append(n[tAtual])

    if(step % interval == 0):
        #print("salvando step " + str(step));        
        color_index += 1
        colorVal = scalarMap.to_rgba(color_index)
        ax1.plot(x,b[tAtual], label="{:.0f}".format(t[step]), color=colorVal)
        ax2.plot(x,n[tAtual], label="{:.0f}".format(t[step]), color=colorVal) 

bstate = np.array(bstate)
fstate = np.array(fstate)
nstate = np.array(nstate)

df = pd.DataFrame(bstate, columns = [xarray])
#df.insert(0, 'Time', t)
df.to_csv('b_x.csv', float_format='%.5f', sep=',')

ax1.legend(loc='best')
ax2.legend(loc='best')        
#ax3.legend(loc='best')
       
fig1.savefig('b_difusao_naolinear.svg', format='svg', bbox_inches='tight')
fig2.savefig('n_difusao_naolinear.svg', format='svg', bbox_inches='tight')
#fig3.savefig('f_difusao_naolinear.svg', format='svg', bbox_inches='tight') 