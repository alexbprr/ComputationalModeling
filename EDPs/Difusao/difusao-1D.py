#Referencia https://www.codespeedy.com/random-walk-program-in-python/
import numpy as np
import matplotlib.pyplot as plt

def plot(passo, x, u):
    plt.plot(x,u,linewidth=2)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title("Difusao")
    plt.savefig("figs/u" + str(passo) + ".svg", format="svg", bbox_inches="tight")
    
def Diffusion1D(nsteps, xsize, dx, u, traces): #n here is the no. of steps that we require
    xvalues = np.arange(0,xsize,dx)
    dt = 0.1 
    Du = 1
    plot(0, xvalues, u)

    for x in xvalues:
        traces[x][0] = u[x]
    
    for i in range (1,nsteps): # para cada passo de tempo 
        for x in xvalues: # para cada posicao x do espaco
            coeff = Du/(dx*dx)
            diffu = 0
            if (0 < x and x < xsize-1):
                diffu = (u[x+1] - 2*u[x] + u[x-1])*coeff  
            elif (x < xsize - 1):
                diffu = (u[x+1] - u[x])*coeff #u[x-1] = u[x] para x == 0
            elif (x > 0):
                diffu = (-u[x] + u[x-1])*coeff

            u[x] += diffu*dt 
            traces[x][i] = u[x]

        if (i % 10 == 0):
            plot(i, xvalues, u)    

    plt.figure(figsize=(9,7))
    times = range (0,nsteps)
    plt.plot(times,traces[10][:], color='b', label='u[10]')
    plt.plot(times,traces[20][:], color='r', label='u[20]')
    plt.plot(times,traces[30][:], color='g', label='u[30]')
    plt.plot(times,traces[40][:], color='m', label='u[40]')
    plt.xlabel('t')
    plt.ylabel('u')
    plt.legend(loc='best')
    plt.title("Evolution of u concentration at different points")
    plt.savefig("figs/u_x.svg", format="svg", bbox_inches="tight")


nsteps = 500
xsize = 50 # tamanho do dominio (de 0 a 100)
dx = 1 # discretizacao Espacial 
npoints = int(xsize/dx)
initialConcentration = 1000 # concentracao inicial  

u = np.zeros(npoints)  # concentracao em cada ponto x do espaco   
u[int(xsize/2)] = initialConcentration

traces = np.zeros((npoints,nsteps))

plt.figure(figsize=(9,7))
Diffusion1D(nsteps,xsize,dx,u, traces) 