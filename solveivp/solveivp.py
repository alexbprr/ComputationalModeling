# -*- coding: utf-8 -*-
from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

def odeSystem(t, P, r, k):
    N = P[0] 
    dN_dt = r*N*(1 - N/k) 
    return [dN_dt]  

def solveOde(t, y):
    global model_args
    return odeSystem(t, y, *model_args)

def setInitialCondition():
    N = 10
    P = [N]
    return P

def setParameters():
    r = 0.05
    k = 1000
    model_args = (r, k)
    return model_args

def setTime():
    dt = 0.01    
    tfinal = 20
    nsteps = int(tfinal/dt)
    print('nsteps ' + str(nsteps))
    t=np.linspace(0,tfinal,nsteps)
    return (tfinal,t)

def saveResults(fname, populations, values):
    df = pd.DataFrame(values.y.transpose(), columns = populations)
    df.insert(0, 'Time', values.t)
    df.to_csv(fname, float_format='%.5f', sep=',')

def readFile(filename):
    file = open(filename, 'r')
    if (file):        
        return pd.read_csv(filename)
    else: 
        print('Error: it was not possible to open the file!')
        exit

def createFigure(title):
    fig = plt.figure(figsize=(9,7))    
    plt.tick_params()
    plt.title(title)
    plt.xlabel('time')
    plt.ylabel('[Y]')
    ax = fig.gca()
    return [fig,ax]
    
def initializeGraphParameters():
    global scalarMap, color_index, colorVal
    cm = plt.get_cmap('jet')
    cNorm  = colors.Normalize(vmin=0, vmax=30)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    color_index = 0
    colorVal = scalarMap.to_rgba(color_index)

def plotFig(name, t, P):    
    global color_index, colorVal
    initializeGraphParameters()
    dir_ = ''
    fig,ax = createFigure(name)
    color_index += 50/4
    ax.plot(t, P, label=name, color=colorVal)
    colorVal = scalarMap.to_rgba(color_index)
    ax.legend(loc='best')
    plt.show()
    fig.savefig(dir_ + name + '.svg', format='svg', bbox_inches='tight')
    fig.clear()
    ax.clear()

if __name__ == "__main__":
    global model_args
    (tfinal,t) = setTime()
    P = setInitialCondition()
    model_args = setParameters()     
    results = solve_ivp(solveOde,(0, tfinal), P, t_eval=t, method='Radau')
    saveResults('ode_results_solveivp.csv', ['N'], results)
    odeValues = readFile('ode_results_solveivp.csv')
    plotFig('N',odeValues['Time'],odeValues['N'])