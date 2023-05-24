import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys 
import os

dataPath = 'results/'
figPath = 'figs/'

x = np.loadtxt(dataPath+ 'x.csv')
y = np.loadtxt(dataPath+ 'y.csv')
t = np.loadtxt(dataPath+ 't.csv')

def saveFigure(data, name, title, minDensity, maxDensity):
    cm = plt.get_cmap('viridis')
    fig = plt.figure(figsize=(10,7))
    plt.title(title, fontsize=14)
    plt.xlabel('x',fontsize=14)
    plt.ylabel('y',fontsize=14)
    pcm = plt.pcolormesh(x, y, data, cmap=cm, vmin=minDensity, vmax=maxDensity)
    plt.colorbar(pcm)
    fig.savefig(name + '.svg', format='svg', bbox_inches='tight')
    plt.close()

def getMaxMinValueFrom2DMatrix(u):
    maxv = np.max(np.max(u))
    minv = np.min(np.min(u))
    return [maxv,minv]

inc = int(sys.argv[1])
values = [i*inc for i in range(0,11)]
for v in values: 
    bname = 'bacteria' + str(v) 
    nname = 'neutrophil' + str(v)
    b = np.loadtxt(dataPath + bname + '.csv', delimiter=',')
    n = np.loadtxt(dataPath + nname + '.csv', delimiter=',')
    [maxv,minv] = getMaxMinValueFrom2DMatrix(b)
    saveFigure(b,figPath+bname, 'Bacteria', maxv, minv)
    [maxv,minv] = getMaxMinValueFrom2DMatrix(n)
    saveFigure(n,figPath+nname, 'Neutrophil', maxv, minv)