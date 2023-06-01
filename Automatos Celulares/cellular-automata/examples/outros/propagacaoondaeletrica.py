"""
Created on Wed Feb  8 11:29:59 2017

@author: alexandre
"""
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import numpy as np
import random
from copy import deepcopy

n = 30 # size of space: n x n
p = 0.2 # Probabilidade de uma célula estar excitada inicialmente
pe = 0.9
pr = 0.9
k = 3

random.seed()

pdf = PdfPages("output-propagacaoeletrica-ca.pdf")
fig = plt.figure()

fd = open('output-propagacaoeletrica-ca.txt','w')

def initialize():
    global config, nextconfig
    config = np.zeros([n, n])
    for x in range(0,n):
        for y in range(0,n):
            if (x >= n-20 and x <= n-10 and y >= n/2 and y <= (n/2 + 10)):
                config[x,y] = 1
            #elif (x >= n-10 and x < n and y > n/2+10 and y < n):
            #    config[x,y] = -1
            #config[x, y] = 1 if random.random() < p else 0
    nextconfig = np.zeros([n, n])

def observe(t):
    global config, nextconfig
    cmap = colors.ListedColormap(['bisque', 'red', 'chocolate'])
    #bounds=[0,1,2]
    #norm = colors.BoundaryNorm(bounds, cmap.N)
    #img = plt.imshow(config, vmin = 0, vmax = k, cmap = cmap)

    fig = plt.figure(figsize=(10,7))
    plt.axis([0, n, 0, n])
    plt.pcolor(config, vmin=0, vmax=k)
    plt.colorbar()
    plt.title('Propagacao eletrica' ,fontsize=16)
    plt.xlabel('x',fontsize=15)
    plt.ylabel('y',fontsize=15)

    pdf.savefig(fig)
    plt.close(fig)

    for line in np.matrix(config):
        fd.write(str(line)) #, fmt='%d'
    fd.write('\n')

def update():
    global config, nextconfig
    nextconfig = deepcopy(config)
    for x in range(0,n):
        for y in range(0,n):
            count = 0
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if (config[(x + dx) % n, (y + dy) % n] == 1):
                        count += 1;
            if (config[x,y] == 0 and count >= 2 and random.random() < pe): #normal -> exitado
                nextconfig[x,y] = 1
            elif (config[x,y] == 1): #excitado -> refratário
                nextconfig[x,y] = 2
            elif (config[x,y] == k): #refratário -> normal
                nextconfig[x,y] = 0
            elif (config[x,y] >= 2 and config[x,y] < k): #refratário -> refratário
                nextconfig[x,y] = config[x,y] + 1
    config = deepcopy(nextconfig)

initialize()
observe(0)
for t in range(1,100):
    update()
    observe(t)

pdf.close()
fd.close()
