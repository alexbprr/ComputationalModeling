#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:29:59 2017

@author: alexandre
"""
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

n = 50 # size of space: n x n
p = 0.3 # probability of initially panicky individuals

pdf = PdfPages("output-panic-ca.pdf")
fig = plt.figure()

def initialize():
    global config, nextconfig
    config = np.zeros([n, n])
    for x in range(0,n):
        for y in range(0,n):
            config[x, y] = 1 if np.random < p else 0
    nextconfig = np.zeros([n, n])

def observe(t):
    global config, nextconfig
    #plt.cla()
    img = plt.imshow(config, vmin = 0, vmax = 1, cmap = cm.binary)
    pdf.savefig(fig) #'panic ' + str(t) + '.pdf'
    #plt.show()

def update():
    global config, nextconfig
    for x in range(0,n):
        for y in range(0,n):
            count = 0
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    count += config[(x + dx) % n, (y + dy) % n]
            nextconfig[x, y] = 1 if count >= 4 else 0
    config, nextconfig = nextconfig, config

initialize()
observe(0)
for t in range(1,50):
    update()
    observe(t)

pdf.close()

#import pycxsimulator
#pycxsimulator.GUI().start(func=[initialize, observe, update])
