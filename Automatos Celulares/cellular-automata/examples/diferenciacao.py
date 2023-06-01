try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *

from time import sleep
from random import randint, random, choice
from sys import argv
import random as rd

S0 = 0
S1 = 1

class Differ(CA):
    def rule(self, x, y):
        s = self[x][y]
        n = neighbors8(self, x, y, old=True)                
        
        n0 = len([i for i in neighbors8(self, x, y) if i == S0])
        n1 = len([i for i in neighbors8(self, x, y) if i == S1])

        pdecay = 1 
        pgrow = 0.8
        
        # regras do CA
        if s == self.getOld(x, y):            
            if s == S0:
                if n0 == 8 and random() < pgrow:
                    return S1 
            if s == S1: 
                if n1 >= 1 and random() < pdecay:
                    return S0                         
            return s 

c = Differ(100, values=1, random_values=False)
c.add(1, points=[(10, 10)], size=(10, 10))
c.add(0, points=[(11, 11)], size=(8, 8))

plot(c, N=50, out='ca_diferenciacao_res.pdf',
    vmax=1, graphic=True,
    colors=['white','black'], 
    names=['inhibited', 'active'])
