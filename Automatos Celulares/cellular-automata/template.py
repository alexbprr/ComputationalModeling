from time import sleep
from random import random
from sys import argv

class MeuCA(CA):
    def rule(self, x, y):
    ...
c = MeuCA(30,values=1,random_values=True)

plot(c, N=50, out="meuca.pdf")
