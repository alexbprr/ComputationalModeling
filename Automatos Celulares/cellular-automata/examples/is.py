try:
    from catemp import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from catemp import *

from time import sleep
from random import random
from sys import argv

# 41 -> red    dead
# 42 -> green  bacteria
# 43 -> yellow neutrophil
# 0 dead
# 1 bacteria
# 2 neutrophil
# 3 toxin

#Atualizacao do estado de um vizinho especifico

class IS(CA):
    def rule(self, x, y):
        s = self[x][y]
        n = neighbors8(self, x, y, old=True)
        news = s

        p = 0.5 #probability of replication or IS migration
        pto = 0.3 #probability of phagocytosis or toxin production
        mto = 0.4 #toxin decay
        b = ne = to =  0

        for i in n:
            if i == 1: b += 1
            elif i == 2: ne += 1
            elif i == 3: to += 1

        q = neighbors8_states(self,x,y,False,3)
        b = q[1]
        ne = q[2]
        to = q[3]

        if s == 0: #dead
            if b >= 2:
                if random() < p:
                    news = 1 # replication
                else:
                    news = 2 # migration
        elif s == 1: #bacteria
            if ne >= 1:
                if random() < pto:
                    news = 0 # phagocytosis
                else:
                    news = 3 # toxin production
        elif s == 2: #neutrophil
            if to >= 1:
                news = 0 #neutrophil death

        return news

    def prettyPrint(self, x, y):
        return b"\033[%dm  \033[m" % (41+self[x][y])

c = IS(30, values=(1, 1, 0, 0, 0, 0, 0, 0, 2))
t = 0

plot(c, N=50, colors=['white','red', 'green', 'orange'], out='sir.pdf', graphic=True, vmax=3)
