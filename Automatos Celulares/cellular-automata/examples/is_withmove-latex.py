from random import randint

VOID = 0
B = 1
NE = 2
TO = 3

class IS(CA):
    def rule(self, x, y):
        s = self[x][y]
        n = neighbors8(self, x, y, old=True)
        news = s

        p = 0.4 # probability of replication or IS migration
        pto = 0.4 # probability of toxin production or phagocytosis
        mto = 0.5 # toxin decay
        mne = 0.3 # natural death of neutrophils
        result = self.events([p, pto, mto, mne])

        q = neighbors8_states(self,x,y,old=True,n_states=4)
        b = ne = to =  0
        b = q[1]
        ne = q[2]
        to = q[3]

        if s == self.getOld(x, y):

            if s == 0: # void
                if b >= 1:
                    if result[0]:
                        news = 1 # replication
                    else:
                        news = 2 # migration
                elif ne >= 2:
                    if result[1]:
                        news = 3  # toxin production
            elif s == 2: # neutrophil
                if to >= 1:
                    news = 0 # neutrophil death
                elif result[3]:
                    news = 0
                else: # chemotaxis
                    bv = [i for i in neighbors8(self, x, y, pos=True)
                            if self[i] == B]
                    if len(bv) >= 1:
                        pos = bv[randint(0, len(bv)-1)]
                        self.move(x, y, *pos) # phagocytosis
                        news = 2 # attract more neutrophils
            elif s == 3:
                if result[2]:
                    news = 0
        return news

c = IS(50, values=[0]*50 + [1]*2 + [2])
plot(c, N=80, colors=['white','red', 'green', 'orange'], out='is_withmove.pdf',
    vmax=3, graphic=True, names=['Vazio','Bacteria', 'Neutrophil', 'Toxin'])
