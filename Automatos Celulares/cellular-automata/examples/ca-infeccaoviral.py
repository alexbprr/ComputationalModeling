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

VOID = 0
VIRUS = 1
APC = 2
APM = 3
I = 4
THN = 5
THE = 6
TKN = 7
TKE = 8
C = 9

class CA_COVID(CA):

    def moveTo(self, x, y, target_state, returned_state): 
        n = [i for i in neighbors8(self, x, y, pos=True)
                if self.getOld(*i) == target_state] 
        if len(n) > 0: 
            p = n[randint(0, len(n)-1)]
            self.move(x, y, *p) 
            return returned_state 
        return self[x][y]

    def moveTov2(self, x, y, target_state, default_state, returned_state): 
        n = [i for i in neighbors8(self, x, y, pos=True)
                if self.getOld(*i) == target_state] 
        if len(n) > 0: 
            p = n[randint(0, len(n)-1)]
            self.move(x, y, *p) 
            return returned_state 
        else: 
            n = [i for i in neighbors8(self, x, y, pos=True)
                if self.getOld(*i) == default_state]
            if len(n) > 0:
                p = n[randint(0, len(n)-1)]
                self.move(x, y, *p) 
                return returned_state 
        return self[x][y]
    
    def moveToVirus(self, x, y):
        return moveTov2(self,x,y,VIRUS,VOID,VOID)

    def moveToInfectedCells(self, x, y): 
        return moveTov2(self,x,y,I,VOID,VOID)
    
    def moveToVoid(self, x, y): 
        return moveTo(self,x,y,VOID,VOID)

    def rule(self, x, y):
        s = self[x][y]
        n = neighbors8(self, x, y, old=True)        

        # probabilidades de "criar" novas células 
        p_apc = 0.3
        p_thn = 0
        p_tkn = 0.3

        # outras probabilidades 
        p_apm_death = 0.1
        p_thn_death = 0.1 
        p_tkn_death = 0.1 
        p_c_decay = 0.6
        p_the_death = 0.4 
        p_tke_death = 0.4             

        p_infection = 1
        p_replicacao = 0.4
        p_releasevirus = 0.5 
        p_infec_elim = 0.5
        p_rep_the = 0.15
        p_rep_tke = 0.3
        p_kill_virus = 1

        # quantidade de virus, apm, c, i, the e tke na vizinhança 
        virus = len([i for i in neighbors8(self, x, y) if i == VIRUS])
        apm = len([i for i in neighbors8(self, x, y) if i == APM])
        c = len([i for i in neighbors8(self, x, y) if i == C])
        i = len([i for i in neighbors8(self, x, y) if i == I])
        the = len([i for i in neighbors8(self, x, y) if i == THE])
        tke = len([i for i in neighbors8(self, x, y) if i == TKE])

        # regras do CA
        if s == self.getOld(x, y):
            
            if s == VOID: 
                if c >= 1:
                    if random() < p_apc: 
                        return APC 
                    #elif random() < p_apc + p_thn:  
                    #    return THE 
                    elif random() <= p_apc + p_thn + p_tkn: 
                        return TKE                  

                #if i >= 1 and random() < 0.5: 
                #    return TKE                      
                
                if the >= 1 and random() < p_rep_the:
                    return THE                 
                if tke >= 1 and random() < p_rep_tke:
                    return TKE      
                if ((apm >= 1 or the >= 1 or tke >= 1)) and random() < 0.2:
                    return C                 

                if virus >= 1: 
                    if random() < p_replicacao: 
                        return VIRUS
                    else: 
                        if random() < p_apc: 
                            return APC 
                        #elif random() < p_apc + p_thn:  
                        #    return THE 
                        elif random() <= p_apc + p_tkn: 
                            return TKE             
                
                if i >= 1 and random() < p_releasevirus: 
                    return VIRUS
                         
            if s == VIRUS:
                if (apm >= 1) and random() < p_kill_virus: 
                    return VOID 
                #return self.moveToVoid(x,y)
                return s 
                
            if s == APC:
                if virus >= 1 and random() < 0.75: 
                    return APM 

                if virus >= 1 and random() < p_infection: 
                    return I 
                    
                if random() < p_apm_death:
                    return VOID
                
                return self.moveToVirus(x,y) 
                
            if s == APM: 
                if random() < p_apm_death:
                    return VOID 
                if virus >= 1 and random() < p_infection: 
                    return I 
                
                return self.moveToVirus(x,y)

            if s == I:
                if tke >= 1 and random() < p_infec_elim:
                    return VOID 

            if s == C:
                if random() < p_c_decay:
                    return VOID 
                return self.moveToVoid(x,y)

            if s == THN:
                if apm >= 1: 
                    return THE 
                if random() < p_thn_death:
                    return VOID 
                                
            if s == THE: 
                if random() < p_the_death:
                    return VOID 
                
                return self.moveToVirus(x,y)
                
            if s == TKN:
                if apm >= 1: 
                    return TKE
                if random() < p_tkn_death:
                    return VOID                 
            
            if s == TKE: 
                if random() < p_tke_death:
                    return VOID 
                if virus >= 1 and random() < p_infection:
                    return I 
                
                return self.moveToInfectedCells(x,y)

        return s

c = CA_COVID(80, values=[0]*100 + [1]*5 + [2]*25 + [5]*5 + [7]*5)

#Void,Virus,Apc,Apm,I,Thn,The,Tkn,Tke,C

plot(c, N=150, out='ca-infecacaoviral-res.pdf',
    vmax=9, graphic=True,
    colors=['white','red','white', 'darkblue', 'yellow', 'white', 'darkgreen', 'white', 'purple', 'orange'], 
    names=['Empty','Virus','Apc','Apm','Infected cells', 'Thn','The','Tkn','Tke','Proinflammatory cytokines'])
