from random import randint, random

RH = 0.3 # Reprodução da presa
RP = 0.1 # Reprodução do predador
MP = 0.1 # Morte do predador

VOID = 0
PREY = 1
PREDATOR = 2

class PP(CA):

    def rule(self, x, y):
        s = self[x, y]

        # Evita que posições já atualizadas se re-atualizem (Caso específico para movimentação)
        if s == self.getOld(x, y):

            if  s == VOID:

                # Se há 2 ou mais presas vizinhas, há uma probabilidade de
                # nascer uma terceira nesta posição
                n = [i for i in neighbors8(self, x, y) if i == PREY]
                if len(n) >= 2:
                    if random() < RH: return PREY
                    else: return VOID
                else: return VOID

            elif s == PREY:

                # Se há posições vazias em volta, mova-se para alguma
                n = [i for i in neighbors8(self, x, y, pos=True)
                     if i == VOID]

                if len(n) > 0:
                    p = n[randint(0, len(n)-1)]
                    self.move(x, y, *p)

                    return VOID

                else: return s

            elif s == PREDATOR:

                # Se há alguma presa em volta, prede-a
                n = [i for i in neighbors8(self, x, y, pos=True)
                     if self[i] == PREY]

                if len(n) > 0:
                    p = n[randint(0, len(n)-1)]
                    self.move(x, y, *p)

                    # Se há algum outro predador em volta, há uma chance de
                    # reprodução
                    n = [i for i in neighbors8(self, *p) if i == PREDATOR]
                    if len(n) > 1 and random() < RP: return PREDATOR
                    else: return VOID

                # Mortalidade
                elif random() < MP: return VOID

                # Movimentação aleatória, em último caso
                else:

                    # Se há alguma posição vazia, mova-se para ela
                    n = [i for i in neighbors8(self, x, y, pos=True)
                         if self[i] == VOID]

                    if len(n) > 0:
                        p = n[randint(0, len(n)-1)]
                        self.move(x, y, *p)

                        return VOID

                    else: return s
            else: return s
        else: return s

c = PP(30, values=[0]*50 + [1]*2 + [2])
plot(c, N=70, graphic=True)
