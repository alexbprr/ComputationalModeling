#Referencia https://www.codespeedy.com/random-walk-program-in-python/
import numpy as np
import matplotlib.pyplot as plt
def Randomwalk1D(n): #n here is the no. of steps that we require
    x = 5
    y = 5
    xposition = [x] 
    yposition = [y] 
    for i in range (1,n+1):
        step = np.random.uniform(0,1)
        if step < 0.5: # if step is less than 0.5 we move up               
            y += 1  #moving up in u direction
        if step > 0.5: # if step is greater than 0.5 we move down             
            y += -1 #moving down in y direction
        x += 1
        xposition.append(x)
        yposition.append(y)

    return [xposition,yposition]

Randwalk = Randomwalk1D(1000) #creating an object for the Randomwalk1D class and passing value of n as 100
plt.plot(Randwalk[0],Randwalk[1],'r-', label = "Randwalk1D") # 'r-' makes the color of the path red
plt.xlabel('tempo')
plt.ylabel('posicao')
plt.title("1-D Random Walks")
plt.show()