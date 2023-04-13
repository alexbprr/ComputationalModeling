#Referencia https://www.codespeedy.com/random-walk-program-in-python/
import numpy as np 
import matplotlib.pyplot as plt 
import random 

n = 100 #n is the number of steps(increase in the value of n increses the compelxity of graph) 
x = np.zeros(n) # x and y are arrays which store the coordinates of the position 
y = np.zeros(n) 
direction=["NORTH","SOUTH","EAST","WEST"] # Assuming the four directions of movement.

for i in range(1, n): 
    step = random.choice(direction) #Randomly choosing the direction of movement. 
    if step == "EAST": #updating the direction with respect to the direction of motion choosen.
        x[i] = x[i - 1] + 1
        y[i] = y[i - 1] 
    elif step == "WEST": 
        x[i] = x[i - 1] - 1
        y[i] = y[i - 1] 
    elif step == "NORTH": 
        x[i] = x[i - 1] 
        y[i] = y[i - 1] + 1
    else: 
        x[i] = x[i - 1] 
        y[i] = y[i - 1] - 1

plt.title("Random Walk 2-D")
plt.plot(x, y) #plotting the walk.
plt.show()