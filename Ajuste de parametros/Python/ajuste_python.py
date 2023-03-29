import numpy as np
from scipy.integrate import  solve_ivp
from scipy.optimize import differential_evolution
import math

path = 'data/'

def odeSystem(t, u, r, k):

    N = u[0]
    dN_dt = r*N*(1 - N/k)
 
    return [dN_dt]  

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    global data, reference_times
    dt = 0.01
    tfinal = 50
    times = np.arange(0,tfinal+dt,dt)

    N0 = x[0]
    u = [N0]

    r = x[1]
    k = x[2]
    params = (r, k)
    
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), u, t_eval=times, method='Radau')    

    u = results.y[0,:]
    error = 0
    sumobs = 0
    i = 0
    j = 0    
    for t in times:
        if isReferenceTime(reference_times,t):
            Ndata = data[i][1]
            error += (u[j] - Ndata)*(u[j] - Ndata) 
            sumobs += Ndata*Ndata
            i = i + 1
        j = j + 1

    error = math.sqrt(error/sumobs) #Erro norma 2                  
    return error 

if __name__ == "__main__":
    global data, reference_times 
    data = np.loadtxt(path+'data.csv', delimiter=',')
    reference_times = data[:,0]
    bounds = [
        (1, 200), (0.01, 1), (1,200)
    ]

    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, bounds, strategy='best1bin', maxiter=30,popsize=80,atol=10**(-3),tol=10**(-3),
                                     mutation=0.8, recombination=0.5,disp=True, workers=4)
    
    print(solucao.x)
    print(solucao.success)
    #saving the best offspring...
    np.savetxt('solucao_ajuste.txt',solucao.x, fmt='%.2f')        
    best=solucao.x
    error = solve(best)
    print("ERROR ", error)