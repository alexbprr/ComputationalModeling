#Referencia: https://medium.com/geekculture/runge-kutta-numerical-integration-of-ordinary-differential-equations-in-python-9c8ab7fb279c
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 

# Runge-Kutta (RK4) Numerical Integration for System of First-Order Differential Equations

def ode_system(t, y):
    """
    system of first order differential equations
    t: discrete time step value
    y: state vector (vetor das variaveis/populacoes do modelo)
    """
    r = 0.2
    k = 100
    #dydt = y[0]*np.exp(a*t)
    dndt = r*y[0]*(1 - y[0]/k)
    #dhdt = 0.1*y[0] - 0.02*y[1]*y[0] #EDO da presa (h)
    #dpdt = 0.02*y[1]*y[0] -0.08*y[1] #EDO do predador (p)
    return np.array([dndt])


def rk4(func, tk, _yk, _dt=0.01, **kwargs):
    """
    single-step fourth-order numerical integration (RK4) method
    func: system of first order ODEs
    tk: current time step
    _yk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """

    # evaluate derivative at several stages within time interval
    f1 = func(tk, _yk, **kwargs)
    f2 = func(tk + _dt / 2, _yk + (f1 * (_dt / 2)), **kwargs)
    f3 = func(tk + _dt / 2, _yk + (f2 * (_dt / 2)), **kwargs)
    f4 = func(tk + _dt, _yk + (f3 * _dt), **kwargs)

    # return an average of the derivative over tk, tk + dt
    return _yk + (_dt / 6) * (f1 + (2 * f2) + (2 * f3) + f4)


# ==============================================================
dt = 0.01  #Passo de tempo
tfinal = 200 
time = np.arange(0, tfinal + dt, dt) #Valores no tempo nos quais a EDO será calculada 

# second order system initial conditions [y1, y2] at t = 0
y0 = np.array([10]) #Condicao inicial 

# ==============================================================
# propagate state

# simulation results
state_history = []

# initialize yk
yk = y0

# intialize time
t = 0
# Resolve o sistema de EDOs (loop temporal)
for t in time:
    state_history.append(yk)
    yk = rk4(ode_system, t, yk, dt)

# convert list to numpy array
state_history = np.array(state_history)
print(f'y evaluated at time t = {t} seconds: {yk[0]}')

#save results to csv file 
df = pd.DataFrame(state_history, columns = ['N'])
df.insert(0, 'Time', time)
df.to_csv('ode_results.csv', float_format='%.5f', sep=',') #salva os valores com cinco casas decimais 
# ==============================================================
# plot history

fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.plot(time, state_history[:, 0])
#ax.plot(time, state_history[:, 1])
ax.set(xlabel='time (days)', ylabel='[Y]', title='Modelo logistico')
plt.legend(['N'], loc='best')
ax.grid()
fig.savefig('logistico.svg', format='svg')
plt.show()