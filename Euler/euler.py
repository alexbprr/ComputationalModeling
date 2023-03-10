import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 

# Euler Numerical Integration for System of First-Order Differential Equations
def ode_system(t, y):
    """
    system of first order differential equations
    t: discrete time step value
    y: state vector (vetor das variaveis/populacoes do modelo)
    """
    r = 0.2
    dndt = r*y[0]
    return np.array([dndt])


def euler(func, tk, _yk, _dt=0.01, **kwargs):
    """
    single-step euler method
    func: system of first order ODEs
    tk: current time step
    _yk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    return _yk + func(tk, _yk, **kwargs)*_dt


# ==============================================================
dt = 0.01  #Passo de tempo
tfinal = 100 
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
    yk = euler(ode_system, t, yk, dt)

# convert list to numpy array
state_history = np.array(state_history)
print(f'y evaluated at time t = {t} seconds: {yk[0]}')

# save results to csv file 
df = pd.DataFrame(state_history, columns = ['N'])
df.insert(0, 'Time', time)
df.to_csv('ode_results.csv', float_format='%.5f', sep=',') 
# ==============================================================
# plot history

fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
ax.plot(time, state_history[:, 0])
ax.set(xlabel='time (days)', ylabel='[Y]', title='Modelo exponencial')
plt.legend(['N'], loc='best')
ax.grid()
fig.savefig('exponencial.svg', format='svg')
plt.show()