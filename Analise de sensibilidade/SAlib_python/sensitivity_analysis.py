#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from scipy.integrate import solve_ivp
from tqdm import tqdm
from uqtools import *
from SALib.sample import saltelli
from SALib.analyze import sobol

#TEMPO
nsteps = 100
t=range(nsteps+1)

def odeSystem(u, t, r, k):
    global model_args
    N = u[0]
    dNdt = r*N*(1 - N/k)
    return [dNdt]

def condicao_inicial():
    N = 90
    u = (N)
    return u

def eval_model(params):
    global model_args
    r = params[0]
    k = params[1]

    model_args = (r, k)
    u0 = condicao_inicial()
    y = odeint(odeSystem, u0, t, args=model_args)

    return y[:,0]

if __name__ == "__main__":
    opt_uq = False
    opt_loocv = False
    opt_save_evals = True

    opt_sobol_salib = True
    opt_sobol_pce = False

    label_param = ['r', 'k']

    # Sensitivity analysis
    if(opt_sobol_salib):
        min_bound = 0.5  
        max_bound = 1.5 

        #Valores iniciais (valores de referências)
        r = 0.1
        k = 100 
        model_params = (r, k)
        vbounds = [] #Vetor de limites inferior e superior pra cada parâmetro
        for i in range(len(model_params)):
            vbounds.append([model_params[i]*min_bound, model_params[i]*max_bound])

        #Pra conseguir gerar a amostra, precisa dessa estrutura
        Nvars = 2 #numero de parametros
        # Define the model inputs
        problem = {
            'num_vars': Nvars,
            'names': label_param, 
            'bounds': vbounds
        }

        # Generate samples
        nsobol = 1024 #numeros de amostras

        param_values = saltelli.sample(problem, nsobol, calc_second_order=False)

        # Run model (example)
        k = 0
        nexec = np.shape(param_values)[0] #numero de linhas da matriz com as amostras
        lpts = range( nexec )
        #Indice sobol principal
        sm_N = []
        print("evaluating samples for SA: ")
        #tqdm exibe a barra de progresso
        #Uma amostra é um conjunto de valores de parâmetros
        for i in tqdm(lpts,bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}'):
            s = param_values[k,:] #samp[k] #acessar todas as colunas d alinha k
            s = np.array(s)
            k = k+1
            N = eval_model(s) #v, P #resolve as EDOs para a amostra na linha k
            #ignora o 1o passo de tempo porque é a condição inicial e não tem variancia (pode levar a um NaN)
            sm_N.append(N[1:])
            
        sm_N = np.array(sm_N, dtype="object")

        #Indice S1 para cada passo de tempo
        Si_N = np.zeros((nsteps,Nvars))

        # Perform analysis
        for i in range(nsteps): #para cada valor no tempo
            sob_N = sobol.analyze(problem, sm_N[:,i], calc_second_order=False, 
                                print_to_console=False)
            Si_N[i,:] = sob_N['S1']
            
        print("salvando arquivos Sobol")
        if(opt_save_evals):
            np.savetxt('ouput_sobol_v.txt',Si_N)

        plt.figure()
        plt.title("N")
        plot_sensitivity_mc(plt, t[1:], Si_N.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_N.pdf')
# Fim
print('done')
