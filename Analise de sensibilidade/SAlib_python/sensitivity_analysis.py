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

def odeSystem(u, t, b, k, a, p, k2, s, mT, lambdaAC, mTreg, betaAC, betaC, mAC, mC, dT, aTreg):
    global model_args
    E, I, V, T, Treg, AC, C, D = u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7]
    #b, k, a, p, k2, s, mT, lambdaAC, mTreg, betaAC, betaC, mAC, mC, E0, dT, aTreg =
    # model_args[0], model_args[1], model_args[2], model_args[3]...
    dD_dt = dT*(T + I)
    dE_dt = -b*E*V
    dI_dt = b*E*V - k*I*T - a*I
    dV_dt = p*a*I - k2*V
    dT_dt = (s*(1+C) / (1+lambdaAC*AC)) - aTreg*D*T - mT*T
    dTreg_dt = aTreg*D*T - mTreg*Treg
    dAC_dt = betaAC*Treg - mAC*AC
    dC_dt = ((betaC*T)/(1+lambdaAC*AC)) - mC*C

    return [dE_dt, dI_dt, dV_dt, dT_dt, dTreg_dt, dAC_dt, dC_dt, dD_dt]

"""
def solveOde(t, y):
    global model_args
    print(y)
    return odeSystem(y, *model_args)
"""

def condicao_inicial():
    E = 5000
    I = 0
    V = 10
    T = 10 #5
    Treg = 5 #10
    AC = 0
    C = 0
    D = 0
    u = (E, I, V, T, Treg, AC, C, D)
    return u

def eval_model(params):
    global model_args
    b = params[0]
    k = params[1]
    a = 0.1
    p = 5
    k2 = 0.4
    s = 0.2
    mT = 0.005
    lambdaAC = params[2]
    mTreg = 0.05
    betaAC = params[3]
    betaC = 0.5
    mAC = 0.5
    mC = 0.1
    dT = 0.005
    aTreg = params[4]
    model_args = (b, k, a, p, k2, s, mT, lambdaAC, mTreg, betaAC, betaC, mAC, mC, dT, aTreg)
    P0 = condicao_inicial()
    y = odeint(odeSystem, P0, t, args=model_args)

    return y[:,0], y[:,1], y[:,2], y[:,3], y[:,4], y[:,5], y[:,6], y[:,7] 

if __name__ == "__main__":
    opt_uq = False
    opt_loocv = False
    opt_save_evals = True

    opt_sobol_salib = True
    opt_sobol_pce = False

    #label_param = ('b', 'k', 'a', 'p', 'k2', 's', 'mT', 'lambdaAC', 'mTreg', 'betaAC', 'betaC', 'mAC', 'mC', 'E0', 'dT', 'aTreg')
    label_param = ['b', 'k','lambdaAC', 'betaAC','aTreg']

    # Sensitivity analysis
    if(opt_sobol_salib):
        min_bound = 0.5  
        max_bound = 1.5 

        #Valores iniciais (valores de referências)
        b = 0.00002
        k = 0.1 # 0.2
        a = 0.1
        p = 5
        k2 = 0.4
        s = 0.2
        mT = 0.005
        lambdaAC = 0.3
        mTreg = 0.05
        betaAC = 0.05
        betaC = 0.5
        mAC = 0.5
        mC = 0.1
        dT = 0.005 
        aTreg = 0.01
        #k = 100 # Valor de referência para gerar as amostras
        model_params = (b, k, lambdaAC, betaAC, aTreg)
        vbounds = [] #Vetor de limites inferior e superior pra cada parâmetro
        for i in range(len(model_params)):
            vbounds.append([model_params[i]*min_bound, model_params[i]*max_bound])

        #Pra conseguir gerar a amostra, precisa dessa estrutura
        Nvars = 5 #numero de parametros
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
        #sm_N = []
        #Analisar para V, T, Treg e AC
        sm_V = []
        sm_T = []
        sm_Treg = []
        sm_AC = []
        sm_C = []
        sm_I = []
        print("evaluating samples for SA: ")
        #tqdm exibe a barra de progresso
        #Uma amostra é um conjunto de valores de parâmetros
        for i in tqdm(lpts,bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}'):
        #for s in samples.T:
            s = param_values[k,:] #samp[k] #acessar todas as colunas d alinha k
            s = np.array(s)
            k = k+1
            #RECLAMOU AO COLOCAR SÓ V, T, Treg, AC
            E, I, V, T, Treg, AC, C, D = eval_model(s) #v, P #resolve as EDOs para a amostra na linha k
            #Plotar curva da população  (if (k % 10) == 0)
            #ignora 1o passo de tempo porque é a condição inicial e não tem variancia, levando a um NaN
            sm_V.append(V[1:])
            sm_T.append(T[1:])
            sm_Treg.append(Treg[1:])
            sm_AC.append(AC[1:])
            sm_C.append(C[1:])
            sm_I.append(I[1:])

        sm_V = np.array(sm_V, dtype="object")
        sm_T = np.array(sm_T, dtype="object")
        sm_Treg = np.array(sm_Treg, dtype="object")
        sm_AC = np.array(sm_AC, dtype="object")
        sm_C = np.array(sm_C, dtype="object")
        sm_I = np.array(sm_I, dtype="object")


        #Indice S1 para cada passo de tempo
        Si_V = np.zeros((nsteps,Nvars))
        Si_T = np.zeros((nsteps,Nvars))
        Si_Treg = np.zeros((nsteps,Nvars))
        Si_AC = np.zeros((nsteps,Nvars))
        Si_C = np.zeros((nsteps,Nvars))
        Si_I = np.zeros((nsteps,Nvars))

        # Perform analysis
        for i in range(nsteps): #para cada valor no tempo
            sob_V = sobol.analyze(problem, sm_V[:,i], calc_second_order=False, print_to_console=False)
            sob_T = sobol.analyze(problem, sm_T[:,i], calc_second_order=False, print_to_console=False)
            sob_Treg = sobol.analyze(problem, sm_Treg[:,i], calc_second_order=False, print_to_console=False)
            sob_AC = sobol.analyze(problem, sm_AC[:,i], calc_second_order=False, print_to_console=False)
            sob_C = sobol.analyze(problem, sm_C[:,i], calc_second_order=False, print_to_console=False)
            sob_I = sobol.analyze(problem, sm_I[:,i], calc_second_order=False, print_to_console=False)
            #Si_N[i,:] = sob_N['S1']
            Si_V[i,:] = sob_V['S1']
            Si_T[i,:] = sob_T['S1']
            Si_Treg[i,:] = sob_Treg['S1']
            Si_AC[i,:] = sob_AC['S1']
            Si_C[i,:] = sob_C['S1']
            Si_I[i,:] = sob_I['S1']

        print("salvando arquivos Sobol")
        if(opt_save_evals):
            np.savetxt('ouput_sobol_v.txt',Si_V)
            np.savetxt('ouput_sobol_T.txt',Si_T)
            np.savetxt('ouput_sobol_Treg.txt',Si_Treg)
            np.savetxt('ouput_sobol_AC.txt',Si_AC)
            np.savetxt('ouput_sobol_C.txt',Si_AC)
            np.savetxt('ouput_sobol_I.txt',Si_I)

        # plt.title("N")
        # plot_sensitivity_mc(plt, t[1:], Si_N.T, label_param)
        # plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        # plt.tight_layout()
        # plt.savefig('output_sens_N.pdf')

        plt.figure()
        plt.title("V")
        plot_sensitivity_mc(plt, t[1:], Si_V.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_V.pdf')

        plt.figure()
        plt.title("T")
        plot_sensitivity_mc(plt, t[1:], Si_T.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_T.pdf')

        plt.figure()
        plt.title("Treg")
        plot_sensitivity_mc(plt, t[1:], Si_Treg.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_Treg.pdf')

        plt.figure()
        plt.title("AC")
        plot_sensitivity_mc(plt, t[1:], Si_AC.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_AC.pdf')

        plt.figure()
        plt.title("C")
        plot_sensitivity_mc(plt, t[1:], Si_C.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_C.pdf')

        plt.figure()
        plt.title("I")
        plot_sensitivity_mc(plt, t[1:], Si_I.T, label_param)
        plt.legend(loc=(1.04,0),prop={'size': 8}, ncol=2)
        plt.tight_layout()
        plt.savefig('output_sens_I.pdf')


# Fim
print('done')
