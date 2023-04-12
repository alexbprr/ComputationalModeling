import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

def simulate(initials, propensities, stoichiometry, duration):
    """
    Run a simulation with given model.

    :param initials: List of initial population counts.
    :param propensities: List of functions that take population counts and give transition rates.
    :param stoichiometry: List of integers, how the population counts change per transition.
    :param duration: Maximum simulation time.
    :return: Two lists: The time points and population counts per time point.
    """

    # initial values
    times = [0.0]
    counts = [initials]

    # while finish time has not been reached
    while times[-1] < duration:
        # get current state
        state = counts[-1]

        # calculate rates with respective propensities
        rates = [prop(*state) for prop in propensities]

        # stop loop if no transitions available
        if all(r == 0 for r in rates):
            break

        # randomly draw one transition
        transition = random.choices(stoichiometry, weights=rates)[0]
        next_state = [a + b for a, b in zip(state, transition)]

        # draw next time increment from random exponential distribution
        # dt = math.log(1.0 / random.random()) / sum(weights)
        dt = random.expovariate(sum(rates))

        # append new values
        times.append(times[-1] + dt)
        counts.append(next_state)

    return times, counts

if __name__ == '__main__':
    N = 500  # whole population
    beta = 2  # transmission rate
    gamma = 0.5  # recovery rate
    t = 15  # duration

    I0 = 3
    initials = [N - I0, I0, 0]  # S, I, R

    propensities = [lambda s, i, r: beta * s * i / N,  # S -> I, Propensity: b * S(t) * I(t) / N
                    lambda s, i, r: gamma * i]  # I -> R Propensity: g * I(t)

    stoichiometry = [[-1, 1, 0],  # S -> I, Population change: S-1, I+1, R+0
                     [0, -1, 1]]  # I -> R Population change: S+0, I-1, R+1

    pdf = PdfPages('sir-gillespiev2.pdf')
    
    fig1 = plt.figure(figsize=(6,4))
    plt.title('Suscetíveis')
    plt.xlabel('tempo')
    plt.ylabel('população')
    ax1 = fig1.gca()

    fig2 = plt.figure(figsize=(6,4))
    plt.title('Infectados')
    plt.xlabel('tempo')
    plt.ylabel('população')
    ax2 = fig2.gca()

    fig3 = plt.figure(figsize=(6,4))
    plt.title('Recuperados')
    plt.xlabel('tempo')
    plt.ylabel('população')
    ax3 = fig3.gca()
    

    nRuns = 10
    k = 0
    colorindex = 0
    while k < nRuns:                  
        t = 15
        t, SIR = simulate(initials, propensities, stoichiometry, t)
        S, I, R = zip(*SIR)

        ax1.plot(t,S, color=(colorindex,colorindex,colorindex), label="S");
        ax2.plot(t,I, color=(colorindex,colorindex,colorindex), label="I")
        ax3.plot(t,R, color=(colorindex,colorindex,colorindex), label="R")

        k += 1
        colorindex += 0.1

    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    pdf.close()
