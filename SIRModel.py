import matplotlib.pyplot as plt
import numpy as np

def forward_euler(N, S, I, R, beta, gamma):
    # We are assuming that every interaction occurs at the exact same time
    # i.e. people who just got infected this time step are not included in people who can recover
    new_infected = (beta * (S/float(N)) * (I/float(N))) * t
    new_recovered = (gamma * (I/float(N))) * t
    return S - new_infected, I + new_infected - new_recovered, R + new_recovered

# params to change
t = 100
iterations = 500
beta_arr = [1, 1.5, 2]
gamma = 0.5

for beta in beta_arr:
    N = 1000
    S = 999
    I = 1
    R = 0

    S_arr = [S]
    I_arr = [I]
    R_arr = [R]
    for i in range(iterations - 1):
        S, I, R = forward_euler(N, S, I, R, beta, gamma)
        S_arr.append(S)
        I_arr.append(I)
        R_arr.append(R)

    x = list(range(iterations))
    fig, ax = plt.subplots()
    ax.plot(x, S_arr, 'b', label='S')
    ax.plot(x, I_arr, 'r', label='I')
    ax.plot(x, R_arr, 'k', label='R')
    ax.set_xlabel('Time')
    ax.set_ylabel('People')
    ax.set_title('SIR with beta=' + str(beta) + ' and gamma=' + str(gamma) + ' (Kevin Ash)')

    legend = ax.legend(loc='right', fontsize='large')

    plt.savefig('SIR_beta_' + str(beta) + '.png')