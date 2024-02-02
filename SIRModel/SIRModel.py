import matplotlib.pyplot as plt
import numpy as np

def forward_euler(N, S, I, R, beta, gamma, t):
    """
    forward_euler computes one step of forward Euler 
    :param N: Total number of people in population
    :param S: Total number of susceptibles in population
    :param I: Total number of infected people in population
    :param R: Total number of recovered people in population
    :param beta: Transmission rate constant
    :param gamma: Recovery rate
    :param t: time step value
    :returns: new values for S, I, R (in that order)
    """
    # We are assuming that every interaction occurs at the exact same time
    # i.e. people who just got infected this time step are not included in people who can recover
    S_dot = -((beta * S * I)/float(N))
    I_dot = ((beta * S * I)/float(N)) - gamma * I
    R_dot = gamma * I
    return S + (S_dot * t), I + (I_dot * t), R + (R_dot * t)

# params to change
t = 0.1
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
        S, I, R = forward_euler(N, S, I, R, beta, gamma, t)
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

##################################
# Q3.d, beta=1, gamma=0.5, R_0=2 #
##################################
    
N = 1000
S = 999
I = 1
R = 0

S_arr = [S]
I_arr = [I]
R_arr = [R]
for i in range(iterations - 1):
    S, I, R = forward_euler(N, S, I, R, 1, 0.5, t)
    S_arr.append(S)
    I_arr.append(I)
    R_arr.append(R)

x = list(range(iterations))
fig, ax = plt.subplots()
ax.plot(x, S_arr, 'b', label='S')
ax.plot(x, I_arr, 'r', label='I')
ax.plot(x, R_arr, 'k', label='R')
ax.plot(x, [N*0.7968]*len(x), 'g--', label='r_inf')
ax.set_xlabel('Time')
ax.set_ylabel('People')
ax.set_title('SIR with beta=1, gamma=0.5, and r_inf=0.7968 (Kevin Ash)')

legend = ax.legend(loc='right', fontsize='large')

plt.savefig('SIR_r_inf.png')