import matplotlib.pyplot as plt
import numpy as np

def forward_euler_ANM(N, S, I, R, beta, gamma, t):
    """
    forward_euler_ANM computes one step of forward Euler for all or nothing model
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

def forward_euler_leaky(N, S, V, I, R, beta, gamma, VE, t):
    """
    forward_euler_leaky computes one step of forward Euler for leaky model
    :param N: Total number of people in population
    :param S: Total number of susceptibles in population
    :param V: Total number of vaccinated people in population
    :param I: Total number of infected people in population
    :param R: Total number of recovered people in population
    :param beta: Transmission rate constant
    :param gamma: Recovery rate
    :param t: time step value
    :returns: new values for S, V, I, R (in that order)
    """
    # We are assuming that every interaction occurs at the exact same time
    # i.e. people who just got infected this time step are not included in people who can recover
    S_dot = -((beta * S * I)/float(N))
    V_dot = -(((beta * V * I)/float(N)) * (1-VE))
    I_dot = ((beta * S * I)/float(N)) + (((beta * V * I)/float(N)) * (1-VE)) - gamma * I
    R_dot = gamma * I
    return S + (S_dot * t), V + (V_dot * t), I + (I_dot * t), R + (R_dot * t)

def all_or_nothing(t, iterations, beta, gamma):
    """
    all_or_nothing computes the all or nothing model, simplifies our main function
    :param t: Time step value
    :iterations: Total number of iterations for the simulation
    :param beta: Transmission rate constant
    :param gamma: Recovery rate
    """
    N = 300_000
    # S = 50% plus the people who were not effected by the vaccine
    S = 149_999 + (150_000 * 0.2)
    # V = 50% times the effecitveness of the vaccine
    V = 150_000 * 0.8
    I = 1
    R = 0

    S_arr = [S]
    I_arr = [I]
    R_arr = [R]
    V_arr = [V] * iterations
    for i in range(iterations - 1):
        S, I, R = forward_euler_ANM(N, S, I, R, beta, gamma, t)
        S_arr.append(S)
        I_arr.append(I)
        R_arr.append(R)
    return S_arr, V_arr, I_arr, R_arr

def leaky(t, iterations, beta, gamma):
    """
    leaky computes the leaky model, simplifies our main function
    :param t: Time step value
    :iterations: Total number of iterations for the simulation
    :param beta: Transmission rate constant
    :param gamma: Recovery rate
    """
    N = 300_000
    S = 149_999
    V = 150_000
    I = 1
    R = 0

    S_arr = [S]
    I_arr = [I]
    R_arr = [R]
    V_arr = [V]
    for i in range(iterations - 1):
        S, V, I, R = forward_euler_leaky(N, S, V, I, R, beta, gamma, 0.8, t)
        S_arr.append(S)
        I_arr.append(I)
        R_arr.append(R)
        V_arr.append(V)
    return S_arr, V_arr, I_arr, R_arr


# params to change
t = 0.1
iterations = 500
beta_arr = [3, 4, 5]
# gamma remain the same bc 14 day recovery for all variants
gamma = 1


for beta in beta_arr:
    S_arr_ANM, V_arr_ANM, I_arr_ANM, R_arr_ANM = all_or_nothing(t, iterations, beta, gamma)
    S_arr_leaky, V_arr_leaky, I_arr_leaky, R_arr_leaky = leaky(t, iterations, beta, gamma)
    x = list(range(iterations))
    fig, ax = plt.subplots()
    # sorry for anyone who has to read how I'm graphing these 
    # i hate Python and I plan on making that everyone else's problem
    total_infected_ANM = [S_arr_ANM[0]] * iterations
    total_infected_ANM = np.subtract(total_infected_ANM, S_arr_ANM)
    total_infected_leaky = [S_arr_leaky[0] + V_arr_leaky[0]] * iterations
    total_infected_leaky = np.subtract(total_infected_leaky, S_arr_leaky)
    total_infected_leaky = np.subtract(total_infected_leaky, V_arr_leaky)
    ax.plot(x, total_infected_ANM, 'r', label='ANM I')
    ax.plot(x, total_infected_leaky, 'g', label='Leaky I')
    ax.set_xlabel('Time')
    ax.set_ylabel('People')
    ax.set_title('SVIR total infections with beta=' + str(beta) + ' and gamma=' + str(gamma) + ' (Kevin Ash)')

    legend = ax.legend(loc='right', fontsize='large')

    plt.savefig('SVIR_total_infected_beta_' + str(beta) + '.png')

    ax.clear()

    ax.plot(x, I_arr_leaky, 'g', label='Leaky I')
    ax.plot(x, I_arr_ANM, 'r', label='ANM I')
    ax.plot(x, [max(I_arr_ANM)] * len(x), 'b', label='ANM peak I')
    ax.plot(x, [max(I_arr_leaky)] * len(x), 'k', label='Leaky peak I')

    ax.set_xlabel('Time')
    ax.set_ylabel('People')
    ax.set_title('SVIR infections with beta=' + str(beta) + ' and gamma=' + str(gamma) + ' (Kevin Ash)')


    legend = ax.legend(loc='right', fontsize='large')

    plt.savefig('SVIR_infected_beta_' + str(beta) + '.png')

    print("infected ratio for beta=" + str(beta) + " is: " + str(max(total_infected_ANM) / max(total_infected_leaky)))
