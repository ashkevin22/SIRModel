import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def forward_euler(S, I, R, P, C, omega, gamma, t):
    """
    forward_euler computes one step of forward Euler 
    :param S: Vector with total number of susceptibles in each population
    :param I: Vector with total number of infected people in each population
    :param R: Vector with total number of recovered people in each population
    :param P: Vector with probability of getting infected given contact for each pop.
    :param C: Contact matrix
    :param omega: Vector with population proportion for each population
    :param gamma: Vector of recovery rates for each population
    :param t: time step value

    :returns: new values for S, I, and R (in that order)
    """
    # \dot{s} = -p(D_S * C * D^{-1}_w)i
    # \dot{i} = p(D_S * C * D^{-1}_w) - \gamma * i
    # Separating the equation out bc it is totally incomprehensible in one line
    step_1 = np.dot(np.diag(S), np.diag(P))
    step_2 = np.dot(step_1, C)
    multiplier = np.dot(step_2, np.linalg.inv(np.diag(omega)))

    S_dot = - np.dot(multiplier, I)
    I_dot = np.dot(multiplier, I) - np.dot(np.diag(gamma), I)
    R_dot = np.dot(np.diag(gamma), I)

    return S + (S_dot * t), I + (I_dot * t), R + (R_dot * t)

# params to change
iterations = 400
t = 0.01
S = np.array([0.999, 0.999, 0.999, 0.999])
I = np.array([0.001, 0.001, 0.001, 0.001])
R = np.array([0] * 4)
P = np.array([1, 2, 3, 4])
# python wizardry (I understand it, I just think it's dumb)
omega = np.array([0.25] * 4)
gamma = np.array([3] * 4)
# :(
C = np.array([[float(9)/20] * 4] * 4)

S_arr = S
I_arr = I
R_arr = R
p_bar_arr = []
for i in range(iterations - 1):
    S, I, R = forward_euler(S, I, R, P, C, omega, gamma, t)
    S_arr = np.vstack((S_arr, S))
    I_arr = np.vstack((I_arr, I))
    R_arr = np.vstack((R_arr, R))
    
    numerator = 0
    denom = 0
    for i, s in enumerate(S):
        numerator += s * P[i]
        denom += s
    p_bar_arr.append(numerator/denom)


#Thanks to Parker for helping me w/ the graphing bc I could not figure it out for the life of me
num_lines = len(I_arr[0])
t_arr = np.linspace(0, iterations * t, iterations)

colors = sns.color_palette("Reds", n_colors=num_lines)

# Actually plotting
for ind in range(num_lines):
    i = I_arr[:, ind]
    plt.plot(t_arr, i, label=f"Group {ind + 1}", color=colors[ind])

plt.xlabel("Time")
plt.ylabel("Proportion of population")
plt.title("Infected Proportions for Each Population (Kevin Ash)")
plt.legend()
plt.savefig('SIRModelPopStruct.png')

plt.cla()

colors = sns.color_palette("Blues", n_colors=num_lines)

# Actually plotting
for ind in range(num_lines):
    s = S_arr[:, ind]
    plt.plot(t_arr, s, label=f"Group {ind + 1}", color=colors[ind])
plt.xlabel("Time")
plt.ylabel("Proportion of population")
plt.title("Susceptible Proportions for Each Population (Kevin Ash)")
plt.legend()
plt.savefig('SIRModelPopStructSus.png')

plt.cla()

plt.plot(t_arr[1:], p_bar_arr, label="p_bar", color="black")
plt.xlabel("Time")
plt.ylabel("Relative Susceptibility")
plt.title("Relative Susceptibility Among Susceptibles (Kevin Ash)")
plt.legend()
plt.savefig('SIRModelPopStructRelSus.png')