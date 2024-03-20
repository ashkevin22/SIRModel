from scipy.stats import nbinom
import matplotlib.pyplot as plt

# Modifiable parameters
K = [1.0]
R0 = 3
G = 10
num_iters = 100_000

def simulate_epidemic(G, k, R0):
    '''
    Takes in generations, k-value, and R_0 value
    Returns true if epidemic dies, false otherwise
    '''
    prev_gen = 1
    max_infected = 0
    for _ in range(G):
        mean = R0
        variance = mean + (mean**2)/k
        p = mean/variance
        n = mean**2 / (variance - mean)
        # draw = nbinom.rvs(n=n,p=p)
        draws = nbinom.rvs(n=n,p=p,size=prev_gen)
        prev_gen = sum(draws)
        max_infected = max(prev_gen, max_infected)
        if prev_gen == 0:
            return True, max_infected
        if prev_gen >= 1000:
            return False, -1
    return False, -1

for k in K:
    num_dies = 0
    max_inf_arr = []
    for _ in range(num_iters):
        dies, max_inf = simulate_epidemic(G, k, R0)
        if dies:
            num_dies += 1
            max_inf_arr.append(max_inf)
    plt.hist(max_inf_arr, align='mid', rwidth=1, bins=range(1,7))
    plt.savefig('Branching_max_inf.png')
    print(f"k: {k}, q: {num_dies/float(num_iters)}")