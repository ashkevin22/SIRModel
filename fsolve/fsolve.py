import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import fsolve


def g_r_inf(r_inf, r_0):
    return 1 - (math.e ** ((-r_0)*(r_inf)))

def f_r_inf(r_inf):
    return r_inf

r_inf = np.arange(0.,0.5,0.01)
R_0 = [0.9,1.0,1.1,1.2]

for r_0 in R_0:
    # referece: https://glowingpython.blogspot.com/2011/05/hot-to-find-intersection-of-two.html
    intersect = fsolve(lambda x : g_r_inf(x, r_0) - f_r_inf(x), 1)
    print(intersect)

    g_r_inf_values = [g_r_inf(x, r_0) for x in r_inf]
    f_r_inf_values = [f_r_inf(x) for x in r_inf]
    # g_r_inf_values = map(g_r_inf, r_inf, )
    # f_r_inf_values = map(f_r_inf, r_inf)

    plt.cla()
    plt.xlabel('r_infinity')
    plt.title('fsolve with R_0='+str(r_0))
    plt.plot(r_inf, f_r_inf_values, 'k', r_inf, g_r_inf_values, 'r')
    plt.plot(intersect, f_r_inf(intersect), 'bo', markersize=10)
    plt.savefig('fsolve_r_0_'+str(r_0)+'.png')