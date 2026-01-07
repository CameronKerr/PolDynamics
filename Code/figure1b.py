#####################################
# Figure 1b - Numerical simulations #
#####################################

## INSTRUCTIONS ##
# Run this script to load in the required function 'gen_fig1b'
# Running 'gen_fig1b' with the appropriate parameters will regenerate the associated numerical trajectory from Fig1b

# Load required libraries
import numpy as np
import matplotlib.pyplot as plt

# Fixed parameters used in figures
tau = 1.5
p = 0.5
r = 0.2
gamma = 0.2
a0 = 0.5
b0 = -0.5

## gen_fig1 ##
# Returns the numerical trajectory for the model
# To regenerate a specific plot, use the values for alpha, beta, lambda, and omega listed in Fig 1b
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig1b(alpha, beta, lam, omega, form):
    # Get trajectory
    iterations = 50
    a = [a0]
    b = [b0]
    e = [0]
    # Iterate
    for i in range(iterations):
        a_next, b_next, e_next = HK_prejudice_map(a[-1], b[-1], e[-1], alpha, beta, lam, omega)
        a.append(a_next)
        b.append(b_next)
        e.append(e_next)
    
    # 'show' plotting
    if form == 'show':
        plt.subplot(211)
        plt.ylim(0, alpha/gamma+1)
        plt.xlim(0, iterations + 1)
        plt.ylabel('Environmental pollution - e(t)')
        plt.plot(list(range(iterations + 1)), e, 'green', alpha = 0.5)
        
        plt.subplot(212)
        plt.xlabel('Time')
        plt.ylabel('Opinion - a(t), b(t)')
        plt.ylim(-1.1, 1.1)
        plt.xlim(0, iterations + 1)
        plt.plot(list(range(iterations + 1)), a, 'blue', alpha=0.5)
        plt.plot(list(range(iterations + 1)), b, 'red', alpha=0.5)
        plt.show() 
    
    # 'export' plotting
    if form == 'export':
        plt.figure(figsize=(4, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-0.1, alpha/gamma + 1)
        plt.xlim(0, 50)
        plt.plot(list(range(iterations + 1)), e, 'green', alpha=0.5)
        plt.savefig('Figure1b_e.pdf')
        
        plt.figure(figsize=(4, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.xlim(0, 50)
        plt.plot(list(range(iterations + 1)), a, 'blue', alpha=0.5)
        plt.plot(list(range(iterations + 1)), b, 'red', alpha=0.5)
        plt.savefig('Figure1b_o.pdf')        
        
## HK_prejudice_map ##
# Returns the next iteration of the model

def HK_prejudice_map(a, b, e, alpha, beta, lam, omega):
    e_next = (1 - gamma)*e + 1/2*alpha*(1 - (p*a + (1-p)*b))
    
    if abs(a-b) > r:
        a_next = (1 - 1/tau)*a + 1/tau*(lam*((1-omega)*a + omega*a0) + (1-lam)*np.sign(e - beta))
        b_next = (1 - 1/tau)*b + 1/tau*(lam*((1-omega)*b + omega*b0) + (1-lam)*np.sign(e - beta))
    else:
        a_next = (1 - 1/tau)*a + 1/tau*(lam*((1-omega)*(p*a + (1-p)*b) + omega*a0) + (1-lam)*np.sign(e - beta))
        b_next = (1 - 1/tau)*b + 1/tau*(lam*((1-omega)*(p*a + (1-p)*b) + omega*b0) + (1-lam)*np.sign(e - beta))    
            
    return a_next, b_next, e_next