################################################
# Figure 4a-c - Bifurcation diagram over theta #
################################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_fig4ac'
# Running this function with the appropriate parameters will regenerate the associated bifurcation diagrams from Fig4a-c

# Load required libraries
import numpy as np
import matplotlib.pyplot as plt

# Fixed parameters used in figures
tau = 1.5
p = 0.1
r = 0.2
m = 5 # Magnitude in the (beta, alpha) parameter space
gamma = 0.2
a0 = 1
b0 = -1

## gen_fig4ac ##
# Returns the bifurcation diagram over the relative cost of mitigation (theta)
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig4ac(q, form):
    # Setting up influence coefficients
    lam = 1 - (1-q)*(1-0.9)
    omega = q*(1-0.9)/lam
    
    # Setting up sampling grid and plotting window
    theta_values = np.linspace(0, np.pi, 10000)
    iterations = 1000
    last = 100
    # Initialize
    a = np.zeros(len(theta_values))
    b = np.zeros(len(theta_values))
    e = np.zeros(len(theta_values))
    a[:] = a0
    b[:] = b0
    e[:] = 0
    a_result = np.empty((last, len(theta_values)))
    b_result = np.empty((last, len(theta_values)))
    e_result = np.empty((last, len(theta_values)))
    
    # Iterate
    for i in range(iterations):
        alpha_values = m*np.sin(theta_values)*gamma
        beta_values = m*np.cos(theta_values)
        a, b, e = HK_prejudice_map(a, b, e, alpha_values, beta_values, lam, omega)
    
        if i >= (iterations - last):
            a_result[i - (iterations - last)] = a
            b_result[i - (iterations - last)] = b
            e_result[i - (iterations - last)] = e
   
    if form == 'show':
        
        plt.figure()
        plt.suptitle('Proportion of prejudice = ' + str((omega*lam)/(1-lam*(1-omega))))
        plt.subplot(311)
        plt.ylim(-0.1, 5.1)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('Environmental pollution')
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, e_result.T, ',k', alpha = 0.25)
        
        
        plt.subplot(312)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('a - opinion')
        plt.ylim(-1.1, 1.1)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
        
        plt.subplot(313)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('b - opinion')
        plt.ylim(-1.1, 1.1)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
        plt.show()     
    
    if form == 'export':
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-0.1, 5.1)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, e_result.T, ',k', alpha = 0.15)
        plt.savefig("Figure4bif_e.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure4bif_a.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4))        
        plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure4bif_b.png")         
         
## HK_prejudice_map ##
# Returns the next iteration in the model

def HK_prejudice_map(a, b, e, alpha, beta, lam, omega):
    e_next = (1 - gamma)*e + 1/2*alpha*(1 - (p*a + (1-p)*b))
    # Initialize empty iteration
    a_next = np.empty(len(a))
    b_next = np.empty(len(b))
    
    for i in range(len(a)):
        if abs(a[i] - b[i]) > r:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega)*a[i] + omega*a0) + (1-lam)*np.sign(e[i] - beta[i]))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega)*b[i] + omega*b0) + (1-lam)*np.sign(e[i] - beta[i]))
        else:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega)*(p*a[i] + (1-p)*b[i]) + omega*a0) + (1-lam)*np.sign(e[i] - beta[i]))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega)*(p*a[i] + (1-p)*b[i]) + omega*b0) + (1-lam)*np.sign(e[i] - beta[i]))
            
    return a_next, b_next, e_next

