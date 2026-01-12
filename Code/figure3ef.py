########################################################
# Figure 3ef - Cost and emissions bifurcation diagrams #
########################################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_fig3e' and 'gen_fig3f'
# Running this function with the appropriate parameters will regenerate the associated bifurcation diagrams from Fig3ef

# Load required libraries
import numpy as np
import matplotlib.pyplot as plt

# Fixed parameters used in figures
tau = 1.5
p = 0.1
r = 0.2
gamma = 0.2
a0 = 1
b0 = -1
lam = 0.95
omega = 0.05/0.95

## gen_fig3e ##
# Returns the bifurcation diagram over beta
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig3e(form):
    # Fixing remaining parameters 
    alpha = 1
    
    # Setting up sampling grid and plotting window
    beta_values = np.linspace(-10, 10, 10000)
    iterations = 1000
    last = 100
    # Initialize
    a = np.zeros(len(beta_values))
    b = np.zeros(len(beta_values))
    e = np.zeros(len(beta_values))
    a[:] = a0
    b[:] = b0
    e[:] = 0
    a_result = np.empty((last, len(beta_values)))
    b_result = np.empty((last, len(beta_values)))
    e_result = np.empty((last, len(beta_values)))
    
    # Iterate
    for i in range(iterations):
        a, b, e = HK_prejudice_map_beta(a, b, e, beta_values, alpha)
        if i >= (iterations - last):
            a_result[i - (iterations - last)] = a
            b_result[i - (iterations - last)] = b
            e_result[i - (iterations - last)] = e
        
    if form == 'show':
        plt.figure()
        plt.subplot(311)
        plt.xlabel('Net cost of mitigation')
        plt.ylabel('Environmental pollution')
        plt.ylim(-0.1, 5.1)        
        plt.plot(beta_values, e_result.T, ',k', alpha = 0.25)
        
        plt.subplot(312)
        plt.xlabel('Net cost of mitigation')
        plt.ylabel('a - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(beta_values, a_result.T, ',k', alpha = 0.25)        
        
        plt.subplot(313)
        plt.xlabel('Net cost of mitigation')
        plt.ylabel('b - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(beta_values, b_result.T, ',k', alpha = 0.25)
        
        plt.show()
    
    if form == 'export':
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-0.1, 5.1)
        plt.plot(beta_values, e_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3e_e.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(beta_values, a_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3e_a.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(beta_values, b_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3e_b.png")         
        
        
## HK_prejudice_map_beta ##
# Returns the next iteration of the model for the beta bifurcation diagram

# Map returns the next iteration in the trajectory 
def HK_prejudice_map_beta(a, b, e, beta, alpha):
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

## gen_fig3f ##
# Returns the bifurcation diagram over alpha
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig3f(form):
    # Fixing remaining parameters 
    beta = 3
    
    # Setting up sampling grid and plotting window
    alpha_values = np.linspace(0, 2, 10000)
    iterations = 1000
    last = 100
    # Initialize
    a = np.zeros(len(alpha_values))
    b = np.zeros(len(alpha_values))
    e = np.zeros(len(alpha_values))
    a[:] = a0
    b[:] = b0
    e[:] = 0
    a_result = np.empty((last, len(alpha_values)))
    b_result = np.empty((last, len(alpha_values)))
    e_result = np.empty((last, len(alpha_values)))
    
    # Iterate
    for i in range(iterations):
        a, b, e = HK_prejudice_map_alpha(a, b, e, beta, alpha_values)
        if i >= (iterations - last):
            a_result[i - (iterations - last)] = a
            b_result[i - (iterations - last)] = b
            e_result[i - (iterations - last)] = e
            
    if form == 'show':
        
        plt.figure()
        plt.subplot(311)
        plt.xlabel('Rate of emissions')
        plt.ylabel('Environmental pollution')
        plt.ylim(-0.1, 5.1)        
        plt.plot(alpha_values, e_result.T, ',k', alpha = 0.25)
        
        plt.subplot(312)
        plt.xlabel('Rate of emissions')
        plt.ylabel('a - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(alpha_values, a_result.T, ',k', alpha = 0.25)
        
        plt.subplot(313)
        plt.xlabel('Rate of emissions')
        plt.ylabel('b - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(alpha_values, b_result.T, ',k', alpha = 0.25)
        plt.show()  
        
    if form == 'export':
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-0.1, 5.1)
        plt.plot(alpha_values, e_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3f_e.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(alpha_values, a_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3f_a.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(alpha_values, b_result.T, ',k', alpha = 0.25)
        plt.savefig("Figure3f_b.png")           
        
## HK_prejudice_map_alpha ##
# Returns the next iteration of the model for the alpha bifurcation diagram

# Map returns the next iteration in the trajectory 
def HK_prejudice_map_alpha(a, b, e, beta, alpha):
    e_next = (1 - gamma)*e + 1/2*alpha*(1 - (p*a + (1-p)*b))
    # Initialize empty iteration
    a_next = np.empty(len(a))
    b_next = np.empty(len(b))
    
    for i in range(len(a)):
        if abs(a[i] - b[i]) > r:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega)*a[i] + omega*a0) + (1-lam)*np.sign(e[i] - beta))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega)*b[i] + omega*b0) + (1-lam)*np.sign(e[i] - beta))
        else:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega)*(p*a[i] + (1-p)*b[i]) + omega*a0) + (1-lam)*np.sign(e[i] - beta))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega)*(p*a[i] + (1-p)*b[i]) + omega*b0) + (1-lam)*np.sign(e[i] - beta))
            
    return a_next, b_next, e_next