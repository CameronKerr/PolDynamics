#########################################
# Figure 3a-c - Cost-emissions heatmaps #
#########################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_fig3ac'
# Running this function with the appropriate parameters will regenerate the associated heatmap from Fig3a-c

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

## gen_fig3ac ##
# Returns the cost-emissions heatmap for given value of q
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig3ac(q, form):
    
    # Fixing remaining parameters
    lam = 1 - (1-q)*(1-0.9)
    omega = q*(1-0.9)/lam
    
    # Setting up sampling grid and averaging window
    # Setting up bifurcation diagram loop
    beta_values = np.linspace(-10, 10, 100)
    alpha_values = np.linspace(0, 2, 100)
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
    
    alpha_averages = []
    beta_averages = []
    c_averages = []
    e_averages = []
    
    # Iterate
    for alpha in alpha_values:
        for i in range(iterations):
            a, b, e = HK_prejudice_map(a, b, e, beta_values, alpha, lam, omega)
            if i >= (iterations - last):
                a_result[i - (iterations - last)] = a
                b_result[i - (iterations - last)] = b
                e_result[i - (iterations - last)] = e
        for k in range(len(beta_values)):
            e_averages.append(np.average(e_result[:, k]))
            c_averages.append(np.average(np.concatenate((a_result[:, k], b_result[:, k]))))
            alpha_averages.append(alpha)
            beta_averages.append(beta_values[k]) 
            
    if form == 'show':
        plt.figure()
        plt.suptitle('Proportion of prejudice = ' + str((omega*lam)/(1-lam*(1-omega))))
        plt.subplot(211)
        plt.ylabel('Pollutant emission rate')
        plt.scatter(beta_averages, alpha_averages, c=e_averages, vmin=0, vmax=10) 
        # Add colorbar (gradient legend)
        cbar = plt.colorbar()
        cbar.set_label('Pollutant level')
        plt.subplot(212)
        plt.xlabel('Net cost of mitigation')
        plt.ylabel('Pollutant emission rate')
        plt.scatter(beta_averages, alpha_averages, c=c_averages, vmin=-1, vmax=1) 
        cbar = plt.colorbar()
        cbar.set_label('Level of mitigation')
        plt.show()
        
    if form == 'export':
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.scatter(beta_averages, alpha_averages, c=e_averages, vmin=0, vmax=10)
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.savefig("Figure3heatmap_e.pdf") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.scatter(beta_averages, alpha_averages, c=c_averages, vmin=-1, vmax=1) 
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.savefig("Figure3heatmap_o.pdf")         

## HK_prejudice_map ##
# Returns the next iteration of the model

# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e, beta, alpha, lam, omega):
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