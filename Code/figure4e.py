################################################################
# Figure 4e - Heatmap for the prejudice-to-objectivity ratio q #
################################################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_fig4e'
# Running this function with the appropriate parameters will regenerate the associated heatmaps from Fig4e

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

## gen_fig4e ##
# Returns the heatmap over the prejudice-to-objectivity ratio (q)
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig4e(form):
    
    # Setting up sampling grid and plotting window
    q_values = np.linspace(0, 1, 100)
    theta_values = np.linspace(0, np.pi, 100)
    alpha_values = m*np.sin(theta_values)*gamma
    beta_values = m*np.cos(theta_values)
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
    
    q_averages = []
    theta_averages = []
    c_averages = []
    e_averages = []
    
    # Iterate
    for q in q_values:
        lam = 1 - (1-q)*(1-0.25)
        omega = q*(1-0.25)/lam
        for i in range(iterations):
            a,b,e = HK_prejudice_map(a, b, e, lam, omega, alpha_values, beta_values)
            if i >= (iterations - last):
                a_result[i - (iterations - last)] = a
                b_result[i - (iterations - last)] = b
                e_result[i - (iterations - last)] = e  
        for k in range(len(theta_values)):
            e_averages.append(np.average(e_result[:, k]))
            c_averages.append(np.average(np.concatenate((a_result[:, k], b_result[:, k]))))
            q_averages.append(q)
            theta_averages.append(theta_values[k])  
            
    if form == 'show':
        plt.figure()
        plt.subplot(211)
        plt.xlabel('Prejudice-to-Objectivity ratio (q)')
        plt.ylabel('Relative cost of mitigation (theta)')
        plt.scatter(q_averages, theta_averages, c=e_averages, vmin=0, vmax=5) 
        # Add colorbar (gradient legend)
        cbar = plt.colorbar()
        cbar.set_label('Pollutant level')
        plt.subplot(212)
        plt.xlabel('Prejudice-to-Objectivity ratio (q)')
        plt.ylabel('Relative cost of mitigation (theta)')
        plt.scatter(q_averages, theta_averages, c=c_averages, vmin=-1, vmax=1) 
        cbar = plt.colorbar()
        cbar.set_label('Level of mitigation')
        plt.show() 
    
    if form == 'export':
        plt.figure(figsize=(6.15, 6.15))
        plt.subplot(211)
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.yticks(np.arange(0, (np.pi + 0.1), step=(np.pi/4)))
        plt.scatter(q_averages, theta_averages, c=e_averages, vmin=0, vmax=5) 
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.subplot(212)
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.yticks(np.arange(0, (np.pi + 0.1), step=(np.pi/4)))
        plt.scatter(q_averages, theta_averages, c=c_averages, vmin=-1, vmax=1) 
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.savefig("Figure4e.png")         
         
## HK_prejudice_map ##
# Returns the next iteration in the model

def HK_prejudice_map(a, b, e, lam, omega, alpha, beta):
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
