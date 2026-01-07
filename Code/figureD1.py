################################################################
# Figure D.1 - Bifurcation diagram over theta with predictions #
################################################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_figd1'
# Running this function with the appropriate parameters will regenerate the associated bifurcation diagrams from FigD.1

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

## gen_figd1 ##
# Returns the bifurcation diagram over the relative cost of mitigation (theta) with analytical predictions
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_figd1(q, form):
    # Setting values for influence coefficients
    lam = 1 - (1-q)*(1-0.9)
    omega = q*(1-0.9)/lam
    
    # Setting up influence coefficients
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
    
    # Getting predictions
    q = lam*omega/(1-lam*(1-omega))
    
    # Boundary predictions
    if q != 0:
        theta_L = np.arctan(2/(2-q*(1+p*a0+(1-p)*b0)))
        theta_H = np.arctan(2/(q*(1-p*a0-(1-p)*b0)))
    else:
    ## If q = 0, adjust to non-prejudiced boundaries
        theta_L = np.pi/4
        theta_H = np.pi/2
    
    m_theta = [theta for theta in theta_values if theta > theta_H]
    n_theta = [theta for theta in theta_values if theta < theta_L]
    c_theta = [theta for theta in theta_values if (theta not in m_theta) and (theta not in n_theta)]
    
    # Opinion profile predictions
    if q*(a0-b0) > r:
        a_m = q*a0 + (1-q)
        b_m = q*b0 + (1-q)
        a_n = q*a0 - (1-q)
        b_n = q*b0 - (1-q)
        a_c_values = 1 - 2/np.tan(c_theta) + (1-p)*q*(a0-b0)
        b_c_values = 1 - 2/np.tan(c_theta) - p*q*(a0-b0)    
    else:
        a_m = q*a0 + (1-q) - lam*(1-omega)*(1-p)*q*(a0-b0)
        b_m = q*b0 + (1-q) + lam*(1-omega)*p*q*(a0-b0)
        a_n = q*a0 - (1-q) - lam*(1-omega)*(1-p)*q*(a0-b0)
        b_n = q*b0 - (1-q) + lam*(1-omega)*p*q*(a0-b0)   
        a_c_values = 1 - 2/np.tan(c_theta) + (1-p)*lam*omega*(a0-b0)
        b_c_values = 1 - 2/np.tan(c_theta) - p*lam*omega*(a0-b0)    
    
    # Environmental state predictions
    alphas = m*np.sin(theta_values)*gamma
    e_m_values = m*np.sin(m_theta)*gamma/(2*gamma)*q*(1-p*a0 - (1-p)*b0)
    e_n_values = m*np.sin(n_theta)*gamma/(2*gamma)*(2-q*(1+p*a0 + (1-p)*b0))
    e_c_values = m*np.cos(c_theta)
   
    if form == 'show':
        plt.figure()
        plt.suptitle('Proportion of prejudice = ' + str((omega*lam)/(1-lam*(1-omega))))
        plt.subplot(311)
        plt.ylim(-0.1, 5.1)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('Environmental pollution')
        plt.plot(theta_values, e_result.T, ',k', alpha = 0.25)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.plot(m_theta, e_m_values, color='blue')
        plt.plot(n_theta, e_n_values, color='red')
        plt.plot(c_theta, e_c_values, color='orange')
        plt.vlines(theta_L, ymin=-0.1, ymax=5.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-0.1, ymax=5.1, color='blue',linestyle='dotted')
        
        plt.subplot(312)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('a - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.hlines(a_m, xmin=theta_H, xmax=np.pi, color='blue')
        plt.hlines(a_n, xmin=0, xmax=theta_L, color='red')
        plt.plot(c_theta, a_c_values, color='orange')
        plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
        
        plt.subplot(313)
        plt.xlabel('Relative cost of mitigation')
        plt.ylabel('b - opinion')
        plt.ylim(-1.1, 1.1)
        plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.hlines(b_m, xmin=theta_H, xmax=np.pi, color='blue')
        plt.hlines(b_n, xmin=0, xmax=theta_L, color='red')
        plt.plot(c_theta, b_c_values, color='orange')
        plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
        plt.show()
    
    if form == 'export':
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-0.1, 5.1)
        plt.plot(theta_values, e_result.T, ',k', alpha = 0.15)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.plot(m_theta, e_m_values, color='blue')
        plt.plot(n_theta, e_n_values, color='red')
        plt.plot(c_theta, e_c_values, color='orange')
        plt.vlines(theta_L, ymin=-0.1, ymax=5.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-0.1, ymax=5.1, color='blue',linestyle='dotted')
        plt.savefig("FigureD1_e.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.hlines(a_m, xmin=theta_H, xmax=np.pi, color='blue')
        plt.hlines(a_n, xmin=0, xmax=theta_L, color='red')
        plt.plot(c_theta, a_c_values, color='orange')
        plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
        plt.savefig("FigureD1_a.png") 
        
        plt.figure(figsize=(8, 3.15))
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.ylim(-1.1, 1.1)
        plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
        plt.xticks(np.arange(0, np.pi+0.1, step=np.pi/4)) 
        ## Plotting predictions
        plt.hlines(b_m, xmin=theta_H, xmax=np.pi, color='blue')
        plt.hlines(b_n, xmin=0, xmax=theta_L, color='red')
        plt.plot(c_theta, b_c_values, color='orange')
        plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
        plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
        plt.savefig("FigureD1_b.png")         
        
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