####################################
# Figure 2 - Polarization heatmaps #
####################################

## INSTRUCTIONS ##
# Run this script to load in the required functions 'gen_fig2a', 'gen_fig2b', 'gen_fig2c'
# Running each of these functions with the appropriate parameters will regenerate the associated heatmap from Fig2

# Load required libraries
import numpy as np
import matplotlib.pyplot as plt

# Fixed parameters used in figures
tau = 1.5
p = 0.1
a0 = 1
b0 = -1

## gen_fig2a ##
# Returns the polarization heatmap over prejudice and confidence bound
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig2a(form):
    
    # Setting up sampling grid and averaging window
    pred_values = np.linspace(0,1,300)
    social_values = (1-pred_values)/2
    r_values = np.linspace(0,2,300)
    iterations = 1000
    last = 100
    
    # Initialize
    P = np.zeros(len(pred_values))
    P[:] = a0-b0
    P_result = np.empty((last, len(pred_values)))
    
    pred_averages = []
    r_averages = []
    P_averages = []
    
    # Iterate
    for r in r_values:
        for i in range(iterations):
            P = polarization_map(P, pred_values, social_values, r)
            if i >= (iterations - last):
                P_result[i - (iterations - last)] = P
        for k in range(len(pred_values)):
            P_averages.append(np.average(P_result[:, k]))
            r_averages.append(r)
            pred_averages.append(pred_values[k])    
    
    if form == 'show':
        plt.figure()
        plt.xlabel('Prejudice')
        plt.ylabel('Confidence bound')
        plt.scatter(pred_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar(label='Polarization')        
        plt.show()
    
    if form == 'export':
        plt.figure(figsize=(8, 5))
        plt.scatter(pred_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.savefig("Figure2a.png")         

## gen_fig2b ##
# Returns the polarization heatmap over objectivity and confidence bound
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig2b(form):
    
    # Setting up sampling grid and averaging window
    obj_values = np.linspace(0,1,300)
    pred_values = (1-obj_values)/2
    social_values = pred_values
    r_values = np.linspace(0,2,300)
    iterations = 1000
    last = 100
    
    # Initialize
    P = np.zeros(len(obj_values))
    P[:] = a0-b0
    P_result = np.empty((last, len(pred_values)))
    
    obj_averages = []
    r_averages = []
    P_averages = []
    
    # Iterate
    for r in r_values:
        for i in range(iterations):
            P = polarization_map(P, pred_values, social_values, r)
            if i >= (iterations - last):
                P_result[i - (iterations - last)] = P
        for k in range(len(obj_values)):
            P_averages.append(np.average(P_result[:, k]))
            r_averages.append(r)
            obj_averages.append(obj_values[k])
    
    if form == 'show':
        plt.figure()
        plt.xlabel('Objectivity')
        plt.ylabel('Confidence bound')
        plt.scatter(obj_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar(label='Polarization')
        plt.show()
    
    if form == 'export':
        plt.figure(figsize=(8, 5))
        plt.scatter(obj_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.savefig("Figure2b.png")    

## gen_fig2c ##
# Returns the polarization heatmap over the prejudice-to-objectivity ratio and confidence bound
# To show the plots put form = 'show', to export the plots of the exact size used in the figures put form = 'export'

def gen_fig2c(form):
    
    # Setting up sampling grid and averaging window
    q_values = np.linspace(0,1,300)
    social_values = [0.25]*300
    pred_values = q_values*(1-0.25)
    obj_values = (1-q_values)*(1-0.25)
    r_values = np.linspace(0,2,300)
    iterations = 1000
    last = 100
    
    # Initialize
    P = np.zeros(len(obj_values))
    P[:] = a0-b0
    P_result = np.empty((last, len(pred_values)))
    
    # Initialize
    P = np.zeros(len(obj_values))
    P[:] = a0-b0
    P_result = np.empty((last, len(pred_values)))
    
    q_averages = []
    r_averages = []
    P_averages = []
    
    # Iterate
    for r in r_values:
        for i in range(iterations):
            P = polarization_map(P, pred_values, social_values, r)
            if i >= (iterations - last):
                P_result[i - (iterations - last)] = P
        for k in range(len(obj_values)):
            P_averages.append(np.average(P_result[:, k]))
            r_averages.append(r)
            q_averages.append(q_values[k])
    
    if form == 'show':
        plt.figure()
        plt.xlabel('Prejudice-to-Objectivity ratio (q)')
        plt.ylabel('Confidence bound')
        plt.scatter(q_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar(label='Polarization')
        plt.show()        
    
    if form == 'export':
        plt.figure(figsize=(8, 5))
        plt.scatter(q_averages, r_averages, c=P_averages, vmin=0, vmax=2)
        cbar = plt.colorbar()
        plt.setp(cbar.ax.get_yticklabels(), color='white')
        plt.xticks(color='w')
        plt.yticks(color='w')
        plt.savefig("Figure2c.png")        

## polarization_map ##
# Returns the polarization at the next time step

def polarization_map(P, pred, social, r):
    P_next = np.empty(len(P))
    for i in range(len(P)):
        if P[i] > r:
            P_next[i] = (tau - 1 + social[i])/tau*P[i] + pred[i]/tau*(a0-b0)
        else:
            P_next[i] = (tau - 1)/tau*P[i] + pred[i]/tau*(a0-b0)
    return P_next