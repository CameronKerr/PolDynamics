##########################
# Cost-emissions heatmap #
##########################

# Contains code to generate Figure 1

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e, beta):
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

# Fixed all parameters except beta and alpha
tau = 1.5
p = 0.1
r = 0.2
lam = 0.9
gamma = 0.2
omega = 0
a0 = 1
b0 = -1
# a-(0.9, 0), b-(0.95, 0.05/0.95), c-(0.99, 1/11)

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
        a, b, e = HK_prejudice_map(a, b, e, beta_values)
        if i >= (iterations - last):
            a_result[i - (iterations - last)] = a
            b_result[i - (iterations - last)] = b
            e_result[i - (iterations - last)] = e
    for k in range(len(beta_values)):
        e_averages.append(np.average(e_result[:, k]))
        c_averages.append(np.average(np.concatenate((a_result[:, k], b_result[:, k]))))
        alpha_averages.append(alpha)
        beta_averages.append(beta_values[k])

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

# For figure
plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.scatter(beta_averages, alpha_averages, c=e_averages, vmin=0, vmax=10)
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.savefig("Figure1c_e.pdf") 

plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.scatter(beta_averages, alpha_averages, c=c_averages, vmin=-1, vmax=1) 
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.savefig("Figure1c_o.pdf") 