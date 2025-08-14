##################################
# Bifurcation diagram over theta #
##################################

# Contains code to generate Figure 2

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Fix all parameters except, alpha and beta
tau = 1.5
p = 0.1
r = 0.2
m = 5 # Magnitude in (beta, alpha) parameter space
# a-(0.9, 0), b-(0.95, 0.05/0.95), c-(0.99, 1/11)
lam = 0.99
gamma = 0.2
omega = 1/11
a0 = 1
b0 = -1

# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e, alpha, beta):
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

# Setting up bifurcation diagram loop
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
    a, b, e = HK_prejudice_map(a, b, e, alpha_values, beta_values)

    if i >= (iterations - last):
        a_result[i - (iterations - last)] = a
        b_result[i - (iterations - last)] = b
        e_result[i - (iterations - last)] = e

# Plot
plt.figure()
plt.suptitle('Proportion of prejudice = ' + str((omega*lam)/(1-lam*(1-omega))))
plt.subplot(211)
plt.ylim(-0.1, 5.1)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('Environmental pollution')
plt.plot(theta_values, e_result.T, ',k', alpha = 0.25)

plt.subplot(223)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('a - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)

plt.subplot(224)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('b - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)

plt.show()

# For figure
plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-0.1, 5.1)
plt.plot(theta_values, e_result.T, ',k', alpha = 0.15)
plt.savefig("Figure2c_e.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
plt.savefig("Figure2c_a.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
plt.savefig("Figure2c_b.pdf") 






