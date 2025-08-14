############################################
# Bifurcation diagrams over beta and alpha #
############################################

# Contains code to generate Figure 5

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Fix all parameters except, alpha and beta
tau = 1.5
p = 0.1
r = 0.2
# a-(0.9, 0), b-(0.95, 0.05/0.95), c-(0.99, 1/11)
lam = 0.95
gamma = 0.2
omega = 0.05/0.95
a0 = 1
b0 = -1

## Bifurcation over beta ##
alpha = 1 
 
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

# Setting up bifurcation diagram loop
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
    a, b, e = HK_prejudice_map(a, b, e, beta_values)
    if i >= (iterations - last):
        a_result[i - (iterations - last)] = a
        b_result[i - (iterations - last)] = b
        e_result[i - (iterations - last)] = e
    

# Plot
plt.figure()
plt.subplot(211)
plt.xlabel('Net cost of mitigation')
plt.ylabel('Environmental pollution')
plt.plot(beta_values, e_result.T, ',k', alpha = 0.25)


plt.subplot(223)
plt.xlabel('Net cost of mitigation')
plt.ylabel('b - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(beta_values, b_result.T, ',k', alpha = 0.25)

plt.subplot(224)
plt.xlabel('Net cost of mitigation')
plt.ylabel('a - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(beta_values, a_result.T, ',k', alpha = 0.25)

plt.show()

# For figure
plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-0.1, 5.1)
plt.plot(beta_values, e_result.T, ',k', alpha = 0.15)
plt.savefig("Figure4a_e.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(beta_values, a_result.T, ',k', alpha = 0.25)
plt.savefig("Figure4a_a.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(beta_values, b_result.T, ',k', alpha = 0.25)
plt.savefig("Figure4a_b.pdf") 

## Bifurcation over alpha ##
beta = 4

# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e, alpha):
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

# Setting up bifurcation diagram loop
alpha_values = np.linspace(0, 5, 10000)
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
    a, b, e = HK_prejudice_map(a, b, e, alpha_values)
    if i >= (iterations - last):
        a_result[i - (iterations - last)] = a
        b_result[i - (iterations - last)] = b
        e_result[i - (iterations - last)] = e
    

# Plot
plt.figure()
plt.subplot(211)
plt.xlabel('Rate of emissions')
plt.ylabel('Environmental pollution')
plt.plot(alpha_values, e_result.T, ',k', alpha = 0.25)


plt.subplot(223)
plt.xlabel('Rate of emissions')
plt.ylabel('b - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(alpha_values, b_result.T, ',k', alpha = 0.25)

plt.subplot(224)
plt.xlabel('Rate of emissions')
plt.ylabel('a - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(alpha_values, a_result.T, ',k', alpha = 0.25)

plt.show()

# For figure
plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-0.1, 5.1)
plt.plot(alpha_values, e_result.T, ',k', alpha = 0.15)
plt.savefig("Figure5b_e.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(alpha_values, a_result.T, ',k', alpha = 0.25)
plt.savefig("Figure5b_a.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(alpha_values, b_result.T, ',k', alpha = 0.25)
plt.savefig("Figure5b_b.pdf") 