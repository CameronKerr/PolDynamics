###################################################
# Bifurcation diagram over theta with predictions #
###################################################

# Contains code to generate Figure 3

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
lam = 0.9
gamma = 0.2
omega = 0
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
    alpha_values = m*np.sin(theta_values)*gamma
    beta_values = m*np.cos(theta_values)
    a, b, e = HK_prejudice_map(a, b, e, alpha_values, beta_values)

    if i >= (iterations - last):
        a_result[i - (iterations - last)] = a
        b_result[i - (iterations - last)] = b
        e_result[i - (iterations - last)] = e

# Getting predictions
q = lam*omega/(1-lam*(1-omega))

# Boundary predictions
theta_L = np.arctan(2/(2-q*(1+p*a0+(1-p)*b0)))
theta_H = np.arctan(2/(q*(1-p*a0-(1-p)*b0)))
## If q = 0, adjust to non-prejudiced boundaries
#theta_L = np.pi/4
#theta_H = np.pi/2

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

# Plot
plt.figure()
plt.suptitle('Proportion of prejudice = ' + str((omega*lam)/(1-lam*(1-omega))))
plt.subplot(211)
plt.ylim(-0.1, 5.1)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('Environmental pollution')
plt.plot(theta_values, e_result.T, ',k', alpha = 0.25)
## Plotting predictions
plt.plot(m_theta, e_m_values, color='blue')
plt.plot(n_theta, e_n_values, color='red')
plt.plot(c_theta, e_c_values, color='orange')
plt.vlines(theta_L, ymin=-0.1, ymax=5.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-0.1, ymax=5.1, color='blue',linestyle='dotted')

plt.subplot(223)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('a - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
## Plotting predictions
plt.hlines(a_m, xmin=theta_H, xmax=np.pi, color='blue')
plt.hlines(a_n, xmin=0, xmax=theta_L, color='red')
plt.plot(c_theta, a_c_values, color='orange')
plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')

plt.subplot(224)
plt.xlabel('Relative cost of mitigation')
plt.ylabel('b - opinion')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
## Plotting predictions
plt.hlines(b_m, xmin=theta_H, xmax=np.pi, color='blue')
plt.hlines(b_n, xmin=0, xmax=theta_L, color='red')
plt.plot(c_theta, b_c_values, color='orange')
plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')

plt.show()

# For figure
plt.figure(figsize=(8, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-0.1, 5.1)
plt.plot(theta_values, e_result.T, ',k', alpha = 0.15)
## Plotting predictions
plt.plot(m_theta, e_m_values, color='blue')
plt.plot(n_theta, e_n_values, color='red')
plt.plot(c_theta, e_c_values, color='orange')
plt.vlines(theta_L, ymin=-0.1, ymax=5.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-0.1, ymax=5.1, color='blue',linestyle='dotted')
plt.savefig("Figure3c_e.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, a_result.T, ',k', alpha = 0.25)
## Plotting predictions
plt.hlines(a_m, xmin=theta_H, xmax=np.pi, color='blue')
plt.hlines(a_n, xmin=0, xmax=theta_L, color='red')
plt.plot(c_theta, a_c_values, color='orange')
plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
plt.savefig("Figure3c_a.pdf") 

plt.figure(figsize=(3.15, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.plot(theta_values, b_result.T, ',k', alpha = 0.25)
## Plotting predictions
plt.hlines(b_m, xmin=theta_H, xmax=np.pi, color='blue')
plt.hlines(b_n, xmin=0, xmax=theta_L, color='red')
plt.plot(c_theta, b_c_values, color='orange')
plt.vlines(theta_L, ymin=-1.1, ymax=1.1, color='red', linestyle='dotted')
plt.vlines(theta_H, ymin=-1.1, ymax=1.1, color='blue',linestyle='dotted')
plt.savefig("Figure3c_b.pdf") 






