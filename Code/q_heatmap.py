###################
# theta-q heatmap #
###################

# Contains code to generate Figure 4d

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt

# Map returns the next iteration in the trajectory 
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

# Fixed all parameters except theta and q

tau = 1.5
p= 0.1
r= 0.2
social = 0.25
gamma = 0.2
m = 4
a0 = 1
b0 = -1
# Setting up bifurcation diagram loop
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
    lam = 1 - (1-q)*(1-social)
    omega = q*(1-social)/lam
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


plt.figure()
plt.suptitle('')
plt.subplot(211)
plt.ylabel('')
plt.scatter(q_averages, theta_averages, c=e_averages, vmin=0, vmax=4) 
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

plt.figure(figsize=(6.15, 6.15))
plt.subplot(211)
plt.xticks(color='w')
plt.yticks(color='w')
plt.yticks(np.arange(0, (np.pi + 0.1), step=(np.pi/4)))
plt.scatter(q_averages, theta_averages, c=e_averages, vmin=0, vmax=4) 
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.subplot(212)
plt.xticks(color='w')
plt.yticks(color='w')
plt.yticks(np.arange(0, (np.pi + 0.1), step=(np.pi/4)))
plt.scatter(q_averages, theta_averages, c=c_averages, vmin=-1, vmax=1) 
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.savefig("figure3d.pdf") 