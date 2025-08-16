##################################
# Influence coefficients heatmap #
##################################

# Contains code to generate Figure C.2

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt
import mpltern

# Fixed parameters
tau = 1.5
p = 0.5
r = 0.2
alpha = 0.8
beta = 5
gamma = 0.2
a0 = 1
b0 = -1

# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e, omega):
    e_next = (1 - gamma)*e + 1/2*alpha*(1 - (p*a + (1-p)*b))
    # Initialize empty iteration
    a_next = np.empty(len(a))
    b_next = np.empty(len(b))
    
    for i in range(len(a)):
        if abs(a[i] - b[i]) > r:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega[i])*a[i] + omega[i]*a0) + (1-lam)*np.sign(e[i] - beta))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega[i])*b[i] + omega[i]*b0) + (1-lam)*np.sign(e[i] - beta))
        else:
            a_next[i] = (1 - 1/tau)*a[i] + 1/tau*(lam*((1-omega[i])*(p*a[i] + (1-p)*b[i]) + omega[i]*a0) + (1-lam)*np.sign(e[i] - beta))
            b_next[i] = (1 - 1/tau)*b[i] + 1/tau*(lam*((1-omega[i])*(p*a[i] + (1-p)*b[i]) + omega[i]*b0) + (1-lam)*np.sign(e[i] - beta))
            
    return a_next, b_next, e_next

# Setting up bifurcation diagram loop
omega_values = np.linspace(0, 1, 100)
lam_values = np.linspace(0, 1, 100)
iterations = 1000
last = 100
# Initialize
a = np.zeros(len(omega_values))
b = np.zeros(len(omega_values))
e = np.zeros(len(omega_values))
a[:] = a0
b[:] = b0
e[:] = 0
a_result = np.empty((last, len(omega_values)))
b_result = np.empty((last, len(omega_values)))
e_result = np.empty((last, len(omega_values)))

social_averages = []
prejudice_averages = []
truth_averages = []
c_averages = []
e_averages = []

# Iterate
for lam in lam_values:
    for i in range(iterations):
        a, b, e = HK_prejudice_map(a, b, e, omega_values)
        if i >= (iterations - last):
            a_result[i - (iterations - last)] = a
            b_result[i - (iterations - last)] = b
            e_result[i - (iterations - last)] = e
    for k in range(len(omega_values)):
        e_averages.append(np.average(e_result[:, k]))
        c_averages.append(np.average(np.concatenate((a_result[:, k], b_result[:, k]))))
        social_averages.append(lam*(1-omega_values[k]))
        prejudice_averages.append(lam*omega_values[k])
        truth_averages.append(1 - lam)
    

fig = plt.figure()
fig.suptitle("Cost of mitigation: " + str(beta) + " Worst case pollution: " + str(alpha/gamma))
fig.subplots_adjust(wspace=0.3)
ax = fig.add_subplot(1, 2, 1, projection="ternary")
pc = ax.scatter(social_averages, prejudice_averages, truth_averages, c=e_averages, vmin=0, vmax=4)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9])
colorbar = fig.colorbar(pc, cax=cax)
colorbar.set_label("Pollutant level")
ax.set_tlabel("Social susceptability")
ax.set_llabel("Prejudice")
ax.set_rlabel("Truth-seeking")

ax = fig.add_subplot(1, 2, 2, projection="ternary")
pc = ax.scatter(social_averages, prejudice_averages, truth_averages, c=c_averages, vmin=-1, vmax=1)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9])
colorbar = fig.colorbar(pc, cax=cax)
colorbar.set_label("Level of mitigation")
ax.set_tlabel("Social susceptability")
ax.set_llabel("Prejudice")
ax.set_rlabel("Objectivity")
plt.show()

## Figure
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(1, 1, 1, projection="ternary")
pc = ax.scatter(social_averages, prejudice_averages, truth_averages, c=e_averages, vmin=0, vmax=4)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.9])
ax.tick_params(labelbottom=False, labelleft=False)
colorbar = fig.colorbar(pc, cax=cax)
colorbar.ax.tick_params(labelright=False)
plt.savefig("Figure4.pdf") 
