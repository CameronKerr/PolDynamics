#########################
# Numerical simulations #
#########################

# Code to generate Figure 1b

# Load required libraries
import numpy as np
import matplotlib.pyplot as plt

# Fix all parameters
tau = 1.5
p = 0.5
r = 0.2
lam = 2/3
gamma = 0.2
omega = 1/2
a0 = 0.5
b0 = -0.5
alpha = 1.5

# Cycle example
#beta = 4
# Mitigative example
#beta = 1

# Nonmitigative convergence
omega = 0
beta = 7.6


# Map returns the next iteration in the trajectory 
def HK_prejudice_map(a, b, e):
    e_next = (1 - gamma)*e + 1/2*alpha*(1 - (p*a + (1-p)*b))
    
    if abs(a-b) > r:
        a_next = (1 - 1/tau)*a + 1/tau*(lam*((1-omega)*a + omega*a0) + (1-lam)*np.sign(e - beta))
        b_next = (1 - 1/tau)*b + 1/tau*(lam*((1-omega)*b + omega*b0) + (1-lam)*np.sign(e - beta))
    else:
        a_next = (1 - 1/tau)*a + 1/tau*(lam*((1-omega)*(p*a + (1-p)*b) + omega*a0) + (1-lam)*np.sign(e - beta))
        b_next = (1 - 1/tau)*b + 1/tau*(lam*((1-omega)*(p*a + (1-p)*b) + omega*b0) + (1-lam)*np.sign(e - beta))    
            
    return a_next, b_next, e_next

# Get trajectory
iterations = 100
a = [a0]
b = [b0]
e = [0]
# Iterate
for i in range(iterations):
    a_next, b_next, e_next = HK_prejudice_map(a[-1], b[-1], e[-1])
    a.append(a_next)
    b.append(b_next)
    e.append(e_next)

# Plot
plt.figure()
plt.suptitle("")
plt.subplot(211)
plt.ylim(0, alpha/gamma+1)
plt.xlim(0, 101)
plt.ylabel('Environmental pollution - e(t)')
plt.plot(list(range(iterations + 1)), e, 'green', alpha = 0.5)

plt.subplot(212)
plt.xlabel('Time')
plt.ylabel('Opinion - a(t), b(t)')
plt.ylim(-1.1, 1.1)
plt.xlim(0, 101)
plt.plot(list(range(iterations + 1)), a, 'blue', alpha=0.5)
plt.plot(list(range(iterations + 1)), b, 'red', alpha=0.5)
plt.show()

# For figure
plt.figure(figsize=(4, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-0.1, alpha/gamma + 1)
plt.xlim(0, 50)
plt.plot(list(range(iterations + 1)), e, 'green', alpha=0.5)
plt.savefig('Figure2d_e.pdf')

plt.figure(figsize=(4, 3.15))
plt.xticks(color='w')
plt.yticks(color='w')
plt.ylim(-1.1, 1.1)
plt.xlim(0, 50)
plt.plot(list(range(iterations + 1)), a, 'blue', alpha=0.5)
plt.plot(list(range(iterations + 1)), b, 'red', alpha=0.5)
plt.savefig('Figure2d_o.pdf')