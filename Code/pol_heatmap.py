####################################################
# Prejudice-confidence bound polarization heat map #
####################################################

# Code to generate Figure 2

# Load required libraries:
import numpy as np
import matplotlib.pyplot as plt

# Fixed all parameters except, beta
tau = 1.5
p = 0.1
social = 0.5
a0 = 1
b0 = -1

# Map returns the polarization in the next iteration in the trajectory
def polarization_map(P, pred, social):
    P_next = np.empty(len(P))
    for i in range(len(P)):
        if P[i] > r:
            P_next[i] = (tau - 1 + social[i])/tau*P[i] + pred[i]/tau*(a0-b0)
        else:
            P_next[i] = (tau - 1)/tau*P[i] + pred[i]/tau*(a0-b0)
    return P_next

## Changing prejudice ##

# Setting up bifurcation diagram loop
pred_values = np.linspace(0, 1, 300)
social_values = (1-pred_values)/2
r_values = np.linspace(0, 2, 300)
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
        P = polarization_map(P, pred_values, social_values)
        if i >= (iterations - last):
            P_result[i - (iterations - last)] = P
    for k in range(len(pred_values)):
        P_averages.append(np.average(P_result[:, k]))
        r_averages.append(r)
        pred_averages.append(pred_values[k])


plt.figure()
plt.subplot(211)
plt.xlabel('Prejudice')
plt.ylabel('Confidence bound')
plt.plot(pred_values, [pred_values[i]*(a0-b0)/(1-social_values[i]) for i in range(len(pred_values))])
plt.scatter(pred_averages, r_averages, c=P_averages, vmin=0, vmax=2)
plt.show()

# For figure
plt.figure(figsize=(8, 5))
plt.scatter(pred_averages, r_averages, c=P_averages, vmin=0, vmax=2)
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.xticks(color='w')
plt.yticks(color='w')
plt.savefig("PolarizationHeatmap_prej.png") 

## Changing objectivity

# Setting up bifurcation diagram loop
obj_values = np.linspace(0, 1, 300)
pred_values = (1-obj_values)/2
social_values = pred_values
r_values = np.linspace(0, 2, 300)
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
        P = polarization_map(P, pred_values, social_values)
        if i >= (iterations - last):
            P_result[i - (iterations - last)] = P
    for k in range(len(obj_values)):
        P_averages.append(np.average(P_result[:, k]))
        r_averages.append(r)
        obj_averages.append(obj_values[k])


plt.figure()
plt.subplot(211)
plt.xlabel('Objectivity')
plt.ylabel('Confidence bound')
plt.scatter(obj_averages, r_averages, c=P_averages, vmin=0, vmax=2)
plt.show()

# For figure
plt.figure(figsize=(8, 5))
plt.scatter(obj_averages, r_averages, c=P_averages, vmin=0, vmax=2)
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.xticks(color='w')
plt.yticks(color='w')
plt.savefig("PolarizationHeatmap_obj.png") 

## Changing q

# Setting up bifurcation diagram loop
q_values = np.linspace(0,1, 300)
social_values = [0.25]*300
pred_values = q_values*(1-0.25)
obj_values = (1-q_values)*(1-0.25)
r_values = np.linspace(0, 2, 300)
iterations = 1000
last = 100
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
        P = polarization_map(P, pred_values, social_values)
        if i >= (iterations - last):
            P_result[i - (iterations - last)] = P
    for k in range(len(obj_values)):
        P_averages.append(np.average(P_result[:, k]))
        r_averages.append(r)
        q_averages.append(q_values[k])


plt.figure()
plt.subplot(211)
plt.xlabel('Prejudice-to-Objectivity ratio (q)')
plt.ylabel('Confidence bound')
plt.scatter(q_averages, r_averages, c=P_averages, vmin=0, vmax=2)
plt.show()

# For figure
plt.figure(figsize=(8, 5))
plt.scatter(q_averages, r_averages, c=P_averages, vmin=0, vmax=2)
cbar = plt.colorbar()
plt.setp(cbar.ax.get_yticklabels(), color='white')
plt.xticks(color='w')
plt.yticks(color='w')
plt.savefig("PolarizationHeatmap_q.png")