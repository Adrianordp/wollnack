#!/usr/bin/python3.8
import numpy as np
import matplotlib.pyplot as plot

N = 1001

wo = 0.01

x1   = 0
x2   = 0
yLTI1= 0
u1   = 0
u2   = 0

x     = np.zeros(N)
yLPV  = np.zeros(N)
yLTI  = np.zeros(N)
y2LPV = np.zeros(N)
u     = np.zeros(N)
theta = np.zeros(N)
t     = np.linspace(0, 1000, N)
for k in range(N):
    theta[k] = 0.5 * np.cos(3*wo*t[k])
    u[k]     = np.sin(wo*t[k])
    a0       = -0.78 + 0.44*theta[k]
    d0       = 0.3 + 0.9*theta[k]

    #A = a0 + q
    #B = 1
    #C = q
    #D = d0

    # LPV synthesis
    x[k]  = (0.78 - 0.44*theta[k])*x1 + u1
    yLPV[k]  = (0.3 + 0.9*theta[k])*x1
    y2LPV[k] = (0.3 + 0.9*theta[k])/(0.78 - 0.44*theta[k])*(x[k] - u1)

    #LTI synthesis
    yLTI[k] = (0.78 - 0.44*theta[k])*yLTI1 + (0.3 + 0.9*theta[k])*u2

    x2 = x1
    x1 = x[k]
    yLTI1 = yLTI[k]

    u2 = u1
    u1 = u[k]
    
fig = plot.figure()
ax = fig.add_subplot(2, 1, 1)
ax.plot(t, yLPV)
ax.plot(t, yLTI)
# ax.plot(t, y2LPV)
ay = fig.add_subplot(2, 1, 2)
ay.plot(t, theta)
plot.show()