#!/usr/bin/python3.8
import numpy as np
import matplotlib.pyplot as plot

N = 101

wo = 0.01

x1 = 0

x = np.zeros(N)
y = np.zeros(N)
t = np.linspace(0, 10, N)
for k in range(N):
    theta = 0.5 * np.cos(3*wo*t[k])
    u     = np.sin(wo*t[k])
    a0    = -0.78 + 0.44*theta
    d0    = 0.3 + 0.9*theta

    #A = a0 + q
    #B = 1
    #C = q
    #D = d0

    x[k] = 0.78*x1 - 0.44*theta*x1 + u
    y[k] = 0.3*x1 + 0.9*theta*x1

    x1 = x[k]
print(x)
print(y)
fig = plot.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(t,x)
ax.plot(t,y)
plot.show()