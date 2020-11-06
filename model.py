#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plot

N = 1001

wo = 0.01

x1   = 0
yLTI1= 0
u1   = 0
u2   = 0

x     = np.zeros(N)
yLPV  = np.zeros(N)
yLTI  = np.zeros(N)
u     = np.zeros(N)
theta = np.zeros(N)
t     = np.linspace(0, 1000, N)
for k in range(N):
    theta[k] = 0.5 * np.cos(3*wo*t[k])
    u[k]     = np.sin(wo*t[k])
    a0       = -0.78 + 0.44*theta[k]
    d0       = 0.3 + 0.9*theta[k]

    # LPV synthesis
    x[k]     = -a0*x1 + u1
    yLPV[k]  = d0*x1

    #LTI synthesis
    yLTI[k] = -a0*yLTI1 + d0*u2

    #A = a0 + q
    #B = 1
    #C = q
    #D = d0
    A_  = np.array([a0, 1])
    B_  = np.array([1, 0])
    Ak_ = np.array([0, 1])
    Bk_ = np.array([0, 1])
    
    R = np.matrix([np.concatenate((A_, -B_)),np.concatenate((Bk_, Ak_))])
    H = np.matrix([np.zeros(2), Bk_])

    y_ = np.array(yLPV1, yLPV[k])
    u_ = np.array(u1, u[k])


    x1 = x[k]
    yLPV1 = yLPV[k]
    yLTI1 = yLTI[k]

    u2 = u1
    u1 = u[k]
    
fig = plot.figure()
ax = fig.add_subplot(2, 1, 1)
ax.plot(t, yLPV)
ax.plot(t, yLTI)
ay = fig.add_subplot(2, 1, 2)
ay.plot(t, theta)
plot.show()
print(A_)
print(B_)
print(Ak_)
print(Bk_)
print(R)
print(H)