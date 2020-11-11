clc, clear all, close all

N = 1001;

wo = 0.01;

x1   = 0;
x2   = 0;
yLTI1= 0;
yLPV1= 0;
u1   = 0;
u2   = 0;

x     = zeros(N, 1);
yLPV  = zeros(N, 1);
yLTI  = zeros(N, 1);
y2LPV = zeros(N, 1);
u     = zeros(N, 1);
w     = ones(N, 1);
w1    = 0;
theta = zeros(N, 1);
t     = linspace(0, 1000, N);

for k=1:N
    theta(k) = 0.5 * cos(3*wo*t(k));
    u(k)     = sin(wo*t(k));
    a0       = -0.78 + 0.44*theta(k);
    d0       =  0.30 + 0.90*theta(k);

    % LPV synthesis
    x(k)     = -a0*x1 + u1;
    yLPV(k)  =  d0*x1;

    % LTI synthesis
    yLTI(k) = -a0*yLTI1 + d0*u2;

    % LPV synthesis
    x(k)  = -a0*x1 + u1;
    yLPV(k)  = d0*x1;
    y2LPV(k) = d0/-a0*(x(k) - u1);
    
    A_  = [0  1];
    B_  = [d0 0];
    Ak_ = [a0 1];
    Bk_ = [1  0];
    
    R = [A_ -B_
         Bk_ Ak_];
    H = [zeros(size(A_,1),size(Bk_,2))
         Bk_];
    
    w_ = [w1
          w(k)];
    
    eta = pinv(R'*R)*R'*H*w;
    
    x2 = x1;
    x1 = x(k);
    yLPV1 = yLPV(k);
    yLTI1 = yLTI(k);

    u2 = u1;
    u1 = u(k);
    
    w1 = w(k);
end
% ax = fig.add_subplot(2, 1, 1)
plot(t, yLPV, 'LineWidth', 2)
hold on
plot(t, yLTI, 'LineWidth', 2)
grid on
% ax.plot(t, y2LPV)
% ay = fig.add_subplot(2, 1, 2)
% ay.plot(t, theta)
% plot.show()
% print(A_)
% print(B_)
% print(Ak_)
% print(Bk_)
% print(R)
% print(H)