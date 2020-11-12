clc, clear all, close all

N = 1001;

wo = 0.01;

x1   = 0;
x2   = 0;
yLTI1= 0;
yLPV1= 0;
u1   = 0;
u2   = 0;
uMF1 = 0;
yMF1 = 0;

x     = zeros(N, 1);
yLPV  = zeros(N, 1);
yLTI  = zeros(N, 1);
y2LPV = zeros(N, 1);
yMF   = zeros(N, 1);
uMF   = zeros(N, 1);
u     = zeros(N, 1);
w     = ones(N, 1);
w1    = 0;
theta = zeros(N, 1);
t     = linspace(0, 1000, N);

for k=1:N
    theta(k) =  0.5 * cos(3*wo*t(k));
    u(k)     =  sin(wo*t(k));
    a0       = -0.78 + 0.44*theta(k);
    d0       =  0.30 + 0.90*theta(k);

    % LPV synthesis
    x(k)     = -a0*x1 + u1;
    yLPV(k)  =  d0*x1;

    % LTI synthesis
    yLTI(k) = -a0*yLTI1 + d0*u2;

    % LPV synthesis
    x(k)     = -a0*x1 + u1;
    yLPV(k)  =  d0*x1;
    y2LPV(k) =  d0/-a0*(x(k) - u1);
    
    A_  = [0  1];
    B_  = [d0 0];
%     Ak_ = [a0 1];
    Ak_ = [-1 1];
    Bk_ = [d0  0];
%     Bk_ = B_;
    
    w_ = [w1
          w(k)];
    
%     R = [A_ -B_
%          Bk_ Ak_];
%     H = [zeros(size(A_, 1),size(Bk_, 2))
%          Bk_];
    
    R_l = [ A_(end) -B_(end)
           Bk_(end) Ak_(end)];
    H_l = [-A_(1:end-1)  B_(1:end-1) zeros(size(Bk_))
          -Bk_(1:end-1) Ak_(1:end-1) Bk_];
    eta_l = [yMF1
             uMF1
             w_];
    eta = R_l\H_l*eta_l;
    
    yMF(k) = eta(1);
    uMF(k) = eta(2);
    
    x2 = x1;
    x1 = x(k);
    yLPV1 = yLPV(k);
    yLTI1 = yLTI(k);

    u2 = u1;
    u1 = u(k);
    
    uMF1 = uMF(k);
    yMF1 = yMF(k);
    
    w1 = w(k);
end

subplot(2,1,1), plot(t, yLPV, 'LineWidth', 2)
hold on
plot(t, yLTI, 'LineWidth', 2)
hold off
grid on
legend('yLPV','yLTI')

subplot(2,1,2), plot(t, yMF, 'LineWidth', 2)
hold on
plot(t, uMF, 'LineWidth', 2)
plot(t, w, 'LineWidth', 2)
hold off
grid on
legend('yMF','uMF','w')