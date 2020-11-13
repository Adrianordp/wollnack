clc, clear all, close all

Ts = 60;

A = [1 theta_1(100*Ts, Ts) theta_2(100*Ts, Ts)
     1 theta_1(300*Ts, Ts) theta_2(300*Ts, Ts)
     1 theta_1(900*Ts, Ts) theta_2(900*Ts, Ts)];
Bk0 = [4.7456
       3.3421
       2.3070];
Bk1 = [-3.9231
       -2.3846
       -0.1758];
   
% A = [1 theta_1(100) theta_2(100)
%      1 theta_1(700) theta_2(700)
%      1 theta_1(900) theta_2(900)];
% Bk0 = [4.7456
%        2.4737
%        2.3070];
% Bk1 = [-3.9231
%        -0.3077
%        -0.0549];

xBk0 = A\Bk0;
bk00 = xBk0(1);
bk01 = xBk0(2);
bk02 = xBk0(3);

xBk1 = A\Bk1;
bk10 = xBk1(1);
bk11 = xBk1(2);
bk12 = xBk1(3);

disp(['bk0(θ(t)) = ' num2str(bk00) ' +' num2str(bk01) '*θ_1(t) +' num2str(bk02) '*θ_2(t)'])
disp(['bk1(θ(t)) = ' num2str(bk10) ' +' num2str(bk11) '*θ_1(t) ' num2str(bk12) '*θ_2(t)'])

str = '#bbbbbb';
gray = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

% figure
% plot(bk(0:900, bk00, bk01, bk02), 'k', 'LineWidth', 2), axis([0 900 1.5 5.5])
% yyaxis right
% plot(bk(0:900, bk10, bk11, bk12), 'Color', gray, 'LineWidth', 2), axis([0 900 -4.5 0.5])
% title('bk0 & bk1')
% h1 = gcf;
% h1.Position = [0 648 500 155];
% 
% figure
% plot(theta_1(0:900), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
% xlim([0 900])
% title('theta1')
% h2 = gcf;
% h2.Position = [0 1 500 155];
% 
% figure
% plot(theta_2(0:900), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
% xlim([0 900])
% title('theta2')
% h3 = gcf;
% h3.Position = [0 248 500 155];

% figure
disp(C1(500, Ts))
disp(C2n(500, Ts))
disp(theta_2(500, Ts))

figure
plot(bk(0:900*Ts, bk00, bk01, bk02, Ts), 'k', 'LineWidth', 2), axis([0 900*Ts 1.5 5.5])
yyaxis right
plot(bk(0:900*Ts, bk10, bk11, bk12, Ts), 'Color', gray, 'LineWidth', 2), axis([0 900*Ts -4.5 0.5])

figure
subplot(5,1,1)
plot(bk(0:900*Ts, bk00, bk01, bk02, Ts), 'k', 'LineWidth', 2), axis([0 900*Ts 1.5 5.5])
yyaxis right
plot(bk(0:900*Ts, bk10, bk11, bk12, Ts), 'Color', gray, 'LineWidth', 2), axis([0 900*Ts -4.5 0.5])
subplot(5,1,2)
plot(theta_2(0:900*Ts, Ts), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
xlim([0 900*Ts])
subplot(5,1,3)
plot(theta_1(0:900*Ts, Ts), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
xlim([0 900*Ts])
subplot(5,1,4)
plot(C1(0:900*Ts, Ts), 'LineWidth', 2), axis([0 900*Ts 400 1300])
yyaxis right
plot( T(0:900*Ts, Ts), 'LineWidth', 2), axis([0 900*Ts 300  500])
subplot(5,1,5)
% plot(C2n(0:900)), axis([0 900 -1.5 1.5])
plot(C2n(0:900*Ts, Ts)), axis([0 900*Ts 140 280])

% figure
% plot(C1(0:900), 'LineWidth', 2), axis([0 900 400 1300])
% yyaxis right
% plot( T(0:900), 'LineWidth', 2), axis([0 900 300  500])
% title('C1 & T')

% figure
% plot(C2n(0:900)), axis([0 900 -1.5 1.5])
% title('C2n')

Tf = 900*Ts;
t  = 0:Ts:Tf;
Nt = length(t);
C2 = zeros(Nt+1, 1);
Q1 = zeros(Nt, 1);
% Q1 = ones(Nt, 1)*0.009;
e  = zeros(Nt, 1);
r  = ones(Nt, 1)*150;
e1 = 0;
Q1_1 = 0;
C2_1 = 0;
for k = 1:Nt
%     C2(k) = -C2_1*(1500*theta_1(t(k), Ts)-1) + Q1_1*12*(C1(t(k), Ts)-C2_1)
%     C2(k) = -C2_1*(1500*exp(-30/(0.008*480))-1) + Q1_1*12*(500-C2_1);
    C2(k) = -C2_1*(1500*exp(-30/(0.008*350))-1) + Q1_1*12*(1200-C2_1);
    e(k)  = r(k) - C2(k);
%     Q1(k) = Q1_1 + e(k)*bk(t(k), bk10, bk11, bk12, Ts) + e1*bk(t(k), bk00, bk01, bk02, Ts);
%     Q1(k) = Q1_1 + e(k)*0.0001000 - e1*.00001000; % melhor para T=480, C1 = 500
%     Q1(k) = Q1_1 + e(k)*0.0000519 - e1*.00005000; % melhor para T=350, C1 = 1200
    disp(['C2 = ' num2str(C2(k))])
    disp(['e = ' num2str(e(k))])
    disp(['Q1 = ' num2str(Q1(k))])
    disp('---')

%     dC = Q1(k)/5*(C1(t(k)-C2(k))) - 25*exp(-30/(0.008*T(t(k))))*C2(k);
%     C2(k+1) = C2(k) + dC*Ts;
    
    e1 = e(k);
    Q1_1 = Q1(k);
    C2_1 = C2(k);
end

figure
% plot(t, C2(1:end-1))
% axis([0 t(end) 0 220])
plot(C2(1:end-1))
axis([0 Tf/Ts 0 220])
figure
plot(Q1)

function x = theta_1(t, Ts)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = exp(-30/(0.008*T(t(i), Ts)));
    end
end

function x = theta_2(t, Ts)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = C1(t(i), Ts)/5000 - 3*(C2n(t(i), Ts)/250);
%         x(i) = 12*(C1(t(i), Ts) - C2n(t(i), Ts));
    end
end

function x = T(t, Ts)
    N = length(t);
    x = zeros(N, 1);
    vecT = [347 400 480];
    vect = [150 250 420 600]*Ts;
    for i = 1:N
        if t(i) < vect(1)
            x(i) = vecT(1);
        elseif t(i) < vect(2)
            A = [vect(1) 1
                 vect(2) 1];
            B = [vecT(1)
                 vecT(2)];
            r = A\B;
            a = r(1);
            b = r(2);
            x(i) = a*t(i)+b;
        elseif t(i) < vect(3)
            x(i) = vecT(2);
        elseif t(i) < vect(4)
            A = [vect(3) 1
                 vect(4) 1];
            B = [vecT(2)
                 vecT(3)];
            r = A\B;
            a = r(1);
            b = r(2);
            x(i) = a*t(i)+b;
        else
            x(i) = vecT(3);
        end
    end
end

function x = C1(t, Ts)
    N = length(t);
    x = zeros(N, 1);
    vect  = [150 250 400 500 600 800]*Ts;
    vecC1 = [1200 800 1000 500];
    for i = 1:N
        if t(i) < vect(1)
            x(i) = vecC1(1);
        elseif t(i) < vect(2)
            A = [vect(1) 1
                 vect(2) 1];
            B = [vecC1(1)
                 vecC1(2)];
            r = A\B;
            a = r(1);
            b = r(2);
            x(i) = a*t(i)+b;
        elseif t(i) < vect(3)
            x(i) = vecC1(2);
        elseif t(i) < vect(4)
            A = [vect(3) 1
                 vect(4) 1];
            B = [vecC1(2)
                 vecC1(3)];
            r = A\B;
            a = r(1);
            b = r(2);
            x(i) = a*t(i)+b;
        elseif t(i) < vect(5)
            x(i) = vecC1(3);
        elseif t(i) < vect(6)
            A = [vect(5) 1
                 vect(6) 1];
            B = [vecC1(3)
                 vecC1(4)];
            r = A\B;
            a = r(1);
            b = r(2);
            x(i) = a*t(i)+b;
            
        else
            x(i) = vecC1(4);
        end
    end
end

function x = C2n(t, Ts)
    N = length(t);
    x = zeros(N, 1);
    vect  = [200 400 550 650 800]*Ts;
%     vecC2 = [  1, -1,  1, -1,  1, -1];
    vecC2 = [  270, 150,  270, 150,  270, 150];
    for i = 1:N
        if t(i) <= 0
            x(i) = 0;
        elseif t(i) < vect(1)
            x(i) = vecC2(1);
        elseif t(i) < vect(2)
            x(i) = vecC2(2);
        elseif t(i) < vect(3)
            x(i) = vecC2(3);
        elseif t(i) < vect(4)
            x(i) = vecC2(4);
        elseif t(i) < vect(5)
            x(i) = vecC2(5);
        else
            x(i) = vecC2(6);
        end
    end
end

function x = bk(t, bk10, bk11, bk12, Ts)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = bk10 + bk11*theta_1(t(i), Ts) + bk12*theta_2(t(i), Ts);
    end
end







