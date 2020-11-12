clc, clear all, close all

A = [1 theta_1(100) theta_2(100)
     1 theta_1(300) theta_2(300)
     1 theta_1(900) theta_2(900)];
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
disp(C1(500))
disp(C2n(500))
disp(theta_2(500))

figure
subplot(5,1,1)
plot(bk(0:900, bk00, bk01, bk02), 'k', 'LineWidth', 2), axis([0 900 1.5 5.5])
yyaxis right
plot(bk(0:900, bk10, bk11, bk12), 'Color', gray, 'LineWidth', 2), axis([0 900 -4.5 0.5])
subplot(5,1,2)
plot(theta_2(0:900), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
xlim([0 900])
subplot(5,1,3)
plot(theta_1(0:900), 'LineWidth', 2)%, axis([0 900 -1.5 1.5])
xlim([0 900])
subplot(5,1,4)
plot(C1(0:900), 'LineWidth', 2), axis([0 900 400 1300])
yyaxis right
plot( T(0:900), 'LineWidth', 2), axis([0 900 300  500])
subplot(5,1,5)
% plot(C2n(0:900)), axis([0 900 -1.5 1.5])
plot(C2n(0:900)), axis([0 900 140 280])




% figure
% plot(C1(0:900), 'LineWidth', 2), axis([0 900 400 1300])
% yyaxis right
% plot( T(0:900), 'LineWidth', 2), axis([0 900 300  500])
% title('C1 & T')

% figure
% plot(C2n(0:900)), axis([0 900 -1.5 1.5])
% title('C2n')

function x = theta_1(t)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = exp(-30/(0.008*T(t(i))));
    end
end

function x = theta_2(t)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
%         x(i) = C1(t(i))/5000 - 3*(C2n(t(i))/250);
        x(i) = 12*(C1(t(i)) - C2n(t(i)));
    end
end

function x = T(t)
    N = length(t);
    x = zeros(N, 1);
    vecT = [347 400 480];
    vect = [150 250 420 600];
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

function x = C1(t)
    N = length(t);
    x = zeros(N, 1);
    vect  = [150 250 400 500 600 800];
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

function x = C2n(t)
    N = length(t);
    x = zeros(N, 1);
    vect  = [200 400 550 650 800];
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

function x = bk(t, bk10, bk11, bk12)
    N = length(t);
    x = zeros(N, 1);
    for i = 1:N
        x(i) = bk10 + bk11*theta_1(t(i)) + bk12*theta_2(t(i));
    end
end


















