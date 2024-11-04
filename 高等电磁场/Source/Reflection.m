%% 计算反射系数随入射角的变化
clc; 
clear; 
close all;

%% n1/n2 = 1/1.5
n1 = 1; % 媒质1的折射率
n2 = 1.5; %媒质2的折射率
mur1 = 1;
mur2 = 1;
epsr1 = n1^2 / mur1; % 媒质1的相对介电常数
epsr2 = n2^2 / mur2; % 媒质2的相对介电常数
thetai = linspace(0, pi/2, 360);

%% TE波入射
GammaTE = (cos(thetai) - sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2)) ...
    ./ (cos(thetai) + sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));
TranTE = 2 * cos(thetai) ./ (cos(thetai) + sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));

%% TM波入射
GammaTM = (cos(thetai) - sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2)) ...
    ./ (cos(thetai) + sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));
TranTM = 2 * cos(thetai) ./ (cos(thetai) + sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));


% 绘图
figure("Name", 'n1/n2=1/1.5 反射系数模')
hold on
plot(thetai*180/pi, GammaTE, 'linewidth', 1.5)
plot(thetai*180/pi, GammaTM, 'linewidth', 1.5)
xlabel('thetai')
ylabel('Gamma')
title('n1/n2=1/1.5')
grid on
legend('TE', 'TM')

figure("Name", 'n1/n2=1/1.5 反射系数模')
hold on
plot(thetai*180/pi, abs(GammaTE), 'linewidth', 1.5)
plot(thetai*180/pi, abs(GammaTM), 'linewidth', 1.5)
xlabel('thetai')
ylabel('|Gamma|')
title('n1/n2=1/1.5')
grid on
legend('TE', 'TM')

figure("Name", 'n1/n2=1/1.5 反射系数相角')
hold on
plot(thetai*180/pi, angle(GammaTE) * 180 / pi, 'linewidth', 1.5)
plot(thetai*180/pi, angle(GammaTM) * 180 / pi, 'linewidth', 1.5)
ylim([-180 180])
xlabel('thetai')
ylabel('ang(Gamma)')
grid on
title('n1/n2=1/1.5')
legend('TE', 'TM')

%% n1/n2 = 1.5/1
n1 = 1.5; % 媒质1的折射率
n2 = 1; %媒质2的折射率
mur1 = 1;
mur2 = 1;
epsr1 = n1^2 / mur1; % 媒质1的相对介电常数
epsr2 = n2^2 / mur2; % 媒质2的相对介电常数
thetai = linspace(0, pi/2, 360);

%% TE波入射
GammaTE = (cos(thetai) - sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2)) ...
    ./ (cos(thetai) + sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));
TranTE = 2 * cos(thetai) ./ (cos(thetai) + sqrt(epsr2 * mur1 / epsr1 / mur2) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));

%% TM波入射
GammaTM = (cos(thetai) - sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2)) ...
    ./ (cos(thetai) + sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));
TranTM = 2 * cos(thetai) ./ (cos(thetai) + sqrt(epsr1 * mur2 / epsr2 / mur1) * sqrt(1 - epsr1 * mur1 * sin(thetai).^2 / epsr2 / mur2));


% 绘图

% 修正
GammaTE(abs(GammaTE)>0.999) = 1;
GammaTM(abs(GammaTM)>0.999) = 1;
figure("Name", 'n1/n2=1.5/1 反射系数')
hold on
plot(thetai*180/pi, GammaTE, 'linewidth', 1.5)
plot(thetai*180/pi, GammaTM, 'linewidth', 1.5)
xlabel('thetai')
ylabel('Gamma')
title('n1/n2=1.5/1')
grid on
legend('TE', 'TM')

figure("Name", 'n1/n2=1.5/1 反射系数模')
hold on
plot(thetai*180/pi, abs(GammaTE), 'linewidth', 1.5)
plot(thetai*180/pi, abs(GammaTM), 'linewidth', 1.5)
xlabel('thetai')
ylabel('|Gamma|')
title('n1/n2=1.5/1')
grid on
legend('TE', 'TM')

figure("Name", 'n1/n2=1.5/1 反射系数相角')
hold on
plot(thetai*180/pi, angle(GammaTE) * 180 / pi, 'linewidth', 1.5)
plot(thetai*180/pi, angle(GammaTM) * 180 / pi, 'linewidth', 1.5)
ylim([-180 180])
xlabel('thetai')
ylabel('ang(Gamma)')
grid on
title('n1/n2=1.5/1')
legend('TE', 'TM')
