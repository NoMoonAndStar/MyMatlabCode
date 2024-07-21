%% 复现3.1和3.2的图
%初始化
clear
clc
close all

%四个电路参数
LR = 1e-8; %右手单位长度电感
CL = 1e-6; %左手倍长度电容
CR = 2e-8; %右手单位长度电感
LL = 5e-6; %左手倍长度电感

%中间量
omegaR = sqrt(1 / (LR * CR));
omegaL = sqrt(1 / (LL * CL));
k = LR * CL + LL * CR;
omegase = sqrt(1 / (LR * CL));
omegash = sqrt(1 / (LL * CR));

%%fig3.3(a)
Omega = 3e7;
omega = 0:Omega / 200000:Omega; %扫频范围
beta = zeros(1, length(omega));

for i = 1:length(omega)
    temp = (omega(i) / omegaR) ^ 2 - k * omegaL ^ 2 + (omegaL / omega(i)) ^ 2; %计算相位常数

    if temp < 0 %纯实属，阻带
        continue
    end

    if omega(i) < min(omegase, omegash)
        beta(i) = -sqrt(temp); %这里我没懂为什么要取负号，因为是左手的吗
    elseif omega(i) > max(omegase, omegash)
        beta(i) = sqrt(temp);
    end

end

%剔除虚数解
idx = beta ~= 0;
beta = beta(idx);
omegaCRLH = omega(idx);

%绘制色散曲线
h1 = figure;
hold on
plot(beta, omegaCRLH, '.', 'Color', 'r')
%画取向相反时的色散曲线
plot(-beta, omegaCRLH, '.', 'Color', 'b')
%画坐标轴
plot([0 0], [0, Omega], "LineWidth", 0.5, "Color", "k")
xlim([-1 1])
xlabel('beta')
ylabel('omega')
title('Dispersion relation')

%%fig3.3(b)
%PRH
betaPRH = omega * sqrt(LR * CR);
%PLH
betaPLH = -omegaL ./ omega;

h2 = figure;
hold on
plot(beta, omegaCRLH, '.', 'Color', 'b')
plot(betaPRH, omega, '.', 'Color', 'r')
plot(betaPLH, omega, '.', 'Color', 'g')
plot([0 0], [0, Omega], "LineWidth", 0.5, "Color", "k")
xlim([-0.6 0.6])
xlabel('beta')
ylabel('omega')
title('Dispersion relation')
legend('CRLH', 'PRH', 'PLH')
