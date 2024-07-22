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
omega = 0:Omega / 20000:Omega; %扫频范围
domega = omega(2) - omega(1);
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

%%fig3.4
Zc = sqrt(((omega / omegase) .^ 2 - 1) ./ ((omega / omegash) .^ 2 - 1)); %特性阻抗
ZL = sqrt(LL / CL);
redx = imag(Zc) == 0;
imdx = imag(Zc) ~= 0;
OmegaZc = omega(redx);
OmegaXc = omega(imdx);
Xc = imag(Zc(imdx));
Zc = Zc(redx);

h3 = figure;
hold on
plot(Zc / ZL, OmegaZc, '.', 'Color', 'b')
plot(Xc / ZL, OmegaXc, '--', 'Color', 'r')
plot([omegash / omegase, omegash / omegase], [0, Omega], 'LineWidth', 0.5, 'Color', 'k')
plot([0 5], [omegash, omegash], 'LineWidth', 0.5, 'Color', 'k')
xlabel('Zc/ZL')
ylabel('omega')
legend('Zc', 'Xc')
xlim([0 5])
title('Characteristic impedance')

%%fig3.5
vpCRLH = omegaCRLH ./ beta;
vpPRH = omega ./ betaPRH;
vpPLH = omega ./ betaPLH; %相速度
vgCRLH = domega ./ gradient(beta);
vgPRH = omega ./ betaPRH;
vgPLH = omega .^ 2 / omegaL; %群速度
h4 = figure;
hold on
plot(omega, vpPLH, '.', 'Color', 'b')
plot(omega, vgPLH, '.', 'Color', 'r')
xlabel('omega')
ylabel('vp, vg')
legend('vpPLH', 'vgPLH')
xlim([0 Omega])
title('Phase velocity and group velocity')

h5 = figure;
hold on
plot(omegaCRLH, vpCRLH / omegaR, '.', 'Color', 'b')
plot(omegaCRLH, vgCRLH / omegaR, '.', 'Color', 'r')
plot([0 Omega], [0 0], 'LineWidth', 0.5, 'Color', 'k')
plot([min(omegase, omegash) min(omegase, omegash)], [-1, 0], '--', 'Color', 'k')
plot([max(omegase, omegash) max(omegase, omegash)], [0, 2], '--', 'Color', 'k')
xlabel('omega')
ylabel('vp/omegaR, vg/omegaR')
legend('vpCRPLH', 'vgCRLH')
xlim([0 Omega])
ylim([-1 2])
title('Phase velocity and group velocity of CRLH')

%%fig3.7
eps0 = 8.85418782e-12; %真空介电常数
mu0 = 4 * pi * 1e-7; %真空磁导率
epsr = (CR - 1 ./ (omega .^ 2 * LL)) / eps0; %等效介电常数
mur = (LR - 1 ./ (omega .^ 2 * CL)) / mu0; %等效磁导率

h6 = figure;
hold on
subplot(2, 1, 1)
hold on
plot(omega, epsr, '.', 'Color', 'b')
plot([0 Omega], [0 0], 'LineWidth', 0.5, 'Color', 'k')
plot([0 Omega], [CR / eps0 CR / eps0], '--', 'LineWidth', 0.5, 'Color', 'k')
plot([omegash omegash], [-1e3 0], '--', 'LineWidth', 0.5, 'Color', 'k')
xlabel('omega')
ylabel('epsr')
title('Relative permittivity and relative permeability')
xlim([0 Omega])
ylim([-1e3 3e3])
subplot(2, 1, 2)
hold on
plot(omega, mur, '.', 'Color', 'r')
plot([0 Omega], [0 0], 'LineWidth', 0.5, 'Color', 'k')
plot([0 Omega], [LR / mu0 LR / mu0], '--', 'LineWidth', 0.5, 'Color', 'k')
plot([omegase omegase], [-1e3 0], '--', 'LineWidth', 0.5, 'Color', 'k')
xlabel('omega')
ylabel('mur')
xlim([0 Omega])
ylim([-0.01 0.01])

eta = sqrt(epsr .* mur); %等效折射率
n0 = 3e8 * sqrt(LR * CR);
LHidx = omega < min(omegase, omegash);
eta(LHidx) = -eta(LHidx);
realidx = imag(eta) == 0;
imagidx = imag(eta) ~= 0;
realeta = eta(realidx);
imageta = -imag(eta(imagidx));
omegarealteta = omega(realidx);
omegaimageta = omega(imagidx);

h7 = figure;
hold on
plot(omegarealteta, realeta, '.', 'Color', 'b')
plot(omegaimageta, imageta, '--', 'Color', 'r')
plot([0 Omega], [0 0], 'LineWidth', 0.5, 'Color', 'k')
plot([omegash omegash], [-1e3 0], '--', 'LineWidth', 0.5, 'Color', 'k')
plot([omegase omegase], [-1e3 0], '--', 'LineWidth', 0.5, 'Color', 'k')
plot([0 Omega], [n0 n0], '--', 'LineWidth', 0.5, 'Color', 'k')
ylim([-10 5])
xlabel('omega')
ylabel('eta')
legend('realeta', 'imageta')
title('Effective Refractive index')
