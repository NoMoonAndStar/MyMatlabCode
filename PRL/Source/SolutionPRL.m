%%复现论文

%% Initialization
clear
clc
close all

%% Parameters Definition
g = 1; %空间调制波数
c = 3e8; %光速
Omega = c * g; %时间调制频率
eps0 = 8.85e-12; %真空介电常数
mu0 = 4 * pi * 1e-7; %真空磁导率
N = 15; %平面波展开阶数
alpha = 0.06; %调制强度

%% Preprocessing
Nt = (-N:N)'; %展开阶数向量
Nnum = length(Nt); %展开后总波数
OMEGA = linspace(-2 * Omega, 2 * Omega, 10000); %扫描频率
omegaNum = length(OMEGA); %扫描频率总数
kv = zeros(2 * Nnum, omegaNum); %真空中特征值矩阵
Hv = zeros(2 * Nnum, omegaNum); %真空中特征向量
Mv = zeros(2 * Nnum, 2 * Nnum, omegaNum); %真空中特征向量组
kv_forward = zeros(2 * Nnum, omegaNum); %前向波特征值
kv_backward = zeros(2 * Nnum, omegaNum); %后向波特征值
km = zeros(2 * Nnum, omegaNum); %超材料中特征值
Hm = zeros(2 * Nnum, 2 * Nnum, omegaNum); %超材料中特征向量
evtra = zeros(Nnum, 1); %透射波系数向量

%% Caculating eigenvalue
loop = 1; %循环次数

for omega = OMEGA

    %% 真空中的特征值
    kv(1:Nnum, loop) = -Nt * g - sqrt(c ^ -2 * (omega + Nt * Omega) .^ 2); %前半部分特征值
    kv(Nnum + 1:2 * Nnum, loop) = -Nt * g + sqrt(c ^ -2 * (omega + Nt * Omega) .^ 2); %后半部分特征值
    Hv(1:Nnum, loop) =- (Nt * g + kv(1:Nnum, loop)) ...
        ./ (mu0 * (omega + Nt * Omega)); %前半部分特征向量
    Hv(Nnum + 1:2 * Nnum, loop) =- (Nt * g + kv(Nnum + 1:2 * Nnum, loop)) ...
        ./ (mu0 * (omega + Nt * Omega)); %后半部分特征向量
    kv_forward(:, loop) = kv(:, loop);
    kv_forward(Hv(:, loop) >= 0, loop) = NaN; %前向波特征值
    kv_backward(:, loop) = kv(:, loop);
    kv_backward(Hv(:, loop) < 0, loop) = NaN; %后向波特征值
    Mv(:, :, loop) = [diag(ones(1, Nnum), 0), diag(ones(1, Nnum), 0); ...
                          diag(Hv(1:Nnum, loop), 0), diag(Hv(Nnum + 1:2 * Nnum, loop), 0)]; %真空中特征向量组

    %% 超材料中特征值
    M = HE(alpha, omega, N, g, Omega, eps0, mu0); %系数矩阵
    [Vers, Vals] = eig(M);
    k_eig = diag(Vals); %超材料特征值
    real_km = real(k_eig);
    [Y, I] = sort(real_km); %特征值排序，Y为排序后的特征值，I为排序序列
    k_eig = k_eig(I);
    Vers = Vers(:, I); %重排特征向量
    km(:, loop) = k_eig;
    Hm(:, :, loop) = Vers;

    loop = loop + 1;
end

%% Plot k-omega relation in metamaterial
omegaPlot = OMEGA / Omega;
kPlot = km;
kImag = real(kPlot) .* (imag(kPlot) ~= 0);
kImag(kImag == 0) = NaN;
h1 = figure;
hold on
grid on
plot(real(kPlot), omegaPlot, 'g', 'LineWidth', 1); %批量绘图
plot(kImag, omegaPlot, 'ro', 'LineWidth', 0.1); %特别显示出虚部不为0的部分
xlabel('beta')
ylabel('Normalizedomega')
title('Dispersion in Luminal Material')

%% Plot k-omega relation in Vaccumn
h2 = figure;
hold on
grid on
plot(kv, omegaPlot, 'ko', 'LineWidth', 1);
plot(kv_forward, omegaPlot, 'r', 'LineWidth', 1);
plot(kv_backward, omegaPlot, 'g', 'LineWidth', 1);
xlabel('beta')
ylabel('Normalizedomega')
title('Dispersion in Vaccumn')

%% Caculating transparent field
NormOmega = 0.69; %入射频率
Normomega = OMEGA / Omega;
[~, omegaIndex] = min(abs(Normomega - NormOmega));
omega = Normomega(omegaIndex) * Omega;

%处理真空中的场
Mv_ = Mv(:, :, omegaIndex);
Hv_ = Hv(:, omegaIndex);

if sum(isnan(Hv_))
    Hv_(isnan(Hv_)) = (-1) .^ (1:sum(isnan(Hv_))) * 0.0027;
end

Mvinc = Mv_(:, Hv_ < 0);
Mvref = Mv_(:, Hv_ > 0);
%再处理超材料中的场
km_ = km(:, omegaIndex);
Mm = Hm(:, :, omegaIndex);

d = [0, 15, 30, 45, 60]; %超材料的长度
P = zeros(2 * Nnum, 2 * Nnum, length(d)); %代表超材料内的相移
A = zeros(2 * Nnum, Nnum, length(d)); %辅助矩阵A
B = zeros(2 * Nnum, Nnum, length(d)); %辅助矩阵B
etra = zeros(Nnum, length(d)); %透射场幅度向量
evref = zeros(Nnum, length(d)); %反射场幅度向量
Em2 = zeros(2 * Nnum, length(d)); %超材料右侧幅度向量

for k = 1:length(d)
    P(:, :, k) = diag(diag(exp(1i * diag(km_, 0) * d(k))));
    A(:, :, k) = Mm * ((Mm * P(:, :, k)) \ Mvinc);
    B(:, :, k) = -Mvref;
end

AB = [A, B];

evinc = [zeros(N, 1); 1; zeros(N, 1)];

for k = 1:length(d)
    ev = AB(:, :, k) \ (Mvinc * evinc); %求解向量
    etra(:, k) = ev(1:Nnum); %透射场幅度等于前半部分向量
    evref(:, k) = ev(Nnum + 1:2 * Nnum); %反射场幅度等于后半部分向量
    Em2(:, k) = (Mm * P(:, :, k)) \ (Mvinc * etra(:, k)); %用边界条件，电场切向分量连续E
end

%% Combine time domain waveform
Ts = 1 / (1200 * Omega);
t = 0:Ts:6 * pi / Omega - Ts;
tNum = length(t);
EHi = Mvinc * evinc; %入射波系数
EHr = Mvref * evref; %反射波系数
EHt = zeros(2 * Nnum, length(d)); %透射波系数
EHtra = zeros(2 * Nnum, length(d), tNum); %谐波分量大小

for k = 1:length(d)
    EHt(:, k) = Mvinc * etra(:, k);
end

Etras = zeros(tNum, length(d)); %透射电场总大小
k1 = omega / c;
loop = 1;

for k = t
    et = exp(1i * ((k1 + Nt * g) * d - (omega + Nt * Omega) * k)); %谐波
    ett = [et; et];
    EHtra(:, :, loop) = EHt .* ett;
    Etras(loop, :) = sum(EHtra(1:Nnum, :, loop));
    loop = loop + 1;
end

Eout = real(Etras);
Eout2 = abs(Etras') .^ 2;
NormalizedEout2 = Eout2 / max(max(Eout2));

%% Draw time domain waveform
h3 = figure;
plot(t * Omega, NormalizedEout2, 'LineWidth', 1)
legend('d=0', 'd=10', 'd=20', 'd=30')
title('NormalizedEout2 N=12')
xlabel('Omegat')

epsilonr = 1 + 2 * alpha * cos((g * d)' - Omega * t);
h4 = figure;
plot(Omega * t, epsilonr, 'LineWidth', 1)
title('epsilonr')
xlabel('Omegat')

%% Draw the spectrogram
fs = 1 / Ts;
fshift = (-tNum / 2:tNum / 2 - 1) * (fs / tNum);
m = length(d) / 2;

if mod(length(d), 2) ~= 0
    m = ceil(m);
end

for i = 1:length(d)
    y = fft(Eout(:, i));
    L = length(y);
    P = abs(y / L);
    yshift = P(1:L / 2 + 1);
    yshift(2:end) = 2 * yshift(2:end);
    f = fs / L * (0:(L / 2));
    subplot(m, 2, i);
    plot(f * 2 * pi / Omega, abs(yshift))
    xlabel('omega/Omega')
    ylabel('Magnitude')
    legend(strcat('d=', num2str(d(i))))
    xlim([0 1000])
end
