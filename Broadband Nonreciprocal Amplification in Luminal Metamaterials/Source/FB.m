%% 用F-B方法计算场和能量，与JOSB方法比较
clear;
clc;
close all;

%% 结构参数
epslion0 = 8.85e-12; % 真空介电常数
mu0 = 4 * pi * 1e-7; % 真空磁导率
epslionr = 1; % 未调制时的相对介电常数
c0 = 3e8; % 光速

% 调制参数
OMEGA = 0.07 * 1e9;
delta = 0; % 速度差
cg = (1 + delta) * c0 / sqrt(epslionr); % 光栅速度
g = OMEGA / cg; % 调制空间频率
alpha = 0.05; % 调制强度
d = 100; % 介质长度

%% 计算色散曲线
Omega = 1e9; % 扫频范围
OmegaNum = length(Omega); % 扫频点数
Order = 50; % 展开阶数
HarmonicNum = 2 * Order + 1; % 谐波数目（包含正负）
HarmonicSequence = -Order:Order; % 阶次序号

MEE = zeros(HarmonicNum); % 初始化特征矩阵的电场-电场部分
MEH = zeros(HarmonicNum); % 初始化特征矩阵的电场-磁场部分
MHE = zeros(HarmonicNum); % 初始化特征矩阵的磁场-电场部分
MHH = zeros(HarmonicNum); % 初始化特征矩阵的磁场-磁场部分
LuminWaveNumber = zeros(OmegaNum, 2 * HarmonicNum); %初始化相应波数
LuminEigenVec = zeros(2 * HarmonicNum, 2 * HarmonicNum, OmegaNum); % 初始化超材料中特征矩阵，每一页表示一个频点

VacEigenValue = zeros(OmegaNum, 2 * HarmonicNum); % 真空中特征值向量
VacEigenMatrix = zeros(2 * HarmonicNum, 2 * HarmonicNum, OmegaNum); % 初始化真空中特征矩阵
VacEigenH = zeros(OmegaNum, 2 * HarmonicNum); % 初始化真空中特征向量中磁场部分
VacEigenVec = zeros(2 * HarmonicNum, 2 * HarmonicNum, OmegaNum); % 初始化真空中特征向量组矩阵
ForwardVacEigenValue = zeros(OmegaNum, 2 * HarmonicNum); % 真空前向波部分
BackwardVacEigenValue = zeros(OmegaNum, 2 * HarmonicNum); % 真空后向波部分

%% 扫频计算
for loop = 1:OmegaNum
    
    % 处理真空中特征值
    VacEigenValue(loop, 1:HarmonicNum) = -HarmonicSequence * g - sqrt(c0 ^ -2 * (Omega(loop) + HarmonicSequence * OMEGA) .^ 2);
    VacEigenValue(loop, HarmonicNum + 1:2 * HarmonicNum) = -HarmonicSequence * g + sqrt(c0 ^ -2 * (Omega(loop) + HarmonicSequence * OMEGA) .^ 2);
    % 真空特征向量
    VacEigenH(loop, 1:HarmonicNum) = -(HarmonicSequence * g + VacEigenValue(loop, 1:HarmonicNum)) ./ (mu0 * (Omega(loop) + HarmonicSequence * OMEGA));
    VacEigenH(loop, HarmonicNum + 1:2 * HarmonicNum) = -(HarmonicSequence * g + VacEigenValue(loop, HarmonicNum + 1:2 * HarmonicNum)) ./ (mu0 * (Omega(loop) + HarmonicSequence * OMEGA));
    % 分离前向波和后向波部分（取决于k+ng的正负）
    ForwardVacEigenValue(loop, :) = VacEigenValue(loop, :);
    ForwardVacEigenValue(loop, find(VacEigenH(loop, :) >= 0)) = NaN; % 根据E, H, S构成右手系确定前向波
    BackwardVacEigenValue(loop, :) = VacEigenValue(loop, :);
    BackwardVacEigenValue(loop, find(VacEigenH(loop, :) < 0)) = NaN; % 根据E, H, S构成右手系确定后向波
    VacEigenVec11 = diag(ones(HarmonicNum, 1), 0);
    VacEigenVec12 = diag(ones(HarmonicNum, 1), 0);
    VacEigenVec21 = diag(VacEigenH(loop, 1:HarmonicNum), 0);
    VacEigenVec22 = diag(VacEigenH(loop, HarmonicNum + 1:2 * HarmonicNum), 0);
    VacEigenVec(:, :, loop) = [VacEigenVec11, VacEigenVec12; VacEigenVec21, VacEigenVec22];
    % 处理Lumin介质中的特征值
    % 构造Lumin介质的特征矩阵
    for i = 1:HarmonicNum
        for j = 1:HarmonicNum
            if i == j
                MEE(i, j) = - HarmonicSequence(i) * g;
                MEH(i, j) = - mu0 * (Omega(loop) + HarmonicSequence(i) * OMEGA);
                MHE(i, j) = - epslionr * epslion0 * (Omega(loop) + HarmonicSequence(i) * OMEGA);
                MHH(i, j) = - HarmonicSequence(i) * g; 
            elseif i - j == -1 % 主对角线下方
                MHE(i, j) = - epslionr * epslion0 * ((Omega(loop) + HarmonicSequence(i) * OMEGA) * alpha + OMEGA * alpha);
            elseif i - j == 1 % 主对角线上方
                MHE(i, j) = - epslionr * epslion0 * ((Omega(loop) + HarmonicSequence(i) * OMEGA) * alpha - OMEGA * alpha);
            end
        end
    end
    MHE = MHE'; % ??
    LuminEigenMatrix = [MEE, MEH; MHE, MHH]; % 特征矩阵
    [V, D] = eig(LuminEigenMatrix); % 解特征矩阵的特征值和特征向量
    k = diag(D);
    Rek = real(k);
    [Y, I] = sort(Rek);
    k = k(I); % 按顺序重排特征值
    V = V(:, I); % 重排特征向量，每一列对应于一个特征值
    LuminWaveNumber(loop, :) = k;
    LuminEigenVec(:, :, loop) = V;
end

%% 绘制Lumin介质能带图
if false
    figure('name', '能带');
    hold on
    % 每一行对应于一个频率，每一列对应于一条线
    for band = LuminWaveNumber
        for i = 1:OmegaNum
            if imag(band(i)) == 0
                plot(Omega(i)/OMEGA, real(band(i)), 'b.') % 波数是纯实数，属于通带
            else
                plot(Omega(i)/OMEGA, real(band(i)), 'r.') % 波数不是纯实数，属于禁带
            end
        end
    end
    xlabel('Omega/OMEGA')
    ylabel('k')
    title(strcat('Lumin介质能带，delta=', num2str(delta)))
end

%%--------------------------------分界线---------------------------------------%%
%% 传输矩阵法解场分布，此处考虑一个有限长度，无限宽度的TEM波导
%% 输入参数
OmegaInput = 1e9;
NormOmegaInput = OmegaInput / OMEGA; % 归一化入射角频率
% 找出对应的波数
[~, OmegaInputIndex] = min(abs(Omega/OMEGA-NormOmegaInput)); % 该输入频率对应的索引
OmegaInput = Omega(OmegaInputIndex); % 实际使用入射角频率

% 真空的特征向量矩阵
VacEigenVec = VacEigenVec(:, :, OmegaInputIndex);

% Lumin介质中的特征向量
LuminWaveNumber = LuminWaveNumber(OmegaInputIndex, :);
LuminEigenVec = LuminEigenVec(:, :, OmegaInputIndex);
Mm = LuminEigenVec; % 与论文中符号一致

%% 提取传输矩阵
% 提取真空正向和反向传播的特征向量
Mvinc = VacEigenVec(:, find(VacEigenH(OmegaInputIndex, :) < 0)); % 正向波的特征向量矩阵
Mvref = VacEigenVec(:, find(VacEigenH(OmegaInputIndex, :) > 0)); % 反向波的特征向量矩阵
VIncIndex = find(VacEigenH(OmegaInputIndex, :) < 0); % 正向波原来的索引
HValue = (VacEigenH(OmegaInputIndex, VIncIndex(1)));
[row, col] = find(abs(VacEigenVec - HValue) < 1e-16);
FundamentalIndex = find(VIncIndex == col(find(row == HarmonicNum + (HarmonicNum + 1) / 2))); % 正向基波的在原特征向量矩阵中的索引

% 传输矩阵
P = diag(exp(1i * LuminWaveNumber * d));

% 辅助矩阵
evinc = zeros(size(Mvinc, 2), 1); % 初始化入射相位向量
evinc(FundamentalIndex) = 1; % 基波幅度为1
A = Mm * (Mm * P) ^ -1 * Mvinc;
B = -Mvref;

% 解边界条件方程得相位向量
evtraref = ([A, B]) ^ -1 * (Mvinc * evinc);
evtra = evtraref(1:HarmonicNum);
evref = evtraref(HarmonicNum + 1:2 * HarmonicNum);

%%--------------------------------分界线---------------------------------------%%
%% 获取场分布
k0 = OmegaInput / c0;
Einc = Mvinc * evinc; % 入射场复矢量
Einc = Einc(1:HarmonicNum); % 只取电场部分
Eref = Mvref * evref; % 反射场复矢量
Eref = Eref(1:HarmonicNum);
Etra = Mvinc * evtra; % 透射场复矢量
Etra = Etra(1:HarmonicNum);

% 时间因子
tnum = 500;
Omegat = linspace(0, 2 * pi, tnum)'; % 观察2个调制周期

% 合成场
Einct = zeros(1, tnum);
Ereft = zeros(1, tnum);
Etrat = zeros(1, tnum);
for i = 1:tnum
    Einct(i) = sum(Einc' .* exp(-1i * (OmegaInput / OMEGA + HarmonicSequence) * Omegat(i)));
    Ereft(i) = sum(Eref' .* exp(-1i * (OmegaInput / OMEGA + HarmonicSequence) * Omegat(i)));
    Etrat(i) = sum(Etra' .* exp(-1i * ((k0 + HarmonicSequence * g) * d - (OmegaInput / OMEGA + HarmonicSequence) * Omegat(i))));
end

Einct2 = abs(Einct) .^ 2;
Ereft2 = abs(Ereft) .^ 2;
Etrat2 = abs(Etrat) .^ 2;

h1 = figure('Name', '入射，透射和反射场');
hold on
plot(Omegat, Einct, 'b', 'linewidth', 2);
plot(Omegat, Ereft, 'k', 'linewidth', 2);
plot(Omegat, Etrat, 'r', 'linewidth', 2);
xlabel('t/T')
ylabel('V/m')
legend('入射场', '反射场', '透射场')
hold off

h2 = figure('Name', '入射，透射和反射场模平方');
hold on
plot(Omegat, Einct2 / max(Einct2), 'b', 'linewidth', 2);
plot(Omegat, Ereft2 / max(Ereft2), 'k', 'linewidth', 2);
plot(Omegat, Etrat2 / max(Etrat2), 'r', 'linewidth', 2);
xlabel('t/T')
ylabel('(V/m)^2')
legend('入射场', '反射场', '透射场')
hold off

%% 超材料内部场
xnum = 5000;
em = (Mm ^ -1) * (Mvinc * evinc + Mvref * evref);
x = linspace(0, d, xnum)'; % 距离因子
Px = exp(1i * LuminWaveNumber .* x); % 相移矩阵，每一行代表一个空间点
Em = zeros(xnum, HarmonicNum); % 超材料内电场复矢量，每一行代表一个空间点
Hm = zeros(xnum, HarmonicNum); % 超材料内电场复矢量，每一行代表一个空间点
Elumin = zeros(xnum, tnum);
Hlumin = zeros(xnum, tnum);
eps = epslionr * (1 + 2 * alpha * cos(Omegat))';
for i = 1:xnum
    EHtemp = Mm * diag(Px(i, :)) * em;
    Em(i, :) = EHtemp(1:HarmonicNum);
    Hm(i, :) = EHtemp(HarmonicNum + 1:2 * HarmonicNum);
end

for i = 1:tnum
    Elumin(:, i) = sum(Em .* exp(-1i * ((OmegaInput / OMEGA + HarmonicSequence) * Omegat(i))), 2); % 超材料内的场
    Hlumin(:, i) = sum(Hm .* exp(-1i * ((OmegaInput / OMEGA + HarmonicSequence) * Omegat(i))), 2); % 超材料内的场
end

%% 绘制动画

% h3 = figure('name', '不同空间位置处的场');
% xlabel('x')
% ylabel('V/m')
% for i = 1:xnum
%     plot(Omegat/2/pi, real(Elumin(i, :)), 'b', 'LineWidth', 2)
%     xlabel('Omegat')
%     ylabel('V/m')
%     title(strcat('x=', num2str(x(i))))
%     pause(0.01)
% end

% h3 = figure('name', '不同空间位置处的场');
% for i = 1:xnum
%     subplot(1,2,1)
%     plot(Omegat/2/pi, 1 / 2 * eps * epslion0 .* abs(Elumin(i, :)).^2, 'b', 'LineWidth', 2)
%     xlabel('Omegat')
%     ylabel('J/m^3')
%     title(strcat('x=', num2str(x(i))))
%     subplot(1,2,2)
%     plot(Omegat/2/pi, 1 / 2 * mu0 * abs(Hlumin(i, :)).^2, 'r', 'LineWidth', 2)
%     xlabel('Omegat')
%     ylabel('J/m^3')
%     title(strcat('x=', num2str(x(i))))
%     pause(0.01)
% end
% 

% h4 = figure('name', '不同时间的空间场');
% xlabel('x')
% ylabel('V/m')
% for i = 1:tnum
%     plot(x, Elumin(:, i), 'b', 'LineWidth', 2)
%     xlabel('x')
%     ylabel('V/m')
%     title(strcat('Omegat=', num2str(Omegat(i))))
%     ylim([-10, 10])
%     pause(0.01)
% end

%% 复现图

% d1 = 37.5; % 对应于tau
% d2 = 150;
% eps = epslionr * (1 + 2 * alpha * cos(Omegat))';
% Et1 = Elumin(round(d1 / (d / xnum)), :);
% Et2 = Elumin(round(d2 / (d / xnum)), :);

% load 'UA.mat'
% load 'UB.mat'
% h4 = figure;
% hold on
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(Et1).^2, 'k--', 'LineWidth', 2)
% plot((1:length(UA))/length(UA), UA, 'r-', 'LineWidth', 2)
% xlabel('gX')
% ylabel('U(X)')
% grid on
% title('\tau=37.5')
% 
% 
% h5 = figure;
% hold on
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(Et2).^2, 'k--', 'LineWidth', 2)
% plot((1:length(UB))/length(UB), UB, 'r-', 'LineWidth', 2)
% xlabel('gX')
% ylabel('U(X)')
% grid on
% title('\tau=150')

% h6 = figure;
% hold on
% plot((1:size(UD,2))/size(UD,2), UD, 'r-', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(dm005).^2, 'k--', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(dm0025).^2, 'k--', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(d0).^2, 'k--', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(d0025).^2, 'k--', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, eps.*abs(d005).^2, 'k--', 'LineWidth', 2)
% xlim([0.5 1])
% xlabel('gX')
% ylabel('U(X)')
% grid on
% title('不同速度差下的能量密度')

% load 'E.mat'
% h7 = figure;
% hold on
% plot((1:length(E))/length(E)+0.003, real(E)/max(real(E)), 'r-', 'LineWidth', 2)
% plot((- Omegat + 2 * pi)/2/pi, real(Et2)/max(real(Et2)), 'k--', 'LineWidth', 2)
% xlim([0.65 0.85])
% ylim([-1 1])
% xlabel('gX')
% ylabel('V/m')
% title('F-B')

%% 与理论方法比较
% t = d / (c0 / sqrt(epslionr));
% FirstU = exp(-2 * alpha * OMEGA * t * sin(Omegat));
% SecondU = exp(-2 * alpha * OMEGA * t * sin(Omegat) - ((alpha * OMEGA * t) ^ 2) * (cos(Omegat) .^ 2));
% CorFirstU = exp(- alpha * OMEGA * t * sin(Omegat));
% CorSecondU = exp(-alpha * OMEGA * t * sin(Omegat) - 1 / 2 *((alpha * OMEGA * t) ^ 2) * (cos(Omegat) .^ 2));
% 
% h7 = figure;
% hold on
% plot((- Omegat + 2 * pi)/2/pi, FirstU, 'b-', 'LineWidth', 1) % 把X坐标变为Omegat
% plot((- Omegat + 2 * pi)/2/pi, SecondU, 'g-') % 把X坐标变为Omegat
% plot(Omegat/2/pi, abs(Elumin(end, :)).^2, 'r-', 'LineWidth', 1)
% plot((- Omegat + 2 * pi)/2/pi, CorFirstU, 'b--', 'LineWidth', 1) % 把X坐标变为Omegat
% plot((- Omegat + 2 * pi)/2/pi, CorSecondU, 'g--') % 把X坐标变为Omegat
% legend('未修正零阶解','未修正一阶解','F-B方法', '修正零阶解', '修正一阶解', 'AutoUpdate', 'off')
% plot([0.25 0.25], [0, 6], 'k:', 'LineWidth', 0.75)
% plot([0.75 0.75], [0, 6], 'k-.', 'LineWidth', 1)
% text(0.22, 5.7, 'gain', 'FontSize', 12)
% text(0.72, 0.5, 'loss', 'FontSize', 12)
% xlabel('\Omegat')
% ylabel('|E_o_u_t|^2')
% grid on

%% 放大倍数曲线
PA = max(abs(Elumin).^2, [], 2);
plot(x/cg*OMEGA/2/pi * 2, PA)
xlim([0 7.5])