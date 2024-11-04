function [TNorm, PA] = TheoryAmplified(Cj0, M, Vj0, VStatic, a, EffectiveDielectricr, Z0, ModulateStrength, fm, tA)

    % 调制电压参数
    t = linspace(0,1e-6,1e5); % 调制时间
    VModulateAC = ModulateStrength*VStatic; % 交流调制电压
    OmegaModulate = 2 * pi *fm; % 调制角频率

    % 调制电压
    VModulate = VStatic + VModulateAC * cos(OmegaModulate * t);

    % 调制电容
    CVaractor = Cj0 ./ (1 + VModulate / Vj0).^M;
    CV = Cj0 ./ (1 + VStatic / Vj0).^M * 2;

    % 基板电容
    c = 3e8; % 光速
    Vp = c / sqrt(EffectiveDielectricr);
    CMsub = a / Vp / Z0;

    % 对调制电容变化作单边傅里叶变换
    N = length(t);
    CVaractorFFT = fftshift(fft(CVaractor, N)); % 傅里叶变换
    TF = islocalmax(abs(CVaractorFFT)); 
    CMag = CVaractorFFT(TF); % 取极大值
    CMag = abs(CMag);
    CenterIndex = (length(CMag) + 1) / 2;
    CM0 = CMag(CenterIndex); % 直流幅度
    CM1 = CMag(CenterIndex + 1); % 基频幅度
    CM2 = CMag(CenterIndex + 2); % 二次谐波幅度
    CM3 = CMag(CenterIndex + 3); % 三次谐波幅度
    CM4 = CMag(CenterIndex + 4); % 三次谐波幅度
    CM1 = CM1 / CM0 * 2; % 归一化，因为有两只变容管所以要乘以2
    CM2 = CM2 / CM0 * 2;
    CM3 = CM3 / CM0 * 2;
    CM4 = CM4 / CM0 * 2;

    % 调制强度
    Mc1 = CM1 * CV / (CMsub + CV); % 等效基频调制强度
    Mc2 = CM2 * CV / (CMsub + CV); % 等效二次谐波调制强度
    Mc3 = CM3 * CV / (CMsub + CV); % 等效二次谐波调制强度
    Mc4 = CM4 * CV / (CMsub + CV); % 等效二次谐波调制强度
    disp("基频调制强度为：")
    disp(Mc1)
    % 放大倍数
    TNorm = t * fm; % t/T
    PA = exp((Mc1 + 2*Mc2 + 3*Mc3 + 4 * Mc4) * OmegaModulate * tA);
end