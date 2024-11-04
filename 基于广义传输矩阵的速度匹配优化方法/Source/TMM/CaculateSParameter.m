function S = CaculateSParameter(UnitE, F, VaractorCj0, VaractorM, VaractorVj0, VStatic, Zc1, Zc2, ...
        LSeries1, LSeries2, ZPort1, ZPort2, ZPort3, ZPort4, ZPort5, ZPort6, N, Fmin, Fmax, FSweep)
    %% 根据参数计算S参数
    arguments
        UnitE(1, 1) {mustBeNumeric, mustBeFinite} % 单元电长度，单位：弧度
        F(1, 1) {mustBeNumeric, mustBeFinite} % 规定电长度对应的频率
        VaractorCj0(1, 1) {mustBeNumeric, mustBeFinite} % 变容管的零偏压结电容
        VaractorM(1, 1) {mustBeNumeric, mustBeFinite} % PN结电压斜率系数
        VaractorVj0(1, 1) {mustBeNumeric, mustBeFinite} % PN结内建电压
        VStatic(1, 1) {mustBeNumeric, mustBeFinite} % 变容管静态偏压
        Zc1(1, 1) {mustBeNumeric, mustBeFinite} % 调制路的特性阻抗
        Zc2(1, 1) {mustBeNumeric, mustBeFinite} % 主路的特性阻抗
        LSeries1(1, 1) {mustBeNumeric, mustBeFinite} % 调制路串联集总电感值
        LSeries2(1, 1) {mustBeNumeric, mustBeFinite} % 主路串联集总电感值
        ZPort1(1, 1) {mustBeNumeric, mustBeFinite} % 端口1的特性阻抗
        ZPort2(1, 1) {mustBeNumeric, mustBeFinite} % 端口2的特性阻抗
        ZPort3(1, 1) {mustBeNumeric, mustBeFinite} % 端口3的特性阻抗
        ZPort4(1, 1) {mustBeNumeric, mustBeFinite} % 端口4的特性阻抗
        ZPort5(1, 1) {mustBeNumeric, mustBeFinite} % 端口5的特性阻抗
        ZPort6(1, 1) {mustBeNumeric, mustBeFinite} % 端口6的特性阻抗
        N(1, 1) {mustBeInteger, mustBeFinite} % 单元总数
        Fmin(1, 1) {mustBeNumeric, mustBeFinite} % 扫频最小值
        Fmax(1, 1) {mustBeNumeric, mustBeFinite} % 扫频最大值
        FSweep(1, 1) {mustBeInteger, mustBeFinite} % 扫频点数
    end

    %% 设置默认参数及计算部分参数
    c = 3e8; % 为什么认为传输线本身的相速度是光速，因为这里假设所有介质的epsr都是1
    Freq = linspace(Fmin, Fmax, FSweep); % 扫频范围
    Omega = 2 * pi * Freq;
    CVS = VaractorCj0 / (1 + VStatic / VaractorVj0) ^ VaractorM; % 变容管静态电容
    S = zeros(6, 6, FSweep);

    %% 计算散射矩阵
    loop = 1;

    for omega = Omega
        % 不保存中间结果以节省内存
        ZLs1 = 1i * omega * LSeries1;
        ZLs2 = 1i * omega * LSeries2;
        TUnitLsHorizontal = [eye(3), diag([ZLs1, ZLs2, ZLs1]); % 为什么要除2，因为这只是一个单元一半的
                                           zeros(3), eye(3)];
        theta = UnitE / F * omega / 2 * pi; % 为什么是这个，因为theta=omega/c*l
        TUnitTLHorizontal = [cos(theta / 2) * eye(3), 1i * diag([Zc1, Zc2, Zc1]) * sin(theta / 2);
                                                                 1i * diag([1 / Zc1, 1 / Zc2, 1 / Zc1]) * sin(theta / 2), cos(theta / 2) * eye(3)
                                                                 ];

        TUnitHorizontal1 = TUnitLsHorizontal * TUnitTLHorizontal; % 左边的一半的横向的传输矩阵
        TUnitHorizontal2 = TUnitTLHorizontal * TUnitLsHorizontal; % 右边的一半的横向的传输矩阵
        ZC = 1 / (1i * omega * CVS); % 变容管的阻抗，两个变容管都工作在反偏状态
        Ya = 1 / ZC;
        Yb = -1 / ZC;
        Yc = 2 / ZC;
        X = [Ya, Yb, 0; % 纵向传输矩阵
             Yb, Yc, Yb;
             0, Yb, Ya];
        TUnitVerital = [eye(3), zeros(3);
                        X, eye(3)];
        TUnit = TUnitHorizontal1 * TUnitVerital * TUnitHorizontal2;
        Ttotal = TUnit ^ N; % 周期结构的总传输矩阵

        %% 阶段性debug分界线
        % ---------------------------------------------------------------------------------------------------------------------------------------------

        Ztotal = [Ttotal(1:3, 1:3) * Ttotal(4:6, 1:3) ^ -1, Ttotal(4:6, 1:3) ^ -1; % 从传输矩阵到阻抗矩阵
                  Ttotal(4:6, 1:3) ^ -1, (Ttotal(4:6, 1:3) ^ -1) * Ttotal(4:6, 4:6)];

        ZPortSqrt = diag(1./sqrt([ZPort1, ZPort2, ZPort3, ZPort4, ZPort5, ZPort6])); % [1/sqrt(Z_{oi})]

        ZNormolized = ZPortSqrt * Ztotal * ZPortSqrt;
        S(:, :, loop) = (ZNormolized - eye(6)) * (ZNormolized + eye(6)) ^ -1; % 从归一化阻抗矩阵得到散射矩阵
        loop = loop + 1;
    end

end
