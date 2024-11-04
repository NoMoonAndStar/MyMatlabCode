VaractorCj0 = 8.47e-12;
VaractorM = 70;
VaractorVj0 = 80;
VStatic = 3;
CVS = VaractorCj0 / (1 + VStatic / VaractorVj0) ^ VaractorM; % 变容管静态电容
Z0 = 50;
f = 1e9;
T = [1, 0;
     Z0 / (1i * f * 2 * pi * CVS), 1];
S11 = 1 / sum(T, "all") * (T(1, 1) + T(1, 2) - T(2, 1) - T(2, 2));
S12 = 1 / sum(T, "all") * det(T) * 2;
S21 = 1 / sum(T, "all") * 2;
S22 = 1 / sum(T, "all") * (-T(1, 1) + T(1, 2) - T(2, 1) + T(2, 2));
SParam = {'S11'; 'S12'; 'S21'; 'S22'};
absS = abs([S11, S12, S21, S22])';
angleS = angle([S11, S12, S21, S22])' * 180 / pi;
table(SParam, absS, angleS)
