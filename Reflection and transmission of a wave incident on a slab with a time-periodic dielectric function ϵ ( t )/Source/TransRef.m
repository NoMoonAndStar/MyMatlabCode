%% 解透射场和反射场
clear
clc
close all

%% 定义参数
eps0 = 5.25;
Deltaeps = 0.085;
norepsm1 = 1i * 0.5 * Deltaeps / eps0;
noreps1 = -1i * 0.5 * Deltaeps / eps0;
eps1 = 1;
eps2 = 1;
LN = 8;
N = 2;
Num = 2 * N + 1;
fNum = 4001;
En0 = zeros(fNum, 3);
En1 = zeros(fNum, 3);
Enm1 = zeros(fNum, 3);
En2 = zeros(fNum, 3);
Enm2 = zeros(fNum, 3);

%% 扫频解透射场和反射场的基频幅度
loop = 1;

for nomega = linspace(0, 10, fNum)
    %% 解特征值获得展开系数
    D = zeros(Num);

    for i = 1:Num

        for j = 1:Num

            if i - j == 0
                D(i, j) = (nomega - (i - (Num + 1) / 2)) ^ 2;
            elseif i - j == 1
                D(i, j) = noreps1 * (nomega - (i - (Num + 1) / 2)) ^ 2;
            elseif i - j == -1
                D(i, j) = norepsm1 * (nomega - (i - (Num + 1) / 2)) ^ 2;
            end

        end

    end

    [V, K] = eig(D);
    K = sqrt(diag(K))';

    %第一个方程系数
    p1 = [V, V];
    %第二个方程系数
    p2 = [V, -V];
    p2 = [K, K] * sqrt(eps0) .* p2;
    p2 = 1 ./ (nomega - (-N:N))' .* p2;
    %第三个方程系数
    p3 = [V, V];
    u = [exp(1i * K * LN), exp(-1i * K * LN)];
    p3 = u .* p3;
    %第四个方程系数
    p4 = [V, -V];
    p4 = [K, K] * sqrt(eps0) .* p4;
    p4 = 1 ./ (nomega - (-N:N))' .* p4;
    p4 = u .* p4;
    %消去非齐次的部分
    p01 = (p2 / sqrt(eps1) + p1) / 2;
    p02 = (p4 / sqrt(eps2) - p3);
    p = [p01; p02];
    einc = zeros(Num, 1);
    einc(N + 1) = 1;
    b = [einc; zeros(Num, 1)];
    %解出系数向量
    A = p \ b;
    %代回求反射场
    eref = p1 * A - einc;
    %代回求透射场
    etra = p3 * A;
    %保存基频的幅度
    En0(loop, :) = [nomega, eref(N + 1), etra(N + 1)];
    En1(loop, :) = [nomega, eref(N + 2), etra(N + 2)];
    Enm1(loop, :) = [nomega, eref(N), etra(N)];
    En2(loop, :) = [nomega, eref(N + 3), etra(N + 3)];
    Enm2(loop, :) = [nomega, eref(N - 1), etra(N - 1)];
    loop = loop + 1;
end

%% 画图
h1 = figure;
hold on
plot(En0(:, 1), abs(En0(:, 2)), 'b.')
plot(En0(:, 1), abs(En0(:, 3)), 'r.')
legend('Reflection', 'Transmission')
title('The magnitude of reflection and transmission coefficient (n=0).')

h2 = figure;
hold on
plot(En1(:, 1), abs(En1(:, 2)), 'b.')
plot(En1(:, 1), abs(En1(:, 3)), 'r.')
legend('Reflection', 'Transmission')
title('The magnitude of reflection and transmission coefficient (n=1).')

h3 = figure;
hold on
plot(Enm1(:, 1), abs(Enm1(:, 2)), 'b.')
plot(Enm1(:, 1), abs(Enm1(:, 3)), 'r.')
legend('Reflection', 'Transmission')
title('The magnitude of reflection and transmission coefficient (n=-1).')

h4 = figure;
hold on
plot(En0(:, 1), angle(En0(:, 2)) * 180 / pi, 'b-')
plot(En0(:, 1), angle(En0(:, 3)) * 180 / pi, 'r-')
legend('Reflection', 'Transmission')
title('The phase of reflection and transmission coefficient (n=0).')

h5 = figure;
hold on
plot(En1(:, 1), angle(En1(:, 2)) * 180 / pi, 'b-')
plot(En1(:, 1), angle(En1(:, 3)) * 180 / pi, 'r-')
xlim([0 2])
legend('Reflection', 'Transmission')
title('The phase of reflection and transmission coefficient (n=1).')

h6 = figure;
hold on
plot(Enm1(:, 1), angle(Enm1(:, 2)) * 180 / pi, 'b-')
plot(Enm1(:, 1), angle(Enm1(:, 3)) * 180 / pi, 'r-')
xlim([0 2])
legend('Reflection', 'Transmission')
title('The phase of reflection and transmission coefficient (n=-1).')

h7 = figure;
hold on
plot(En2(:, 1), abs(En2(:, 2)), 'b.')
plot(En2(:, 1), abs(En2(:, 3)), 'r.')
legend('Reflection', 'Transmission')
title('The magnitude of reflection and transmission coefficient (n=2).')

h8 = figure;
hold on
plot(Enm2(:, 1), abs(Enm2(:, 2)), 'b.')
plot(Enm2(:, 1), abs(Enm2(:, 3)), 'r.')
legend('Reflection', 'Transmission')
title('The magnitude of reflection and transmission coefficient (n=-2).')

h9 = figure;
hold on
plot(En2(:, 1), angle(En2(:, 2)) * 180 / pi, 'b-')
plot(En2(:, 1), angle(En2(:, 3)) * 180 / pi, 'r-')
xlim([0 2])
legend('Reflection', 'Transmission')
title('The phase of reflection and transmission coefficient (n=2).')

h10 = figure;
hold on
plot(Enm2(:, 1), angle(Enm2(:, 2)) * 180 / pi, 'b-')
plot(Enm2(:, 1), angle(Enm2(:, 3)) * 180 / pi, 'r-')
xlim([0 2])
legend('Reflection', 'Transmission')
title('The phase of reflection and transmission coefficient (n=-2).')
