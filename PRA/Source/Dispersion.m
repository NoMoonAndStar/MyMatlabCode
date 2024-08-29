%% 解时间光子晶体的色散曲线
clear
clc
close all

%% 定义参数
eps0 = 5.25;
Deltaeps = 0.85;
norepsm1 = 1i * 0.5 * Deltaeps / eps0;
noreps1 = -1i * 0.5 * Deltaeps / eps0;
knum = 1001;
N = 3;
solutions = NaN(knum, N + 1);

%% 扫k解omega 解方程太慢了
% syms noromega
% assume(noromega < 0.5 & noromega > 0)
% loop = 1;
%
% tic
%
% for nork = linspace(0, 3, knum)
%     D = sym(zeros(N));
%
%     for i = 1:N
%
%         for j = 1:N
%
%             if (i - j) == 0
%                 D(i, j) = (noromega - (i - (N + 1) / 2)) ^ 2 - nork ^ 2;
%             elseif (i - j) == 1
%                 D(i, j) = noreps1 * (noromega - (i - (N + 1) / 2)) ^ 2;
%             elseif (i - j) == -1
%                 D(i, j) = norepsm1 * (noromega - (i - (N + 1) / 2)) ^ 2;
%             end
%
%         end
%
%     end
%
%     tempsolve = solve(det(D) == 0, 'Real', true);
%
%     if ~isempty(tempsolve)
%         tempsolve = double(tempsolve);
%
%         if isreal(tempsolve)
%             solutions(loop, 1:length(tempsolve) + 1) = [nork, tempsolve];
%         end
%
%     end
%
%     loop = loop + 1;
% end
%
% toc

%% 扫omega解k
tic
loop = 1;

for nomega = linspace(0, 1, knum)
    D = zeros(N);

    for i = 1:N

        for j = 1:N

            if i - j == 0
                D(i, j) = (nomega - (i - (N + 1) / 2)) ^ 2;
            elseif i - j == 1
                D(i, j) = noreps1 * (nomega - (i - (N + 1) / 2)) ^ 2;
            elseif i - j == -1
                D(i, j) = norepsm1 * (nomega - (i - (N + 1) / 2)) ^ 2;
            end

        end

    end

    K = eig(D);

    if ~isempty(K)
        nork = sqrt(K);
        solutions(loop, 1:length(K) + 1) = [nomega, nork'];
    end

    loop = loop + 1;
end

toc

%% 绘图
h1 = figure;
hold on
% plot(solutions(:, 1), solutions(:, 2), 'b.')
% plot(solutions(:, 1), 1 - solutions(:, 2), 'r.')
for i = 1:N
    plot(solutions(:, i + 1), solutions(:, 1), 'b.')
end

xlim([0 3])
xlabel('Normalized k')
ylabel('Normalized omega')
title(['Dispersion relation of the time-varying crystal with Deltaeps=', num2str(Deltaeps)])
