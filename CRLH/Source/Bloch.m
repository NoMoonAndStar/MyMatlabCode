%%计算空间周期结构的特性

%初始化
clear
clc
close all

%扫频范围
Omega = 2 * pi * 5e10;
omega = 0:Omega / 10000:Omega;

%矩形波导参量
a = 22.86 * 1e-3; %矩形波导的长度
b = 5.08 * 1e-3; %矩形波导的宽度
h = 3.81 * 1e-3; %膜片高度
p = 5 * 1e-3; %膜片间距

c = 3e8; %真空光速
k0 = omega / c; %真空波数

Beta10 = sqrt(k0 .^ 2 - (pi / a) ^ 2); %TE10模相位常数

B = 4 * b * Beta10 * log(sec((pi * h) / (2 * b))) / pi; %膜片的归一化导纳

Beta = zeros(length(omega), 4); %周期结构的相位常数

for i = 1:length(omega)

    if imag(Beta10(i)) ~= 0
        continue
    else
        solution0 = acos(cos(Beta10(i) * p) - B(i) * sin(Beta10(i) * p) / 2);

        if isreal(solution0)
            solutions = [];

            for j = 0:3 % 因为周期是2π，所以4π内有2个完整的周期
                solutions = [solutions solution0 + 2 * pi * j]; %#ok<AGROW>
                solutions = [solutions -solution0 + 2 * pi * (j + 1)]; %#ok<AGROW>
            end

            ValidSolutions = solutions(-4 * pi <= solutions & solutions <= 4 * pi);
            Beta(i, 1:length(ValidSolutions)) = ValidSolutions;
        end

    end

end

%%fig3.33(b)
k = omega' * p / c;
idx = Beta(:, 1) ~= 0;
Beta = reshape(Beta(Beta(:) ~= 0), [], 4);
k = k(idx);

if length(k) < length(Beta)
    Delta = length(Beta) - length(k);
    Beta(end - Delta + 1:end, :) = [];
end

h1 = figure;
hold on
plot([0 0], [0 k(end)], 'LineWidth', 1, 'Color', 'k')
plot(Beta, k, 'b.')
plot(-Beta, k, 'b.')
xlim([-4 * pi, 4 * pi])
ylim([0, k(end)])
xlabel('Beta*p')
ylabel('k*omega/p')
title('Dispersion relation of the periodic structure')
plot([pi pi], [0 k(end)], '--', 'Color', 'r')
plot([-pi -pi], [0 k(end)], '--', 'Color', 'r')
plot([-4 * pi, 4 * pi], [2, 2], '--', 'Color', 'r')
