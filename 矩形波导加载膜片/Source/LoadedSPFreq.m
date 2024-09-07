S = load('..\Data\加载波导S21相位频域仿真.txt');
Freq0 = S(:, 1);
SPhase0 = S(:, 2);
index = find(Freq0 == 2.5);
Freq = Freq0(index:end);
SPhase = SPhase0(index:end);
uS = -unwrap(SPhase, 180);
uS = uS - uS(1);
L = 400 * 1e-3;
Beta = uS * pi / 180 / L;
Beta0 = Beta;

TheoryData = load('..\Data\Disepersion.mat');

index1 = find(Freq == 9.60);
index2 = find(Freq == 15.2);
index3 = find(Freq == 20);
Beta(index1 + 1:index2 - 1) = NaN;
Beta(index2:end) = Beta(index2:end) + 315.536 - Beta(index2);
Beta(index3:end) = NaN;

h1 = figure;
hold on
plot(Beta, Freq, 'b--', 'LineWidth', 3)
plot(TheoryData.Beta(:, 1), TheoryData.k, 'r.', 'MarkerSize', 5)
plot(TheoryData.Beta(:, 2), TheoryData.k, 'r.', 'MarkerSize', 5)
xlabel('Beta(rad/m)')
ylabel('Frequency(GHz)')
ylim([0 20])
grid on
legend('根据相位展开的色散曲线', '理论色散曲线')
title('理论色散曲线和S21相位展开色散曲线对比(频域仿真)')

h2 = figure;
hold on
plot(Freq0, SPhase0, 'b-')
xlabel('Frequency(GHz)')
ylabel('Phase(deg)')
grid on
legend('S21相位')

h3 = figure;
hold on
plot(Beta0, Freq, 'b-')
grid on
ylabel('Frequency(GHz)')
xlabel('Beta(rad/m)')
legend('根据相位展开的色散曲线')
