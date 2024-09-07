S = load('..\Data\空波导S21相位.txt');
Freq = S(:, 1);
SPhase = S(:, 2);
index = find(Freq == 2.5);
Freq = Freq(index:end);
SPhase = SPhase(index:end);
uS = -unwrap(SPhase, 180);
uS = uS - uS(1);
L = 990 * 1e-3;
Beta = uS * pi / 180 / L;

a = 60 * 1e-3;
f = linspace(2.5e9, 5e9, 1000);
k = f * 2 * pi / 3e8;
kc = 2 * pi / 2 / a;
beta = sqrt(k .^ 2 - kc .^ 2);

h1 = figure;
hold on
plot(Beta, Freq, 'b-', 'LineWidth', 2)
plot(beta, f / 1e9, 'r-', 'LineWidth', 2)
xlim([0 100])
ylim([2.5 5])
xlabel('Beta(rad/m)')
ylabel('Frequency(GHz)')
grid on
legend('根据相位展开的色散曲线', '理论色散曲线')
