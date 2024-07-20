d = 0.01;
syms betad
kd = 0:1:10;
k = kd / d;
sd = zeros(length(kd), 2);

for i = 1:length(kd)
    f = cos(betad) - cos(kd(i)) + 2 * kd(i) * sin(kd(i));
    s = solve(f == 0);
    sd(i, :) = double(s);
end

beta = sd / d;
realbeta = real(beta);
plot(k, real(beta)');
