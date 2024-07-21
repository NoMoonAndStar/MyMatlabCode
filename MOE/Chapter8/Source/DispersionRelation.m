d = 0.01;
syms betad
kd = 0:0.001:1;
k = kd / d;
sd = zeros(length(kd), 2);

for i = 1:length(kd)
    f = cos(betad) - cos(kd(i)) + 2 * kd(i) * sin(kd(i));
    s = solve(f == 0);
    sd(i, :) = double(s);
end

beta = sd / d;
[rows, cols] = find(imag(beta) == 0);
realbeta = real(beta(rows, cols));

plot(k, realbeta', '.');
