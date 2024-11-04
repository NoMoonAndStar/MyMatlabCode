figure;
grid on
a = input('暂停')
for i = 1:97
    w = plot(trwave(:,i), 'k-', 'LineWidth', 1.5);
    ylim([-2e-6, 2e-6])
    pause(0.1)
end
