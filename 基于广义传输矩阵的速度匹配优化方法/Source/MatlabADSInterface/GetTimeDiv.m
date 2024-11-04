function [index1, index2] = GetTimeDiv(waveform1, waveform2)
    %% 从两个相邻的波形找到时延
    [~, index1] = max(waveform1); % 找到 waveform1 的最大峰值及其索引
    [~, index2] = max(waveform2); % 找到 waveform2 的最大峰值及其索引
end