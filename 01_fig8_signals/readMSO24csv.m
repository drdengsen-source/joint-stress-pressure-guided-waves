function [t, V, names] = readMSO24csv(file)
%READMSO24CSV  读 Tektronix MSO2x 系列保存的 CSV 波形（多通道）
%   [t, V, names] = readMSO24csv(file)
%   t：时间列（s）；V：各通道电压矩阵（每列一通道）；names：通道名（如 {'CH1','CH2','CH3'}）。
%   格式：前面若干行 "键,值" 元数据，然后一行以 TIME 开头的列名，再是数据。找不到 TIME 行时退化为
%   readmatrix（第一列当时间）。
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
hit = ~cellfun(@isempty, regexp(lines, '^\s*TIME\s*,', 'once', 'ignorecase'));
i0 = find(hit, 1);
if isempty(i0)
    M = readmatrix(file);
    t = M(:,1);  V = M(:,2:end);  names = arrayfun(@(k) sprintf('CH%d', k), 1:size(V,2), 'UniformOutput', false);
    return
end
names = strtrim(strsplit(lines{i0}, ','));
names = names(2:end);
body = lines(i0+1:end);
body = body(~cellfun(@isempty, strtrim(body)));
n = numel(names) + 1;
C = textscan(strjoin(body, newline), repmat('%f', 1, n), 'Delimiter', ',', 'CollectOutput', true);
M = C{1};
t = M(:,1);  V = M(:,2:end);
end
