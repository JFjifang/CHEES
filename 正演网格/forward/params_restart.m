%% 收集 D:\desktop\-20_-40 下所有 a_b 形式子文件夹的参数，写入 params_restart.txt
clear; clc;

% 1. 根目录
rootDir = 'D:\desktop\-20_-40';

% 2. 输出文件
outFile = fullfile(rootDir, 'params_restart.txt');

% 保存所有 [a b]
params = [];   % 每一行：[a b]

% 3. 列出根目录下所有项目
items = dir(rootDir);

for i = 1:numel(items)
    % 只要子文件夹，且排除 '.' '..'
    if ~items(i).isdir || startsWith(items(i).name, '.')
        continue;
    end

    folderName = items(i).name;

    % 4. 解析文件夹名，假设是 "-5_-36" 这样的 a_b 形式，允许负号
    tokens = regexp(folderName, '^(-?\d+)_(-?\d+)$', 'tokens', 'once');
    if isempty(tokens)
        % 如果不是 a_b 格式，可以选择打印提示，也可以直接跳过
        fprintf('跳过命名不符合 a_b 格式的文件夹：%s\n', folderName);
        continue;
    end

    a = str2double(tokens{1});
    b = str2double(tokens{2});

    params = [params; a, b];
end

% 5. 写入 params_restart.txt
fid = fopen(outFile, 'w');
if fid == -1
    error('无法打开输出文件：%s', outFile);
end

% 如果不需要表头，就不要写这一行
% fprintf(fid, 'a\tb\n');

for k = 1:size(params, 1)
    fprintf(fid, '%g\t%g\n', params(k, 1), params(k, 2));  % 每行: a<tab>b
end

fclose(fid);

fprintf('完成！共记录 %d 个文件夹的参数。\n', size(params, 1));
fprintf('结果已写入：%s\n', outFile);
