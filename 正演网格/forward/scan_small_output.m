%% 扫描 CHEES 结果文件，找出 <10 KB 的 csv，并记录 a b c d
clear; clc;

% 1. 根目录路径（根据你的实际路径修改）
rootDir = 'D:\desktop\relargement_seed';

% 结果输出 txt
outFile = fullfile(rootDir, 'small_CHEES_files_a_b_c_d.txt');

% 小文件阈值：10 KB
sizeThreshold = 10 * 1024;  % 单位：字节

% 用来保存所有 (a, b, c, d)
results = [];  % 每行: [a b c d]

% 2. 列出根目录下的所有项
items = dir(rootDir);

for i = 1:numel(items)
    % 跳过非文件夹和 . / ..
    if ~items(i).isdir || startsWith(items(i).name, '.')
        continue;
    end

    folderName = items(i).name;
    folderPath = fullfile(rootDir, folderName);

    % ========= 解析文件夹名得到 a, b =========
    % 假设文件夹名形式为 "-5_-36"，即 a_b，允许负号
    tFolder = regexp(folderName, '^(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)$', 'tokens', 'once');
    if isempty(tFolder)
        % 不符合 a_b 形式的文件夹就跳过
        fprintf('跳过文件夹（命名不符合 a_b 格式）：%s\n', folderName);
        continue;
    end
    a = str2double(tFolder{1});
    b = str2double(tFolder{2});

    % 3. 找到该文件夹下所有 csv 文件
    csvFiles = dir(fullfile(folderPath, '*.csv'));

    for j = 1:numel(csvFiles)
        fileInfo = csvFiles(j);
        filePath = fullfile(folderPath, fileInfo.name);

        % 判断文件大小是否 < 10 KB
        if fileInfo.bytes < sizeThreshold

            % ========= 解析文件名得到 c, d =========
            % 假设文件名形式为 "CHEES_-22_4200.csv"
            % 即 CHEES_c_d.csv，允许负号
            fName = fileInfo.name;
            tFile = regexp(fName, '^CHEES_(-?\d+)_(-?\d+)\.csv$', 'tokens', 'once');

            if isempty(tFile)
                fprintf('文件名不符合 CHEES_c_d.csv 格式，已跳过：%s\n', filePath);
                continue;
            end

            c = str2double(tFile{1});
            d = str2double(tFile{2});

            % 记录一行结果
            results = [results; a, b, c, d];
        end
    end
end

% 4. 把结果写入 txt 文件
fid = fopen(outFile, 'w');
if fid == -1
    error('无法打开输出文件：%s', outFile);
end

% 写表头（如果不需要可以注释掉这一行）
fprintf(fid, 'a\tb\tc\td\n');

for k = 1:size(results, 1)
    fprintf(fid, '%g\t%g\t%g\t%g\n', ...
        results(k, 1), results(k, 2), results(k, 3), results(k, 4));
end

fclose(fid);

fprintf('完成！共找到 %d 个小于 10 KB 的 csv 文件。\n', size(results, 1));
fprintf('结果已写入：%s\n', outFile);
