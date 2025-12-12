clc; clear; close all;

%% 1. 设置路径参数
% 源文件所在的目录 (存放 .csv 文件)
sourceDir = 'D:\desktop\repairing1';

% 目标父目录 (存放各种 a_b 文件夹)
destParentDir = 'D:\desktop\enlargement';

%% 2. 获取源目录下符合模式的文件
% 查找所有以 CHEES_ 开头的 csv 文件
filePattern = fullfile(sourceDir, 'CHEES_*.csv');
csvFiles = dir(filePattern);

% 初始化计数器
movedCount = 0;
skippedCount = 0;

fprintf('开始处理...\n');

%% 3. 遍历文件并移动
for k = 1:length(csvFiles)
    baseFileName = csvFiles(k).name;
    fullSourcePath = fullfile(sourceDir, baseFileName);
    
    % 解析文件名
    % 文件名格式示例: CHEES_-5_-48_-22_4200.csv
    % 我们使用下划线 '_' 进行分割
    parts = strsplit(baseFileName, '_');
    
    % 检查文件名分割后是否有足够的长度，防止文件名格式不对报错
    if length(parts) >= 3
        % 根据格式:
        % parts{1} 是 'CHEES'
        % parts{2} 是 a (例如 '-5')
        % parts{3} 是 b (例如 '-48')
        a = parts{2};
        b = parts{3};
        
        % 构建目标文件夹名称: a_b (例如 '-5_-48')
        targetFolderName = [a, '_', b];
        
        % 构建目标文件夹的完整路径
        targetFolderPath = fullfile(destParentDir, targetFolderName);
        
        % 检查目标文件夹是否存在
        if isfolder(targetFolderPath)
            % 构建目标文件的完整路径
            fullDestPath = fullfile(targetFolderPath, baseFileName);
            
            % 移动文件
            try
                movefile(fullSourcePath, fullDestPath);
                fprintf('已移动: %s  -->  %s\n', baseFileName, targetFolderName);
                movedCount = movedCount + 1;
            catch ME
                fprintf('移动失败: %s (错误: %s)\n', baseFileName, ME.message);
                skippedCount = skippedCount + 1;
            end
        else
            % 如果目标文件夹不存在，打印警告并跳过
            fprintf('跳过 (目标文件夹不存在): %s (目标: %s)\n', baseFileName, targetFolderName);
            skippedCount = skippedCount + 1;
        end
    else
        fprintf('跳过 (文件名格式不符): %s\n', baseFileName);
        skippedCount = skippedCount + 1;
    end
end

%% 4. 结束报告
fprintf('------------------------------------------------\n');
fprintf('处理完成。\n');
fprintf('成功移动文件数: %d\n', movedCount);
fprintf('跳过/失败文件数: %d\n', skippedCount);