function find_missing_folders()

    % 1. 定义目标根目录和参数范围
    root_dir = 'D:\desktop\waiting';
    output_filename = fullfile(root_dir, 'missing_combinations.txt');
    
    % 参数 a 的取值范围
    a_values = [-5 -6 -7 -8 -10 -12 -14 -16 -18 -20 -22 -24 -26 -28 -30];
    
    % 参数 b 的取值范围
    b_values = [-32 -34 -35 -36 -37 -38 -40 -42 -44 -46 -48 -50 -52 -54 -56 -58 -60];
    
    % 检查目录是否存在
    if ~isfolder(root_dir)
        fprintf(2, '错误：指定目录不存在：%s\n', root_dir);
        return;
    end

    fprintf('--- 正在分析目录下的文件夹组合 ---\n');

    % 2. 生成所有期望的文件夹组合
    
    % 初始化一个 Cell 数组来存储所有期望的文件夹名 (例如: '-5_-32')
    expected_folders = {};
    
    % 使用两重循环生成所有可能的 a_b 组合
    for i = 1:length(a_values)
        a = a_values(i);
        for j = 1:length(b_values)
            b = b_values(j);
            % MATLAB 中使用 sprintf 来格式化字符串，生成 'a_b' 格式
            folder_name = sprintf('%d_%d', a, b);
            expected_folders{end+1} = folder_name;
        end
    end
    
    % 将期望的文件夹名转换为 Set 结构，便于快速查找（逻辑上）
    expected_set = containers.Map(expected_folders, num2cell(true(size(expected_folders))));


    % 3. 读取当前目录下已存在的文件夹
    
    % 查找 root_dir 下的所有子文件夹，排除 '.' 和 '..'
    dir_info = dir(root_dir);
    existing_dir_list = dir_info([dir_info.isdir] & ~strcmp({dir_info.name}, '.') & ~strcmp({dir_info.name}, '..'));
    
    % 提取已存在的文件夹名
    existing_folder_names = {existing_dir_list.name};
    
    % 4. 找出缺失的组合
    
    missing_folders = {};
    
    % 遍历所有期望的组合，检查它们是否已存在
    for i = 1:length(expected_folders)
        expected_name = expected_folders{i};
        
        % 检查该名称是否在实际存在的文件夹列表中
        % 使用 ismember 函数或容器可以高效地进行查找
        if ~ismember(expected_name, existing_folder_names)
             missing_folders{end+1} = expected_name;
        end
    end

    % 5. 将缺失的组合写入 TXT 文件
    
    if isempty(missing_folders)
        disp('所有参数组合的文件夹都已存在，没有缺失。');
        
        % 既然没有缺失，为避免产生空文件，可以删除上次生成的同名文件
        if isfile(output_filename)
             delete(output_filename);
             disp('已删除上次生成的空统计文件。');
        end
        return;
    end

    % 打开文件准备写入，'wt' 表示写入文本模式
    fid = fopen(output_filename, 'wt'); 
    
    if fid == -1
        fprintf(2, '错误：无法打开或创建输出文件：%s\n', output_filename);
        return;
    end
    
    fprintf('找到 %d 个缺失的文件夹组合，正在写入文件：%s\n', length(missing_folders), output_filename);
    
    % 遍历缺失的文件夹名
    for i = 1:length(missing_folders)
        folder_name = missing_folders{i};
        
        % 文件夹名格式是 a_b，我们需要提取 a 和 b
        % strsplit 函数可以按指定分隔符（这里是 '_'）分割字符串
        parts = strsplit(folder_name, '_');
        
        if length(parts) == 2
            a_part = parts{1};
            b_part = parts{2};
            
            % 将 a 和 b 部分以空格分隔写入文件
            fprintf(fid, '%s %s\n', a_part, b_part);
        end
    end
    
    % 关闭文件
    fclose(fid);
    
    fprintf('文件写入完成。您可以在以下路径查看结果：\n%s\n', output_filename);
end