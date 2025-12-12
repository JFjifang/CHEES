function delete_incomplete_folders()

    % 1. 定义目标根目录
    root_dir = 'D:\desktop\waiting';

    % 检查目录是否存在
    if ~isfolder(root_dir)
        fprintf(2, '错误：指定目录不存在：%s\n', root_dir);
        return;
    end

    % 查找 root_dir 下的所有子文件夹
    % 'dir' 返回一个结构体数组，其中 'name' 是文件夹/文件名字
    % 'isdir' 字段指示是否是文件夹
    % '.' 和 '..' 是特殊目录，需要排除
    dir_info = dir(root_dir);
    
    % 过滤掉 '.' '..' 和不是文件夹的项目
    folder_list = dir_info([dir_info.isdir] & ~strcmp({dir_info.name}, '.') & ~strcmp({dir_info.name}, '..'));
    
    if isempty(folder_list)
        disp('在指定目录下没有找到任何子文件夹。');
        return;
    end

    % 初始化存储每个文件夹CSV文件数量的数组
    folder_names = {};
    csv_counts = [];

    fprintf('--- 正在统计每个文件夹中的CSV文件数量 ---\n');

    % 2. 统计每个文件夹里的CSV文件数量
    for i = 1:length(folder_list)
        folder_name = folder_list(i).name;
        full_path = fullfile(root_dir, folder_name);
        
        % 查找当前文件夹下的所有 .csv 文件
        % 注意：这里的查找是针对当前目录，不递归子目录
        csv_files = dir(fullfile(full_path, '*.csv'));
        count = length(csv_files);
        
        % 存储结果
        folder_names{end+1} = folder_name;
        csv_counts(end+1) = count;
        
        fprintf('文件夹：%s | CSV文件数量：%d\n', folder_name, count);
    end

    % 找出最大的CSV文件数量
    max_count = 0;
    if ~isempty(csv_counts)
        max_count = max(csv_counts);
    end

    fprintf('\n--- 统计结果 ---\n');
    fprintf('所有文件夹中，最多的CSV文件数量是：%d 个。\n', max_count);

    if max_count == 0
        disp('没有找到任何CSV文件，无法确定最大数量。脚本结束。');
        return;
    end

    % 找出文件数量不足 max_count 的文件夹
    % 逻辑：找出 csv_counts < max_count 的索引
    incomplete_indices = find(csv_counts < max_count);

    if isempty(incomplete_indices)
        disp('所有文件夹都达到了最大文件数，不需要删除任何文件夹。');
        return;
    end

    % 提取不符合要求的文件夹名字
    incomplete_folders = folder_names(incomplete_indices);
    
    fprintf('\n--- 不符合要求（文件数 < %d）的文件夹 ---\n', max_count);
    
    % 打印不符合要求的文件夹名字
    for i = 1:length(incomplete_folders)
        fprintf('* %s (文件数: %d)\n', incomplete_folders{i}, csv_counts(incomplete_indices(i)));
    end

    % 3. 删除这些不符合要求的文件夹
    % 警告提示
    user_input = input(sprintf('\n【警告】即将删除以上 %d 个文件夹及其所有内容。是否继续？(输入 YES 继续)：', length(incomplete_folders)), 's');

    if strcmpi(strtrim(user_input), 'YES')
        fprintf('\n--- 正在执行删除操作 ---\n');
        
        deleted_count = 0;
        for i = 1:length(incomplete_folders)
            folder_to_delete = fullfile(root_dir, incomplete_folders{i});
            
            % 'rmdir' 函数用于删除文件夹。's' 参数表示递归删除（包括子文件和子文件夹）
            [status, message, message_id] = rmdir(folder_to_delete, 's');
            
            if status
                fprintf('✅ 成功删除文件夹：%s\n', incomplete_folders{i});
                deleted_count = deleted_count + 1;
            else
                % 记录删除失败的情况
                fprintf(2, '❌ 删除失败：%s (错误信息: %s)\n', incomplete_folders{i}, message);
            end
        end
        
        fprintf('\n--- 结果总结 ---\n');
        fprintf('总共尝试删除 %d 个文件夹，成功删除 %d 个。\n', length(incomplete_folders), deleted_count);
        
    else
        disp('用户取消了删除操作。脚本已安全退出。');
    end
end