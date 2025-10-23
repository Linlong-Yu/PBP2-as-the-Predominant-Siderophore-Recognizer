function [] = plot_cluster(distm, group, saveLoc)

% 如果只输入两个参数，不保存
if nargin < 3
    saveLoc = '';
end

% 绘图
figure()
imagesc(distm)
colorbar
caxis([0, 1])
axis equal
axis off
hold on

% draw squares which represent groups
if ~isempty(group)
    split_group = zeros(1, length(group));
    for i = 1:length(group)-1
        if group(i+1) == group(i)+1
            split_group(i) = 1;
        end
    end
    
    pos = [0, find(split_group == 1), length(group)];
    lens = diff(pos);
    
    for i = 1:length(lens)
        rectangle('Position', [pos(i)+0.5, pos(i)+0.5, lens(i), lens(i)], ...
                  'EdgeColor', 'r', 'LineWidth', 0.2)
    end
end
hold off

% 保存为 PDF（如果提供了路径）
if ~isempty(saveLoc)
    [filepath, name, ~] = fileparts(saveLoc);
    full_path = fullfile(filepath, name);
    print(gcf, full_path, '-dpdf', '-vector');
end

end
