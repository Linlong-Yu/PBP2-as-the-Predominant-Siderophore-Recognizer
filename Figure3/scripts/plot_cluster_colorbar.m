function [] = plot_cluster(distm, group, labels, saveLoc)
%%% 自动调整legend和标签栏的布局
% 验证参数数量
if nargin < 4
    saveLoc = [];
end
if nargin < 3
    labels = [];
end

% 创建图形
figure()

% 绘制标签栏（如果提供了标签）
if ~isempty(labels)
    % 获取唯一标签并为每个标签分配颜色
    unique_labels = unique(labels, 'stable'); % 保持标签顺序
    num_labels = length(unique_labels);
    cmap = hsv(num_labels); % 生成一个颜色映射

    % 为每个标签创建颜色索引
    label_colors = zeros(length(labels), 3);
    for i = 1:num_labels
        label_colors(strcmp(labels, unique_labels{i}), :) = repmat(cmap(i, :), sum(strcmp(labels, unique_labels{i})), 1);
    end

    % 绘制标签栏
    subplot('Position', [0.12, 0.2, 0.03, 0.6]) % 调整标签栏的宽度和位置
    label_bar = reshape(label_colors, [length(labels), 1, 3]);
    imagesc(label_bar)
    axis equal
    axis off

    % 添加图例
    subplot('Position', [0.9, 0.2, 0.1, 0.6]) % 将图例放到右边
    hold on
    for i = 1:num_labels
        plot(0, 0, 's', 'MarkerSize', 10, 'MarkerEdgeColor', cmap(i, :), 'MarkerFaceColor', cmap(i, :));
    end
    hold off
    lgd = legend(unique_labels, 'Location', 'westoutside', 'Orientation', 'vertical', 'FontSize', 8); % 调整FontSize和位置
    lgd.Box = 'off'; % 去掉图例框
    axis off
end

% 绘制热图
subplot('Position', [0.15, 0.2, 0.65, 0.6]) % 调整热图的位置
imagesc(distm)
colorbar
caxis([0, 1])
axis equal
axis off
hold on

% 画出表示组的方框
if ~isempty(group)
    split_group = zeros(1, length(group));
    for i = 1:length(group)-1
        if group(i+1) == group(i)+1
            split_group(i) = 1;
        end
    end

    pos = [0, find(split_group == 1) length(group)];
    lens = [];
    for i = 2:length(pos)
        lens = [lens, pos(i)-pos(i-1)];
    end

    for i = 1:length(pos)-1
        rectangle('Position', [pos(i)+0.5 pos(i)+0.5 lens(i) lens(i)], 'EdgeColor', 'r', 'LineWidth', 0.4)
        hold on
    end
end
hold off

% 如果提供了保存位置，则保存图形为PDF
if ~isempty(saveLoc)
    % 确保路径以不带扩展名的字符串形式传入
    [filepath, name, ~] = fileparts(saveLoc);
    full_path = fullfile(filepath, name);
    print(gcf, full_path, '-dpdf', '-painters') % 使用矢量格式 PDF 输出
end

end
