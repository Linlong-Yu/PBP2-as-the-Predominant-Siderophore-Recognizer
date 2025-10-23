function plot_heatmap(distm, group, syn_type, save_path)
    
    if size(distm,1) ~= size(distm,2)
        error('distm 必须是 n×n 矩阵');
    end
    n = size(distm,1);
    if length(group) ~= n
        error('group 长度必须等于 distm 尺寸');
    end
    if length(syn_type) ~= n
        error('syn_type 长度必须等于 distm 尺寸');
    end

    [unique_types, ~, type_ids] = unique(syn_type, 'stable');
    cmap_types = lines(length(unique_types));
   
    for i = 1:length(unique_types)
        if strcmpi(unique_types{i}, 'NIS')
            cmap_types(i, :) = [231/255, 36/255, 60/255];  % #9c2232
        elseif strcmpi(unique_types{i}, 'NRPS')
            cmap_types(i, :) = [63/255, 185/255, 74/255];    % #014948
        end
    end
    
    type_colors = cmap_types(type_ids, :);

  
    figure('Position', [100 100 1200 800], 'Color', 'w');
    
    ax1 = subplot('Position', [0.05 0.1 0.03 0.8]);
    image(permute(type_colors, [1 3 2]));
    axis off;
    set(ax1, 'YDir', 'reverse');

    ax2 = subplot('Position', [0.10 0.1 0.65 0.8]);
    imagesc(distm, [0 1]);
    colormap(ax2, parula);
    
    cb = colorbar('Position', [0.77 0.1 0.02 0.8]);
    
    axis square;
    set(ax2, 'XTick', [], 'YTick', []);

    hold on;
   
    split_group = zeros(1, length(group));
    for i = 1:length(group)-1
        if group(i+1) == group(i)+1
            split_group(i) = 1;
        end
    end
    pos = [0, find(split_group == 1), length(group)];
    lens = diff(pos);
    for i = 1:length(lens)
        rectangle('Position', [pos(i)+0.51, pos(i)+0.51, lens(i)-0.02, lens(i)-0.02], ...
                  'EdgeColor', 'r', 'LineWidth', 1.5);
    end
    hold off;

    ax_legend = axes('Position', [0.81 0.65 0.18 0.25]); 
    
    xlim(ax_legend, [0 1]);
    ylim(ax_legend, [0 1]);
    axis(ax_legend, 'off');
    hold(ax_legend, 'on');
    box_h = 0.08;  
    box_w = 0.2;   
    spacing = 0.03;
    text_size = 8;
    total_height = length(unique_types) * (box_h + spacing) - spacing;
    start_y = 0.9;

    for i = 1:length(unique_types)
        y_pos = start_y - (i-1) * (box_h + spacing);
        rectangle('Position', [0.05 y_pos box_w box_h], ...
            'FaceColor', cmap_types(i,:), 'EdgeColor', 'k', 'LineWidth', 0.5, ...
            'Parent', ax_legend);
        text(0.3, y_pos + box_h/2, unique_types{i}, ...
            'Parent', ax_legend, 'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', text_size, 'Interpreter', 'none', ...
            'Color', 'k');
    end
    hold(ax_legend, 'off');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [12 8]);
    set(gcf, 'PaperPosition', [0 0 12 8]);
    set(gcf, 'PaperPositionMode', 'manual');
    
    print(gcf, save_path, '-dpdf', '-vector', '-r600');

    fprintf('已保存 PDF: %s\n', save_path);
end
