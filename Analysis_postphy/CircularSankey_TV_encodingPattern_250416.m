clear; clc; close all;
load('D:\DATA\listTV3_PPC_RT.mat')
target = 'PPC';

data = {};
data{1} = lick.';
data{2} = five.';
data{3} = ten.';
data{4} = hit.';
data{5} = fa.';
data{6} = cr.';
data{7} = []; %miss.';
data{8} = rew.';





clear feature_martix
all_num = sortrows(unique(vertcat(data{:})));
for i = 1:numel(all_num)
    neuron = all_num(i);
    for j = 1:length(data)
        feature_matrix(i,j) = ismember(neuron, data{1,j});
    end
end
clearvars -except feature_matrix target


%%
tv_label = {'lick', 'five', 'ten', 'Hit', 'FA', 'CR', 'miss','Reward'};
num_tv = size(tv_label,2);

clear idx
idx = sum(feature_matrix, 2) >= 2;
multi = feature_matrix(idx, 1:num_tv);

for i = 1:num_tv
    for j = 1:num_tv
        if i ~= j
            multi2(i,j) = sum(multi(:,i) & multi(:,j));
        end
    end
end


[combo_unique, ~, combo_idx] = unique(multi, 'rows');
combo_type = accumarray(combo_idx, 1);

curr = 1;
for i = 1:size(combo_unique,1)
    encodeTV = find(combo_unique(i,:));
    if length(encodeTV) < 2
        continue
    end
    combo_str = strjoin(tv_label(encodeTV), '+');
    temp_counts(curr) = struct('combo', combo_str,'count', combo_type(i),'TV', encodeTV);
    curr= curr +1;
end

[~, idx] = sort([temp_counts.count], 'descend');
combo_counts = temp_counts(idx);

    
    
    
    RadiusSankeyPlot(multi2, tv_label, target, combo_counts);

%% radius sankey plot
function RadiusSankeyPlot(data_matrix, tv_label,target_area, top_combos)
   
figure('Name', ['Multiplex encoding: ', target_area], 'Position', [300 300 600 600])

    num_tv = length(tv_label);
    theta = linspace(0, 2*pi, num_tv+1);
    theta = theta(1:end-1);
    x = cos(theta);
    y = sin(theta);
    max_width = 10;
    max_count = max(data_matrix(:));
    dist = 0.5;
    tt = linspace(0, 1, 100);
    colors = flipud(bone(max(10, max_count)));

    hold on
    for i = 1:num_tv
        plot(x(i), y(i), 'o', 'MarkerSize', 20, 'MarkerFaceColor', 'k')
        text(1.25*x(i), 1.25*y(i), tv_label{i}, 'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    
    for i = 1:num_tv
        for j = i+1:num_tv
            if data_matrix(i,j) > 0
                clear width tv1 tv2 mid_angle xx yy xxx yyy
                mid_angle = (theta(i) + theta(j))/2;

                if abs(theta(i) - theta(j)) > pi
                    mid_angle = mid_angle + pi;
                end

                xx = dist * cos(mid_angle);
                yy = dist * sin(mid_angle);
                xxx = ((1-tt).^2*x(i)) + (2*(1-tt).*tt*xx) + (tt.^2*x(j));
                yyy = ((1-tt).^2*y(i)) + (2*(1-tt).*tt*yy) + (tt.^2*y(j));
                width = max_width * data_matrix(i,j) / max_count;

                plot(xxx, yyy, 'Color', [0.8 0.8 0.8], 'LineWidth', width)
            end
        end
    end
    
    for rank = 1:min(10, length(top_combos))
        temp_tv = top_combos(rank).TV;

        for i = 1:length(temp_tv)
            for j = i+1:length(temp_tv)
                clear width tv1 tv2 mid_angle xx yy xxx yyy
                tv1 = temp_tv(i);
                tv2 = temp_tv(j);
                mid_angle = (theta(tv1) + theta(tv2))/2;

                if abs(theta(tv1) - theta(tv2)) > pi
                    mid_angle = mid_angle + pi;
                end

                xx = dist * cos(mid_angle);
                yy = dist * sin(mid_angle);
                xxx = ((1-tt).^2*x(tv1)) + (2*(1-tt).*tt*xx) + (tt.^2*x(tv2));
                yyy = ((1-tt).^2*y(tv1)) + (2*(1-tt).*tt*yy) + (tt.^2*y(tv2));
                width = max_width * data_matrix(tv1,tv2) / max_count;

                plot(xxx, yyy, 'Color', colors(data_matrix(tv1,tv2),:), 'LineWidth', width)
%                 plot(xxx, yyy, 'Color', [.5,.5,.5], 'LineWidth', width)
                if data_matrix(tv1,tv2) > max_count/10
                    mid_idx = round(length(tt)/2);
                    text(xxx(mid_idx), yyy(mid_idx), num2str(data_matrix(tv1,tv2)),'FontSize', 10, 'HorizontalAlignment', 'center');
                end
            end
        end
    end
    
    axis off
    title([target_area, ': Task Variable Co-encoding'], 'FontSize', 15);
end
