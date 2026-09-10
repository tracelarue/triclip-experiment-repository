%% Spot-check processed valve annotations
clear
close all
clc

[file, folder] = uigetfile('*.mat', 'Select processed annotation file');
if isequal(file, 0)
    return
end

S = load(fullfile(folder, file));

% Find the annotation structure regardless of its variable name.
names = fieldnames(S);
A = [];
for k = 1:numel(names)
    candidate = S.(names{k});
    if istable(candidate)
        candidate = table2struct(candidate);
    end
    if isstruct(candidate) && isfield(candidate, 'AnnulusOutline')
        A = candidate;
        break
    end
end

assert(~isempty(A), ...
    'Select a processed annotation file containing AnnulusOutline.');

for k = 1:numel(A)
    a = A(k);
    V = a.AnnulusOutline.Vertices;

    figure('Name', sprintf('Annotation %d: %s', k, char(a.Test)), ...
        'Color', 'w', 'Position', [100 100 1100 800]);
    ax = axes;
    hold(ax, 'on');

    % Background image and geometry use the same aligned pixel coordinates.
    if ~isempty(a.Frame)
        imshow(a.Frame, a.Ref, 'Parent', ax);
        hold(ax, 'on');
    else
        set(ax, 'YDir', 'reverse');
    end

    plot(ax, a.AnnulusOutline, ...
        'FaceColor', 'none', 'EdgeColor', 'y', 'LineWidth', 2);

    % Commissures
    idx = a.Comm_Idx(:);
    labels = {'AP', 'AS', 'SP'};
    plot(ax, V(idx,1), V(idx,2), 'go', ...
        'MarkerFaceColor', 'g', 'MarkerSize', 7);
    for j = 1:numel(idx)
        text(ax, V(idx(j),1), V(idx(j),2), ['  ' labels{j}], ...
            'Color', 'g', 'FontWeight', 'bold', 'FontSize', 12);
    end

    % Pins
    idx = a.Pin_Idx(:);
    plot(ax, V(idx,1), V(idx,2), 'cs', ...
        'MarkerFaceColor', 'c', 'MarkerSize', 7);
    for j = 1:numel(idx)
        text(ax, V(idx(j),1), V(idx(j),2), sprintf('  P%d', j), ...
            'Color', 'c', 'FontWeight', 'bold', 'FontSize', 12);
    end

    % AP and SL measurement lines
    idx = a.Comm_Idx([1 2]);
    plot(ax, V(idx,1), V(idx,2), 'r--', 'LineWidth', 1.5);
    idx = a.SL_Idx;
    plot(ax, V(idx,1), V(idx,2), 'b--', 'LineWidth', 1.5);

    % Clip axes and centers
    axisFields   = {'ClipAxis_Idx', 'ClipAxis2_Idx'};
    centerFields = {'ClipCenter', 'ClipCenter2'};
    clipColors   = [1 0 1; 1 0.5 0];

    for c = 1:2
        color = clipColors(c,:);
        if isfield(a, axisFields{c}) && ~isempty(a.(axisFields{c}))
            idx = a.(axisFields{c});
            plot(ax, V(idx,1), V(idx,2), 'o--', ...
                'Color', color, 'LineWidth', 2, 'MarkerSize', 6);
        end

        if isfield(a, centerFields{c}) && ~isempty(a.(centerFields{c}))
            p = a.(centerFields{c});
            plot(ax, p(1), p(2), 'x', ...
                'Color', color, 'LineWidth', 3, 'MarkerSize', 12);
            text(ax, p(1), p(2), sprintf('  Clip %d', c), ...
                'Color', color, 'FontWeight', 'bold', 'FontSize', 12);
        end
    end

    plot(ax, 0, 0, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
    axis(ax, 'image');

    title(ax, {
        sprintf('%s | %s | Annotation %d/%d', ...
            char(a.Valve), char(a.Test), k, numel(A))
        sprintf(['Area: %.1f mm^2 | Perimeter: %.1f mm | ' ...
                 'AP: %.1f mm | SL: %.1f mm | Gap: %.1f mm^2'], ...
            a.Area, a.Perimeter, a.AP_Diameter, ...
            a.SL_Diameter, a.Gap_Area)
        'Yellow: outline | Green: commissures | Cyan: pins | Red: AP | Blue: SL'
        }, 'Interpreter', 'none');

    zoom(gcf, 'on');
end