%% load_valve_annotations.m
% Asks user for a directory and loads in all .mat annotations
% Uniformly orients orientations
% Saves out as TestNum_Annotation.mat

% Define:
% AP Axis: AS -> AP commissure
% SL Axis: SP -> Anterior free wall, through centroid and perpendicular to
%                AP axis

clc
clear all
close all

% Get the full path of this script
scriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(scriptPath);
projectRoot = fileparts(scriptDir);

% Build data paths
dataDir = fullfile(projectRoot, 'data');

% ask user to select folder 
% folderPath = uigetdir(dataDir);
folderPath = uigetdir('C:\Users\ch49693\Box\Research\Experiments\Clipping\');

% Create pattern to identify all .mat files
filePattern = fullfile(folderPath, '*.mat');
% Get information about all matching .mat files
matFiles = dir(filePattern);

% Initialize cell arrays to store data and filenames
annotationData = cell(length(matFiles), 1);
traceData = annotationData;
names = cell(length(matFiles), 1);

%% Loop through each .mat file
for k = 1:length(annotationData)
    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folderPath, baseFileName);
    fprintf('Now reading %s\n', fullFileName);
    
    % Read .mat file into table
    data = load(fullFileName,"annotation");
    names{k} = baseFileName;
    annotationData{k} = data;

    % crop and convert everything to CROPPED coords
    crop = data.annotation.cropRect_xywh;   % [x y w h] in ORIGINAL coords
    x = crop(1); y = crop(2); w = crop(3); h = crop(4);
    frame = data.annotation.frame(y:y+h-1, x:x+w-1, :);

    offset = [x-1, y-1];          % original -> cropped: pC = p0 - off
    outline = data.annotation.outline_xy - offset;
    comm    = data.annotation.commissures_xy    - offset;
    pins    = data.annotation.pins_xy    - offset;

    % resample outline
    % ensure closed polygon
    if any(outline(1,:) ~= outline(end,:))
        outline(end+1,:) = outline(1,:);
    end
    
    % Compute cumulative arc length
    d = sqrt(sum(diff(outline).^2,2));
    s = [0; cumsum(d)];
    L = s(end);
    
    % Target parameter locations
    N = 1000;
    s_new = linspace(0, L, N+1)';
    s_new(end) = [];   % avoid duplicating start/end
    
    % Interpolate x and y separately
    x_new = interp1(s, outline(:,1), s_new, 'spline');
    y_new = interp1(s, outline(:,2), s_new, 'spline');
    
    outline_resampled = [x_new, y_new];

    % polyshape + centroid (in CROPPED coords)
    ann = polyshape(outline_resampled);
    [centx, centy] = centroid(ann);
    cent = [centx, centy];

    % identify commissure indicies
    % 1: AP, 2: AS, 3: SP
    ann_comm_dex = knnsearch(ann.Vertices, comm);

    % identify pin indicies
    ann_pin_dex  = knnsearch(ann.Vertices, pins);
    
    % align AP axis with global Y vector, AS at top, septal on left
    vec_AP = ann.Vertices(ann_comm_dex(2),:)-ann.Vertices(ann_comm_dex(1),:);
    theta = (3*pi/2) - atan2(vec_AP(2), vec_AP(1));

    R2 = [cos(theta) -sin(theta);
      sin(theta)  cos(theta)];

    % % Check: after transform, AS should have larger y than AP
    % AP_al = (comm(1,:) - cent) * R2.';   % comm(1,:) = AP
    % AS_al = (comm(2,:) - cent) * R2.';   % comm(2,:) = AS
    %
    % if AS_al(2) > AP_al(2)
    %     theta = theta + pi;  % flip 180 deg
    % end
    % 
    % R2 = [cos(theta) -sin(theta);
    %     sin(theta)  cos(theta)];
    
    % rotate and translate shape
    ann = translate(ann,-[centx, centy]);
    ann = rotate(ann,rad2deg(theta));

    % diameter vectors
    vec_AP = ann.Vertices(ann_comm_dex(2),:)-ann.Vertices(ann_comm_dex(1),:);
    dex_S = knnsearch(ann.Vertices, [min(ann.Vertices(:,1)), 0]);
    dex_L = knnsearch(ann.Vertices, [max(ann.Vertices(:,1)), 0]);
    vec_SL = ann.Vertices(dex_L,:)-ann.Vertices(dex_S,:);

    % clip axis
    if ~isempty(data.annotation.clipAxis_xy)
        clipAxis_xy_cropped = data.annotation.clipAxis_xy - offset;   % offset = [x-1 y-1]
        clipAxis_aligned = (clipAxis_xy_cropped - cent) * R2.';
    end
    if ~isempty(data.annotation.clipAxis_xy_2)
        clipAxis_xy_2_cropped = data.annotation.clipAxis_xy_2 - offset;   % offset = [x-1 y-1]
        clipAxis_2_aligned = (clipAxis_xy_2_cropped - cent) * R2.';
    end

    % Forward map: p_out = (p_in - cent) * R2.'
    M = R2.';                   % 2x2
    t = (-cent) * M;            % 1x2 translation in output coords
    
    A = [ M(1,1) M(1,2) 0;
          M(2,1) M(2,2) 0;
          t(1)   t(2)   1];
    
    tform = affine2d(A);
    
    [frame_aligned, ref_aligned] = imwarp(frame, imref2d([h w]), tform);

    figure('Name',baseFileName)
    imshow(frame_aligned, ref_aligned);
    hold on
    % plot(ann)
    plot(ann.Vertices(:,1),ann.Vertices(:,2),'k.')
    plot(ann.Vertices(ann_comm_dex,1),ann.Vertices(ann_comm_dex,2),'LineStyle','none','Color','g','Marker','o','MarkerFaceColor','g','MarkerSize',5)
    plot(ann.Vertices(ann_pin_dex,1),ann.Vertices(ann_pin_dex,2),'LineStyle','none','Color','c','Marker','o','MarkerFaceColor','c','MarkerSize',5)
    quiver(ann.Vertices(ann_comm_dex(1),1),ann.Vertices(ann_comm_dex(1),2),vec_AP(1),vec_AP(2),'off','LineWidth',3,'Color',[1,0,0])
    quiver(ann.Vertices(dex_S,1),ann.Vertices(dex_S,2),vec_SL(1),vec_SL(2),'off','LineWidth',3,'Color',[0,0,1])

    if ~isempty(data.annotation.clipAxis_xy)
        plot(clipAxis_aligned(:,1), clipAxis_aligned(:,2), 'm-', 'LineWidth', 2);
        plot(clipAxis_aligned(:,1), clipAxis_aligned(:,2), 'mo', 'MarkerFaceColor','m');
    end
    if ~isempty(data.annotation.clipAxis_xy_2)
        plot(clipAxis_2_aligned(:,1), clipAxis_2_aligned(:,2), 'm-', 'LineWidth', 2);
        plot(clipAxis_2_aligned(:,1), clipAxis_2_aligned(:,2), 'mo', 'MarkerFaceColor','m');
    end

    for i = 1 : length(data.annotation.commissure_labels)
        text(gca, ann.Vertices(ann_comm_dex(i),1)+5, ann.Vertices(ann_comm_dex(i),2)+5, char(data.annotation.commissure_labels(i)), ...
            'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
    end

    for i = 1 : length(data.annotation.pin_labels)
        text(gca, ann.Vertices(ann_pin_dex(i),1)+5, ann.Vertices(ann_pin_dex(i),2)+5, char(data.annotation.pin_labels(i)), ...
            'Color', 'c', 'FontSize', 12, 'FontWeight', 'bold');
    end

end

