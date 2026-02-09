%% annotate_valve_frame.m
% Loads a single frame from an .mp4, lets user trace valve outline,
% select 3 commissures + 8 pins, then saves pixel coordinates.

clear; clc;

%% ---- User settings ----
NUM_COMMISSURES = 3;
NUM_PINS = 8;

%% ---- Choose video ----
[vidName, vidPath] = uigetfile({'*.mp4;*.m4v;*.mov','Video Files (*.mp4, *.m4v, *.mov)'}, ...
    'Select a video file');
if isequal(vidName,0)
    error('No video selected.');
end
videoFile = fullfile(vidPath, vidName);

v = VideoReader(videoFile);

%% ---- Choose a frame ----
% Prompt for frame number (1..NumFrames when available). If NumFrames is not
% reliable in your MATLAB version/codecs, we also allow time-based reading.
defaultFrame = 1;
frameStr = inputdlg({sprintf('Enter frame index (>=1). Video duration = %.2f s, FPS = %.3f', ...
    v.Duration, v.FrameRate)}, ...
    'Frame selection', 1, {num2str(defaultFrame)});
if isempty(frameStr)
    error('Frame selection cancelled.');
end
frameIdx = max(1, round(str2double(frameStr{1})));
if isnan(frameIdx) || frameIdx < 1
    error('Invalid frame index.');
end

% Try frame-index based access via CurrentTime
t = (frameIdx-1) / v.FrameRate;
t = min(max(t, 0), max(v.Duration - 1/v.FrameRate, 0));
v.CurrentTime = t;

frame = readFrame(v);

%% ---- Display frame ----
hFig = figure('Name','Valve annotation','Color','w');
imshow(frame);
axis image; hold on;
title({'Trace valve outline, double-click to finish.', ...
       'After outline: select 3 commissures, then 8 pins.'});

%% ---- Trace valve outline ----
% Freehand ROI; user double-clicks to finish
roi = drawfreehand('Color','y','LineWidth',2);  % requires R2018b+
outline_xy = roi.Position;  % Nx2 [x y] in pixels

% Close the outline explicitly (optional)
if ~isempty(outline_xy) && any(outline_xy(1,:) ~= outline_xy(end,:))
    outline_xy(end+1,:) = outline_xy(1,:);
end
plot(outline_xy(:,1), outline_xy(:,2), 'y-', 'LineWidth', 2);

%% ---- Select commissures ----
title(sprintf('Click %d commissure points (in order).', NUM_COMMISSURES));
[comm_x, comm_y] = ginput(NUM_COMMISSURES);
commissures_xy = [comm_x(:), comm_y(:)];
plot(commissures_xy(:,1), commissures_xy(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Label commissures
for i = 1:NUM_COMMISSURES
    text(commissures_xy(i,1)+5, commissures_xy(i,2), sprintf('C%d', i), ...
        'Color','r','FontSize',12,'FontWeight','bold');
end

%% ---- Select pin locations ----
title(sprintf('Click %d pin locations (in order).', NUM_PINS));
[pin_x, pin_y] = ginput(NUM_PINS);
pins_xy = [pin_x(:), pin_y(:)];
plot(pins_xy(:,1), pins_xy(:,2), 'cs', 'MarkerSize', 8, 'LineWidth', 2);

% Label pins
for i = 1:NUM_PINS
    text(pins_xy(i,1)+5, pins_xy(i,2), sprintf('P%d', i), ...
        'Color','c','FontSize',12,'FontWeight','bold');
end

title('Done. Close figure or keep open.');

%% ---- Package results ----
annotation = struct();
annotation.videoFile = videoFile;
annotation.frameIdx  = frameIdx;
annotation.frameTime = t;
annotation.frameSize = size(frame); % [H W C]
annotation.outline_xy = outline_xy;         % Nx2 [x y]
annotation.commissures_xy = commissures_xy; % 3x2 [x y]
annotation.pins_xy = pins_xy;               % 8x2 [x y]

% Optional: store the frame itself (can be big; comment out if undesired)
annotation.frame = frame;

%% ---- Save results ----
% Default save name
[~, base, ~] = fileparts(vidName);
defaultMat = fullfile(vidPath, sprintf('%s_frame%06d_annotation.mat', base, frameIdx));

[saveName, savePath] = uiputfile('*.mat', 'Save annotation as', defaultMat);
if isequal(saveName,0)
    warning('Save cancelled. Results are in workspace variable "annotation".');
else
    saveFile = fullfile(savePath, saveName);
    save(saveFile, 'annotation');
    fprintf('Saved annotation to:\n  %s\n', saveFile);
end

%% ---- (Optional) also save JSON (outline + points only) ----
choice = questdlg('Also save a JSON (outline + points)?', 'Export JSON', 'Yes','No','No');
if strcmp(choice,'Yes')
    jsonStruct = annotation;
    jsonStruct = rmfield(jsonStruct, 'frame'); % don't embed image in JSON

    jsonText = jsonencode(jsonStruct);
    jsonText = prettyjson(jsonText); % helper below

    jsonFile = fullfile(savePath, sprintf('%s_frame%06d_annotation.json', base, frameIdx));
    fid = fopen(jsonFile, 'w');
    fwrite(fid, jsonText, 'char');
    fclose(fid);
    fprintf('Saved JSON to:\n  %s\n', jsonFile);
end

%% ---- Helper: pretty print JSON (simple) ----
function out = prettyjson(in)
% Minimal JSON pretty-printer (good enough for small structs)
indent = 0;
out = "";
i = 1;
while i <= strlength(in)
    ch = extractBetween(in, i, i);
    if ch == "{"
        indent = indent + 1;
        out = out + "{\n" + repmat("  ",1,indent);
    elseif ch == "}"
        indent = max(indent - 1, 0);
        out = out + "\n" + repmat("  ",1,indent) + "}";
    elseif ch == "["
        indent = indent + 1;
        out = out + "[\n" + repmat("  ",1,indent);
    elseif ch == "]"
        indent = max(indent - 1, 0);
        out = out + "\n" + repmat("  ",1,indent) + "]";
    elseif ch == ","
        out = out + ",\n" + repmat("  ",1,indent);
    elseif ch == ":"
        out = out + ": ";
    else
        out = out + ch;
    end
    i = i + 1;
end
out = char(out);
end
