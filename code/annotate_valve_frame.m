function annotate_valve_frame
% annotate_valve_frame
% Workflow:
% 1) Load frame from video
% 2) User zooms by drawing rectangle (crop); save cropOffset to map coords back
% 3) If filename contains "AS" or "AP" or "SP": user draws clip axis line (drawline)
% 4) Select commissures AP, AS, SP (undo/finish)
% 5) Select pins 1..8 (undo/finish)
% 6) Ask if valve closed; if not, trace gap area (drawfreehand)
% 7) Trace valve outline (final step) (drawfreehand)
% 8) Save pixel coordinates in ORIGINAL (uncropped) frame coordinates

clc;

%% ---- Labels ----
COMM_LABELS = ["AP","AS","SP"]; % click order
PIN_LABELS  = string(1:8);     % click order

%% ---- Choose video ----
[vidName, vidPath] = uigetfile({'*.mp4;*.m4v;*.mov','Video Files (*.mp4, *.m4v, *.mov)'}, ...
    'Select a video file');
if isequal(vidName,0)
    error('No video selected.');
end
videoFile = fullfile(vidPath, vidName);
v = VideoReader(videoFile);

lowerName = lower(vidName);

%% ---- Choose a frame ----
defaultFrame = 1;
prompt = sprintf('Enter frame index (# / %d). Duration = %.2f s, FPS = %.3f', v.NumFrames, v.Duration, v.FrameRate);
frameStr = inputdlg({prompt}, 'Frame selection', 1, {num2str(defaultFrame)});
if isempty(frameStr)
    error('Frame selection cancelled.');
end
frameIdx = max(1, round(str2double(frameStr{1})));
if isnan(frameIdx) || frameIdx < 1
    error('Invalid frame index.');
end

t = (frameIdx-1) / v.FrameRate;
t = min(max(t, 0), max(v.Duration - 1/v.FrameRate, 0));
v.CurrentTime = t;
frame = readFrame(v);

%% ---- Figure + UI ----
hFig = figure('Name',['Annotate ', vidName],'Color','w');
set(hFig,'KeyPressFcn',@onKeyPress);

hAx = axes('Parent', hFig);
imshow(frame, 'Parent', hAx);
axis(hAx, 'image'); hold(hAx, 'on');

btnUndo = uicontrol('Style','pushbutton','String','UNDO', ...
    'Units','normalized','Position',[0.01 0.01 0.10 0.06], ...
    'FontSize',11,'Enable','off','Callback',@(~,~)onUndo());
btnFinish = uicontrol('Style','pushbutton','String','FINISH', ...
    'Units','normalized','Position',[0.12 0.01 0.10 0.06], ...
    'FontSize',11,'Enable','off','Callback',@(~,~)onFinish());

txtStatus = uicontrol('Style','text','String','', ...
    'Units','normalized','Position',[0.24 0.01 0.75 0.06], ...
    'HorizontalAlignment','left','BackgroundColor','w','FontSize',11);

%% ---- Step 1: Zoom by rectangle (crop) ----
title(hAx, {'Drag a rectangle to zoom into the valve region.', ...
            'Double-click inside the rectangle to confirm.'});
set(txtStatus,'String','Zoom: drag rectangle (double-click to confirm).');

zoomRect = drawrectangle(hAx, 'Color','g','LineWidth',1.5);
wait(zoomRect); % double-click to confirm

rectPos = round(zoomRect.Position);  % [x y w h]

% Clamp to image bounds
x1 = max(1, rectPos(1));
y1 = max(1, rectPos(2));
x2 = min(size(frame,2), x1 + rectPos(3));
y2 = min(size(frame,1), y1 + rectPos(4));

frameCrop = frame(y1:y2, x1:x2, :);
cropOffset = [x1-1, y1-1]; % add this to map from crop->original

delete(zoomRect);

cla(hAx);
imshow(frameCrop, 'Parent', hAx);
axis(hAx,'image'); hold(hAx,'on');

%% ---- Step 2a: Select clip center ----
needsClipPoint = contains(lowerName, "as") || contains(lowerName, "ap") || contains(lowerName, "sp");
needsClipPoint2 = contains(lowerName, "asap") || contains(lowerName, "spas") || contains(lowerName, "spap");

clip_center_xy = [];
clip_center_xy_2 = [];

if needsClipPoint

title(hAx, 'Click the center of the clip(s)');
set(txtStatus,'String','Select clip center (left click).');

set(btnUndo,'Enable','off');
set(btnFinish,'Enable','off');

[clipCx, clipCy, button] = ginput(1);

if isempty(button) || button ~= 1
    error('Clip center selection cancelled.');
end

clip_center_xy = [clipCx, clipCy];

% Plot marker
plot(hAx, clipCx, clipCy, 'mo', ...
    'MarkerSize',10, ...
    'LineWidth',2, ...
    'MarkerFaceColor','m');

end

if needsClipPoint2

    title(hAx, 'Click the center of the second clip');
    set(txtStatus,'String','Select clip center (left click).');

    set(btnUndo,'Enable','off');
    set(btnFinish,'Enable','off');

    [clipCx, clipCy, button] = ginput(1);

    if isempty(button) || button ~= 1
        error('Clip center selection cancelled.');
    end

    clip_center_xy_2 = [clipCx, clipCy];

    % Plot marker
    plot(hAx, clipCx, clipCy, 'mo', ...
        'MarkerSize',10, ...
        'LineWidth',2, ...
        'MarkerFaceColor','m');

end


%% ---- Step 2b: If filename contains AS/AP/SP, ask for clip axis line ----
needsClipAxis = contains(lowerName, "as") || contains(lowerName, "ap") || contains(lowerName, "sp");
needsClipAxis2 = contains(lowerName, "asap") || contains(lowerName, "spas") || contains(lowerName, "spap");

clipAxis_xy = []; % 2x2 in ORIGINAL coords: [x y; x y]
clipAxis_xy_2 = []; % 2x2 in ORIGINAL coords: [x y; x y]

if needsClipAxis
    title(hAx, {'Draw a line along the clip axis; double-click to finish.'});
    set(txtStatus,'String','Draw clip axis line (double-click to finish).');

    hLine = drawline(hAx, 'Color','m','LineWidth',2);
    wait(hLine); % double-click to finish

    % Line endpoints are in CROPPED coords
    p = hLine.Position; % 2x2 [x y] in cropped
    clipAxis_xy = p + cropOffset; % to original
end

if needsClipAxis2
    title(hAx, {'Draw a line along the second clip axis; double-click to finish.'});
    set(txtStatus,'String','Draw clip axis line (double-click to finish).');

    hLine = drawline(hAx, 'Color','m','LineWidth',2);
    wait(hLine); % double-click to finish

    % Line endpoints are in CROPPED coords
    p = hLine.Position; % 2x2 [x y] in cropped
    clipAxis_xy_2 = p + cropOffset; % to original
end

%% ---- Step 3: Select commissures (AP, AS, SP) ----
commissures_xy_crop = selectLabeledPointsStage(COMM_LABELS, "Commissure", 'r', 'o');
commissures_xy = commissures_xy_crop + cropOffset; % to original

%% ---- Step 4a: Trace Pin Width
pin_width_xy = []; % 2x2

title(hAx, {'Draw a line along the pin width; double-click to finish.'});
set(txtStatus,'String','Draw pin width (double-click to finish).');

hLine = drawline(hAx, 'Color','c','LineWidth',2);
wait(hLine); % double-click to finish

% Line endpoints are in CROPPED coords
p = hLine.Position; % 2x2 [x y] in cropped
pin_width_xy = p + cropOffset; % to original

%% ---- Step 4: Select pins (1..8) ----
pins_xy_crop = selectLabeledPointsStage(PIN_LABELS, "Pin", 'c', 's');
pins_xy = pins_xy_crop + cropOffset; % to original

%% ---- Step 5: Ask if valve closed; if not, trace MULTIPLE gap regions (UNDO supported) ----
gap_regions_xy = {};   % cell array of polygons in ORIGINAL coords
gap_roi_handles = gobjects(0);  % handles for drawn ROIs (for undo)
gap_plot_handles = gobjects(0); % optional overlay plot handles (for undo)

choice = questdlg('Is the valve closed in this frame?', 'Valve closure', 'Yes','No','Yes');
valveClosed = strcmp(choice,'Yes');

if ~valveClosed
    title(hAx, {'Trace GAP region(s). Draw as many as needed.', ...
                'Double-click to finish each region.', ...
                'UNDO (z) removes the last region.', ...
                'Press FINISH (or Enter) when done.'});
    set(txtStatus,'String','Gap: draw region(s). z=undo last region. Press FINISH when done.');

    % Enable FINISH immediately; UNDO enabled only once at least 1 region exists
    set(btnFinish,'Enable','on');
    set(btnUndo,'Enable','off');

    % Arm finish gating + set a gate-specific undo handler
    startFinishGate();
    setGateUndo(@undoLastGapRegion);

    while ~finishGateTriggered()
        % Draw one region (user double-clicks to finish)
        hGap = drawfreehand(hAx, 'Color',[1 0.5 0], 'LineWidth',2); %#ok<NBRAK>
        wait(hGap);

        % If user pressed FINISH while interacting, allow loop to exit after processing
        poly = hGap.Position; % CROPPED coords

        if size(poly,1) >= 3
            % Close polygon
            if any(poly(1,:) ~= poly(end,:))
                poly(end+1,:) = poly(1,:);
            end

            % Optional overlay line (ROI itself is visible, but this is helpful)
            hP = plot(hAx, poly(:,1), poly(:,2), '-', 'LineWidth',2);

            % Save handles for undo
            gap_roi_handles(end+1,1) = hGap; %#ok<AGROW>
            gap_plot_handles(end+1,1) = hP;  %#ok<AGROW>

            % Store in ORIGINAL coords
            gap_regions_xy{end+1,1} = poly + cropOffset; %#ok<AGROW>

            % Enable undo now that at least one region exists
            set(btnUndo,'Enable','on');
        else
            % Too few points; discard
            if isgraphics(hGap), delete(hGap); end
        end

        drawnow;
    end

    % Disarm gate + gate undo
    clearFinishGate();
    clearGateUndo();
else
    gap_regions_xy = {};
end

function undoLastGapRegion()
    if isempty(gap_roi_handles) && isempty(gap_regions_xy)
        beep;
        set(btnUndo,'Enable','off');
        return;
    end

    % Delete most recent ROI + plot
    if ~isempty(gap_roi_handles)
        if isgraphics(gap_roi_handles(end)), delete(gap_roi_handles(end)); end
        gap_roi_handles(end) = [];
    end
    if ~isempty(gap_plot_handles)
        if isgraphics(gap_plot_handles(end)), delete(gap_plot_handles(end)); end
        gap_plot_handles(end) = [];
    end

    % Remove last saved polygon
    if ~isempty(gap_regions_xy)
        gap_regions_xy(end) = [];
    end

    if isempty(gap_regions_xy)
        set(btnUndo,'Enable','off');
        set(txtStatus,'String','Gap: draw region(s). z=undo last region. Press FINISH when done.');
    else
        set(txtStatus,'String',sprintf('Gap: %d region(s) saved. z=undo last region. Press FINISH when done.', numel(gap_regions_xy)));
    end
end


%% ---- Step 6 (final): Trace valve outline, continue only on FINISH ----
outline_xy = []; % Nx2 in ORIGINAL coords

title(hAx, {'FINAL STEP: Trace valve outline.', ...
            'Double-click to finish the outline.', ...
            'Press FINISH (or Enter) to continue.'});
set(txtStatus,'String','Outline: draw outline, double-click. Press FINISH to continue.');

set(btnUndo,'Enable','off');     % UNDO not used here
set(btnFinish,'Enable','on');    % allow immediate finishing after outline is satisfactory

startFinishGate();

% Allow re-draws until FINISH is pressed.
% We keep only the most recent outline.
hOutline = [];
while ~finishGateTriggered()
    if isgraphics(hOutline)
        delete(hOutline);
    end

    hOutline = drawfreehand(hAx, 'Color','y','LineWidth',2);
    wait(hOutline); % double-click ends outline

    poly = hOutline.Position; % CROPPED coords
    if size(poly,1) >= 3
        if any(poly(1,:) ~= poly(end,:))
            poly(end+1,:) = poly(1,:);
        end
        outline_xy = poly + cropOffset; % ORIGINAL coords
    else
        outline_xy = [];
        delete(hOutline);
    end

    set(txtStatus,'String','Outline set. Press FINISH to continue, or redraw the outline.');
    drawnow;
end

clearFinishGate();

if isempty(outline_xy)
    warning('No outline saved (outline was empty).');
end


%% ---- Package results ----
title(hAx, 'Done.');

annotation = struct();
annotation.videoFile = videoFile;
annotation.frameIdx  = frameIdx;
annotation.frameTime = t;

annotation.originalFrameSize = size(frame);
annotation.cropFrameSize = size(frameCrop);
annotation.cropOffset = cropOffset;
annotation.cropRect_xywh = [x1 y1 (x2-x1) (y2-y1)]; % original coords

annotation.clip_center_xy = clip_center_xy;
annotation.clip_center_xy_2 = clip_center_xy_2;
annotation.clipAxis_xy = clipAxis_xy; % [] or 2x2 [x y] in original coords
annotation.clipAxis_xy_2 = clipAxis_xy_2;

annotation.commissure_labels = COMM_LABELS;
annotation.commissures_xy = commissures_xy; % 3x2 original coords

annotation.pin_labels = PIN_LABELS;
annotation.pins_xy = pins_xy; % 8x2 original coords
annotation.pin_width_xy = pin_width_xy; % to get pixels to mm

annotation.valve_closed = valveClosed;
annotation.gap_regions_xy = gap_regions_xy;  % cell array of polygons in original coords

annotation.outline_xy = outline_xy;          % polygon in original coords

% Optional: store the frame (big)
annotation.frame = frame;

%% ---- Save results ----
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

%% ================== NESTED FUNCTIONS ==================
function pts_xy = selectLabeledPointsStage(labels, stageName, colorChar, markerChar)
    % Select points in order defined by labels (string array).
    N = numel(labels);

    pts_xy   = nan(N,2);
    hMarkers = gobjects(N,1);
    hLabels  = gobjects(N,1);
    k = 0;
    stageDone = false;

    % Expose stage callbacks to UI callbacks
    stage.undo   = @undoLocal;
    stage.finish = @finishLocal;
    setappdata(hFig, 'currentStage', stage);

    updateUIStage();

    while ishandle(hFig) && ~stageDone
        [x, y, button] = ginput(1);
        if ~ishandle(hFig), error('Figure closed.'); end
        if isempty(button) || button ~= 1
            continue;
        end

        if k >= N
            beep;
            continue;
        end

        k = k + 1;
        pts_xy(k,:) = [x y];

        hMarkers(k) = plot(hAx, x, y, [colorChar markerChar], ...
            'MarkerSize', 8, 'LineWidth', 2);
        hLabels(k)  = text(hAx, x+5, y, char(labels(k)), ...
            'Color', colorChar, 'FontSize', 12, 'FontWeight', 'bold');

        updateUIStage();
    end

    if isappdata(hFig,'currentStage')
        rmappdata(hFig, 'currentStage');
    end

    function updateUIStage()
        if k < N
            set(txtStatus, 'String', sprintf('Select %s %s (left click).   (z=undo, Enter=finish)', ...
                stageName, char(labels(k+1))));
        else
            set(txtStatus, 'String', sprintf('All %s points selected. Press Enter or click FINISH.', stageName));
        end
        set(btnUndo,   'Enable', ternary(k>0,'on','off'));
        set(btnFinish, 'Enable', ternary(k==N,'on','off'));

        if k < N
            title(hAx, sprintf('Select %s: %s (z=undo, Enter=finish)', stageName, char(labels(k+1))));
        else
            title(hAx, sprintf('All %s selected. Press Enter / FINISH.', stageName));
        end
    end

    function undoLocal()
        if k <= 0, return; end
        if isgraphics(hMarkers(k)), delete(hMarkers(k)); end
        if isgraphics(hLabels(k)),  delete(hLabels(k));  end
        pts_xy(k,:) = [nan nan];
        k = k - 1;
        updateUIStage();
    end

    function finishLocal()
        if k ~= N
            beep;
            set(txtStatus,'String',sprintf('Need %d %s points before finishing.', N, stageName));
            return;
        end
        stageDone = true;
    end
end

function onUndo()
    if ~ishandle(hFig), return; end

    % If a finish gate is active (gap/outline steps), UNDO calls the gate's undo handler (if any)
    if isappdata(hFig,'finish_gate_active') && getappdata(hFig,'finish_gate_active')
        if isappdata(hFig,'gate_undo_fcn')
            f = getappdata(hFig,'gate_undo_fcn');
            f(); % undo last (e.g., last gap region)
        else
            beep;
        end
        return;
    end

    % Otherwise, behave as before (point selection stage)
    if ~isappdata(hFig, 'currentStage'), return; end
    st = getappdata(hFig, 'currentStage');
    st.undo();
end

function onFinish()
    if ~ishandle(hFig), return; end

    % If a "finish gate" is active (gap/outline steps), FINISH triggers the gate
    if isappdata(hFig,'finish_gate_active') && getappdata(hFig,'finish_gate_active')
        setappdata(hFig,'finish_gate_triggered', true);
        return;
    end

    % Otherwise, behave as before (point selection stage)
    if ~isappdata(hFig, 'currentStage'), return; end
    st = getappdata(hFig, 'currentStage');
    st.finish();
end


function onKeyPress(~, event)
    switch event.Key
        case 'z'
            onUndo();
        case {'return','enter'}
            onFinish();
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function startFinishGate()
    % Enables FINISH to act as a "continue" gate (independent of point stages)
    setappdata(hFig, 'finish_gate_active', true);
    setappdata(hFig, 'finish_gate_triggered', false);
end

function tf = finishGateTriggered()
    tf = isappdata(hFig,'finish_gate_triggered') && getappdata(hFig,'finish_gate_triggered');
end

function clearFinishGate()
    if isappdata(hFig,'finish_gate_active'), rmappdata(hFig,'finish_gate_active'); end
    if isappdata(hFig,'finish_gate_triggered'), rmappdata(hFig,'finish_gate_triggered'); end
end

function setGateUndo(undoFcnHandle)
    % undoFcnHandle: function handle with signature () -> void
    setappdata(hFig, 'gate_undo_fcn', undoFcnHandle);
end

function clearGateUndo()
    if isappdata(hFig,'gate_undo_fcn')
        rmappdata(hFig,'gate_undo_fcn');
    end
end

end
