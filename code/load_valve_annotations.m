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
folderPath = uigetdir('C:\Users\ch49693\Box\Shared\In Vitro Heart Experiments\Experimental Data\Videos');

% Create pattern to identify all .mat files
filePattern = fullfile(folderPath, '*frame*.mat');
% Get information about all matching .mat files
matFiles = dir(filePattern);

% Initialize cell arrays to store data and filenames
annotationData = cell(length(matFiles), 1);
traceData = annotationData;
names = cell(length(matFiles), 1);

% Initialize data structure
AnnotatedVideoData.Valve = {};
AnnotatedVideoData.Test = {};
AnnotatedVideoData.Area = [];
AnnotatedVideoData.Perimeter = [];
AnnotatedVideoData.SL_Diameter = [];
AnnotatedVideoData.AP_Diameter = [];

AnnotatedVideoData.AnnulusOutline = [];
AnnotatedVideoData.Comm_Idx = [];
AnnotatedVideoData.SL_Idx = [];
AnnotatedVideoData.Pin_Idx = [];
AnnotatedVideoData.ClipCenter = [];
AnnotatedVideoData.ClipAxis_Idx = [];
AnnotatedVideoData.ClipCenter2 = [];
AnnotatedVideoData.ClipAxis2_Idx = [];

AnnotatedVideoData.Gap_Area = [];
AnnotatedVideoData.Frame = {};
AnnotatedVideoData.Ref = {};
AnnotatedVideoData.mmPerPixel = [];

load_clip_centers = false;

if load_clip_centers
    load(fullfile(folderPath,'ClipCenters.mat'))
end

%% Loop through each .mat file
for k = 1:length(annotationData)
    % Grab Valve Name
    folderFullPath = strsplit(matFiles(k).folder,'\');
    testDate = strsplit(folderFullPath{end},'_');
    AnnotatedVideoData(k).Valve = [testDate{1},'_',testDate{2},'_',testDate{3}];

    baseFileName = matFiles(k).name;
    fullFileName = fullfile(folderPath, baseFileName);
    fprintf('Now reading %s\n', fullFileName);

    testParts = strsplit(baseFileName,'_');
    AnnotatedVideoData(k).Test = testParts{1};
    
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
    x_new = interp1(s, outline(:,1), s_new, 'pchip');
    y_new = interp1(s, outline(:,2), s_new, 'pchip');
    
    outline_resampled = [x_new, y_new];

    % DO NOT RESAMPLE
    % outline_resampled = outline;

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
        [clipAxis_IntersectionPoints, segIdx, bxy] = clipAxisAnnulusIntersection(ann, clipAxis_aligned);
        clipAxis_Idx = knnsearch(ann.Vertices, clipAxis_IntersectionPoints);
    else
        clip_center_xy = [];
        clipAxis_Idx = [];
    end
    if ~isempty(data.annotation.clipAxis_xy_2)
        clipAxis_xy_2_cropped = data.annotation.clipAxis_xy_2 - offset;   % offset = [x-1 y-1]
        clipAxis_2_aligned = (clipAxis_xy_2_cropped - cent) * R2.';
        [clipAxis2_IntersectionPoints, segIdx2, bxy2] = clipAxisAnnulusIntersection(ann, clipAxis_2_aligned);
        clipAxis2_Idx = knnsearch(ann.Vertices, clipAxis2_IntersectionPoints);
    else
        clip_center_xy_2 = [];
        clipAxis2_Idx = [];
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
        plot(ann.Vertices(clipAxis_Idx,1),ann.Vertices(clipAxis_Idx,2),'m--')
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

    % Select clip center if not saved

    if ~isempty(data.annotation.clipAxis_xy)
        if ~load_clip_centers
            [clipCx, clipCy, button] = ginput(1);
            clip_center_xy = [clipCx, clipCy];
            ClipCenter(k,:) = clip_center_xy;
        else
            clip_center_xy = ClipCenter(k,:);
        end
        % Plot marker
        plot(clip_center_xy(1), clip_center_xy(2), 'mo', ...
            'MarkerSize',10, ...
            'LineWidth',2, ...
            'MarkerFaceColor','m');
    end

    if ~isempty(data.annotation.clipAxis_xy_2)
        if ~load_clip_centers
            [clipCx2, clipCy2, button] = ginput(1);
            clip_center_xy_2 = [clipCx2, clipCy2];
            ClipCenter2(k,:) = clip_center_xy_2;
        else
            clip_center_xy_2 = ClipCenter2(k,:);
        end
        % Plot marker
        plot(clip_center_xy_2(1), clip_center_xy_2(2), 'mo', ...
            'MarkerSize',10, ...
            'LineWidth',2, ...
            'MarkerFaceColor','m');
    end

    pixdist = pdist(data.annotation.pin_width_xy);
    mmPERpix = 7/pixdist; % pins are 7 mm in width

    AnnotatedVideoData(k).Area = (mmPERpix^2).*area(ann);
    AnnotatedVideoData(k).Perimeter = mmPERpix.*perimeter(ann);
    AnnotatedVideoData(k).SL_Diameter = mmPERpix.*norm(vec_SL);
    AnnotatedVideoData(k).AP_Diameter = mmPERpix.*norm(vec_AP);

    AnnotatedVideoData(k).AnnulusOutline = ann;
    AnnotatedVideoData(k).Comm_Idx= ann_comm_dex;
    AnnotatedVideoData(k).SL_Idx = [dex_S,dex_L];  
    AnnotatedVideoData(k).Pin_Idx = ann_pin_dex;
    AnnotatedVideoData(k).ClipCenter = clip_center_xy;
    AnnotatedVideoData(k).ClipAxis_Idx = clipAxis_Idx;
    AnnotatedVideoData(k).ClipCenter2 = clip_center_xy_2;
    AnnotatedVideoData(k).ClipAxis2_Idx = clipAxis2_Idx;

    % compute coaptation gap area
    gap_area = 0;
    if ~isempty(data.annotation.gap_regions_xy)
        for i = 1 : length(data.annotation.gap_regions_xy)
            try
                temp_gap_area = area(polyshape(data.annotation.gap_regions_xy{i}));
            catch
                temp_gap_area = 0;
            end
            gap_area = gap_area + (mmPERpix^2).*temp_gap_area;
        end
    end
        
    AnnotatedVideoData(k).Gap_Area = gap_area;

    AnnotatedVideoData(k).Frame = frame_aligned;
    AnnotatedVideoData(k).Ref = ref_aligned;
    AnnotatedVideoData(k).mmPerPixel = mmPERpix;

end

AnnotatedVideoDataTable = struct2table(AnnotatedVideoData);
% 
% save(fullfile(folderPath,'AnnotatedVideoData.mat'),'AnnotatedVideoData')
% save(fullfile(folderPath,'ClipCenters.mat'),'ClipCenter','ClipCenter2')

%% HELPER FUNCTIONS
function [P, segIdx, bxy] = clipAxisAnnulusIntersection(ann, clipAxis_xy, varargin)
% clipAxisAnnulusIntersection
% Intersect a clip axis (line segment) with an annulus polyshape boundary.
%
% Inputs
%   ann           : polyshape (annulus)
%   clipAxis_xy   : 2x2 [x y; x y] endpoints of the clip axis segment
%
% Name-value options
%   'Tol'         : numeric tolerance for dedup (default 1e-6)
%   'KeepTwoExtremes' : true/false (default true) keep the two points farthest
%                      apart along the clip direction if >2 intersections occur
%
% Outputs
%   P      : Kx2 intersection points (in same coords as ann/clipAxis_xy)
%   segIdx : Kx1 boundary segment indices intersected (segment i is bxy(i)->bxy(i+1))
%   bxy    : Mx2 boundary vertices used (closed: bxy(end,:)=bxy(1,:))

p = inputParser;
p.addParameter('Tol', 1e-6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('KeepTwoExtremes', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
tol = p.Results.Tol;
keepTwo = p.Results.KeepTwoExtremes;

% Boundary of polyshape (outer boundary). For annulus you likely have only one boundary.
bxy = ann.Vertices; % polyshape stores boundary vertices here for simple polygons
if isempty(bxy)
    P = zeros(0,2); segIdx = zeros(0,1); return;
end

% Ensure closed boundary
if any(bxy(1,:) ~= bxy(end,:))
    bxy(end+1,:) = bxy(1,:);
end

clipAxis_xy = extendClipAxisToAnnulus(ann, clipAxis_xy, 2); % extend to cross annulus
A = clipAxis_xy(1,:);
B = clipAxis_xy(2,:);

P = zeros(0,2);
segIdx = zeros(0,1);

for i = 1:size(bxy,1)-1
    C = bxy(i,:);
    D = bxy(i+1,:);
    [hit, pint] = segmentIntersectPoint2D(A,B,C,D);
    if hit
        P(end+1,:) = pint; %#ok<AGROW>
        segIdx(end+1,1) = i; %#ok<AGROW>
    end
end

% Deduplicate points (common when line hits a vertex)
if ~isempty(P)
    keep = true(size(P,1),1);
    for i = 1:size(P,1)
        if ~keep(i), continue; end
        d = sqrt(sum((P - P(i,:)).^2,2));
        same = d < tol;
        same(i) = false;
        keep(same) = false;
    end
    P = P(keep,:);
    segIdx = segIdx(keep,:);
end

% If more than 2 intersections (nonconvex boundary), keep the extremes along the clip axis
if keepTwo && size(P,1) > 2
    dir = (B - A);
    dir = dir / norm(dir);
    s = (P - A) * dir(:);        % scalar projection along clip direction
    [~, order] = sort(s);
    use = [order(1); order(end)];
    P = P(use,:);
    segIdx = segIdx(use,:);
end

end

function [hit, P] = segmentIntersectPoint2D(A,B,C,D)
% Segment-segment intersection in 2D returning the point when they intersect.
hit = false;
P = [NaN NaN];

r = B - A;
s = D - C;

den = cross2(r, s);
qmp = (C - A);

if abs(den) < 1e-12
    % Parallel/colinear: ignore for typical clip-axis use
    return;
end

t = cross2(qmp, s) / den;
u = cross2(qmp, r) / den;

if t >= 0 && t <= 1 && u >= 0 && u <= 1
    hit = true;
    P = A + t*r;
end
end

function z = cross2(a,b)
z = a(1)*b(2) - a(2)*b(1);
end

function clipExt = extendClipAxisToAnnulus(ann, clipAxis_xy, marginFactor)
% ann: polyshape
% clipAxis_xy: 2x2 endpoints (often inside annulus)
% marginFactor: e.g. 2 (default) enlarges bbox diagonal

if nargin < 3, marginFactor = 2; end

A = clipAxis_xy(1,:);
B = clipAxis_xy(2,:);

v = B - A;
nv = norm(v);
if nv < 1e-9
    error('Clip axis points are identical or too close.');
end
u = v / nv;

% Bounding box of annulus boundary
vx = ann.Vertices(:,1);
vy = ann.Vertices(:,2);
xmin = min(vx); xmax = max(vx);
ymin = min(vy); ymax = max(vy);

diagLen = hypot(xmax - xmin, ymax - ymin);
L = marginFactor * diagLen;   % extension length each direction

mid = (A + B)/2;
P1 = mid - L*u;
P2 = mid + L*u;

clipExt = [P1; P2];
end

