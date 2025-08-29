function [fig, varargout] = pubPlot(NVArgs)
%PUBPLOT  Format the current figure for publication-quality output.
%
%   fig = PUBPLOT(opts) rescales and trims whitespace around axes; aligns
%   tick labels, axis labels, and titles; and (optionally) exports an EPS
%   with clean fonts and an exact bounding box. Subplots are laid out by
%   inferring a grid from current axes positions and computing per-band
%   margins that account for tick/label sizes.
%
%   INPUTS (Name-Value Pairs passed and stored in NVArgs structure)
%   ----------------------------------------------------------------------
%   Figure              : (figure handle) figure to format. Default: gcf
%   Filename            : (char) "Filename" to export. Default: ''
%   FileExtension   : {char} .FilenameExtension to export. Default: '.eps'
%   Width               : {'single','1.5','double','custom'} column width. Default: 'single'
%   Journal             : {'Elsevier','Springer','Nature','Wiley'} page model. Default: 'Elsevier'
%   Height              : (double) overall figure height [pt]. Default: 235
%   SpacingOffset       : (double) base padding [pt] used around ticks/labels. Default: 1
%   AxisOffset          : (double) [numPlots x 4] of [left, right, bottom,
%                                   top] offsets for each plot
%
%   FontSize            : (double) legend/colorbar/title font size. Default: 12
%   FontName            : (char) text font (e.g., 'Arial'). Default: 'Arial'
%   FontColor           : (3x1 double) RGB for titles/legends/cbars. Default: [0 0 0]
%   AxisFontSize        : (double) axis label and tick font size. Default: 12
%   AxisColor           : (3x1 double) RGB for axes/ticks. Default: [0 0 0]
%   Interpreter         : {'none','tex','latex'} text interpreter. Default: 'tex'
%
%   LineWidth           : (double) plotted line width [pt]. Default: 2
%   MarkerSize          : (double) marker size [pt]. Default: 8
%   AxisLineWidth       : (double) axes line/tick width [pt]. Default: 1.8
%   TickLength          : (2x1 double) desired tick length [pt, pt] (used isotropically). Default: [3, 1.5]
%   BoxColor            : (color) axes background. Default: 'none'
%   MoveTitleAboveAxis  : (logical) if true, title is always above axis line (bottom-aligned). Default: false
%
%   Grid                : ('on'|'off') grid visibility. Default: 'off'
%   GridLineWidth       : (double) grid line width [pt]. Default: 1
%   GridColor           : (3x1 double) grid color. Default: [0.15 0.15 0.15]
%
%   ToolBar             : {'none','auto','figure'} figure toolbar mode. Default: 'none'
%   ExportMessage       : {'on','off'} message after export. Default: 'on'
%   verbosity           : {0..4} log level (≥2 keeps ExportMessage='on'). Default: 2
%
%   figOptions          : Options to control this figures appearance. Will
%                         override everything but Figure, Filename, and Extension
%   CustomWidth         : If Width is set to 'custom', the numeric value of
%                         the width desired
%
%   RETURNS
%   ----------------------------------------------------------------------
%   fig                 : (figure handle) formatted figure
%   
%   OPTIONAL RETURNS
%   ----------------------------------------------------------------------
%   figOptions          : Settings used to generate this figure. Useful
%                         when passing to other figures to keep consistency
%
%   NOTES
%   ----------------------------------------------------------------------
%   * Legends are not repositioned and not recommended.
%   * Only 2-D Cartesian axes are supported.
%   * Axis limits/ticks are preserved; text objects are re-positioned.
%   * Tick labels are redrawn as text in points to control spacing tightly.
%   * Subplots: margins are computed per row/column band to make labels/ticks
%     fit without collisions, then each axes is placed inside its spanned
%     band using those margins.
%
%   Author : Collin Haese
%   Date   : 16-August-2025
%   Ver    : 1.1.3 Alpha
%
% 
% TODO:
% Tex interpreter breaks eps export?
% Only supports RGB triplet colors
% Validate font names
% Subplotting
%   Check support for tiledLayout
% Manually set text sizes to account for extra whitespace
% Force board size on eps
% Exponential Ticks - currently lose the order of magnitude
% Colorbars not supported
% Rotated ticks still seems a bit off (recommended to use 45 deg)


%----------------------------- Arguments ---------------------------------
arguments
    % Figure properties
    NVArgs.Figure (1,1) matlab.ui.Figure = gcf
    NVArgs.Filename (1,:) char = ''
    NVArgs.FileExtension {mustBeMember(NVArgs.FileExtension, {'.jpg','.jpeg','.png','.tif','.tiff','.gif','.svg','.pdf','.emf','.eps'})} = '.eps'
    NVArgs.Width (1,:) char {mustBeMember(NVArgs.Width, {'single','1.5','double','custom'})} = 'single'
    NVArgs.CustomWidth (1,1) double {numericPositiveScalar} = 235
    NVArgs.Journal (1,:) char {journalList} = 'Elsevier'
    NVArgs.Height (1,1) double {numericPositiveScalar} = 235
    NVArgs.SpacingOffset (1,1) double {numericPositiveScalar} = 1
    NVArgs.AxisOffset double = []

    % Font properties
    NVArgs.FontSize (1,1) double {numericPositiveScalar} = 12
    NVArgs.FontName (1,:) char = 'Arial'
    NVArgs.FontColor (3,1) double = [0 0 0]
    NVArgs.AxisFontSize (1,1) double {numericPositiveScalar} = 12
    NVArgs.AxisColor (3,1) double = [0 0 0]
    NVArgs.Interpreter (1,:) char {mustBeMember(NVArgs.Interpreter, {'none','tex','latex'})} = 'tex'

    % Plot properties
    NVArgs.LineWidth (1,1) double {numericPositiveScalar} = 2
    NVArgs.MarkerSize (1,1) double {numericPositiveScalar} = 8
    NVArgs.AxisLineWidth (1,1) double {numericPositiveScalar} = 1.8
    NVArgs.TickLength (2,1) double {mustBeNonnegative} = [3, 1.5]
    NVArgs.BoxColor = 'none'
    NVArgs.MoveTitleAboveAxis (1,1) logical = false

    % Grid properties
    NVArgs.Grid matlab.lang.OnOffSwitchState = 'off'
    NVArgs.GridLineWidth (1,1) double {numericPositiveScalar} = 1
    NVArgs.GridColor (3,1) double = [0.15 0.15 0.15]

    % Optional
    NVArgs.ToolBar (1,:) char {mustBeMember(NVArgs.ToolBar, {'none','auto','figure'})} = 'none'
    NVArgs.ExportMessage (1,:) char {mustBeMember(NVArgs.ExportMessage, {'on','off'})} = 'on'
    NVArgs.verbosity (1,1) {ismember(NVArgs.verbosity, [0,1,2,3,4])} = 2

    % Override defaults
    NVArgs.figOptions struct = struct()
end

if ~isempty(NVArgs.figOptions)
    % Overwrite NVArgs fields with values from figOptions (case-insensitive),
    % except for Figure/Filename/FileExtension. Empty values in figOptions are ignored.
    figOptions = NVArgs.figOptions;

    exclude = lower({'Figure','Filename','FileExtension'});

    fnA = fieldnames(NVArgs);
    fnB = fieldnames(figOptions);

    % Case-insensitive intersection, but keep original names for assignment
    [commonLower, ia, ib] = intersect(lower(fnA), lower(fnB));

    for k = 1 : numel(ia)
        nameA = fnA{ia(k)};    % NVArgs field (original case)
        nameB = fnB{ib(k)};    % figOptions field (original case)

        if ismember(lower(nameA), exclude)
            continue
        end

        val = figOptions.(nameB);
        if ~isempty(val)  % only override when figOptions supplies a value
            try
                NVArgs.(nameA) = val;
            catch ME
                warning('Skipping field "%s": %s', nameA, ME.message);
            end
        end
    end
end

if NVArgs.verbosity < 2
    NVArgs.ExportMessage = 'off';
end

%% ============================ SETUP ====================================
% get figure handle
fig = NVArgs.Figure;

% Collect axes. If multiple, treat as subplots (in any layout).
axList = findall(fig, 'Type', 'axes');
if isempty(axList)
    error('pubPlot:NoAxes','No axes detected in the figure.');
end
isSubplotFigure = numel(axList) > 1; % if a subplot, will have multiple axes

% ensure all axes are cartesian
axesListCart = axList(arrayfun(@(a) isprop(a,'XLim') && isprop(a,'YLim'), axList));
if any(size(axesListCart) ~= size(axList))
    error('Error: Non-Cartesian axes detected.')
end

% Only allow 2-D Cartesian axes
axList = flip(axList); % so first subplot is first element
isCartesian = arrayfun(@(a) isprop(a,'XLim') && isprop(a,'YLim'), axList);
if ~all(isCartesian)
    error('pubPlot:NonCartesian','Only 2-D Cartesian axes are supported.');
end

% Units in points for consistent measurements
for ax = axList'
    if ~isempty(get(ax,'ZLim')) && ~isequal(get(ax,'ZLim'),[-1 1])
        warning('pubPlot:3DDetected');
        disp('3D axes detected—only 2D is supported.')
    end
    ax.Units = 'points';
end

% if custom offset is supplied
if ~isempty(NVArgs.AxisOffset)
    if any(size(NVArgs.AxisOffset) ~= numel(axList),4)
        error('AxisOffsets must by numAxis x 4.')
    end
else
    AxisOffset = [];
end

%---------------------- FIGURE-LEVEL FORMATTING --------------------------

[fig, figW, figH] = formatFigure(fig,NVArgs);

%------------------------- AXES PROCESSING -------------------------------
if isSubplotFigure

    % Tolerance (points) for grouping axes in same row/column band
    tol = 1;

    % Positions (in points)
    lefts   = arrayfun(@(a) a.Position(1), axList);
    bottoms = arrayfun(@(a) a.Position(2), axList);
    widths  = arrayfun(@(a) a.Position(3), axList);
    heights = arrayfun(@(a) a.Position(4), axList);

    % Round to decimal points to keep band edges tidy
    lefts = round(lefts,1);
    bottoms = round(bottoms,1);
    widths = round(widths,1);
    heights = round(heights,1);

    rights = lefts + widths;
    tops   = bottoms + heights;

    % Unique band edges
    xStarts = uniquetol(sort(lefts(:)),  tol/max(lefts(:)));
    xEnds   = uniquetol(sort(rights(:)), tol/max(rights(:)));
    yStarts = uniquetol(sort(bottoms(:)),tol/max(bottoms(:)));
    yEnds   = uniquetol(sort(tops(:)),   tol/max(tops(:)));

    numCols = numel(xStarts); % cell columns
    numRows = numel(yStarts); % cell rows
    
    % Map each axes to start/end band indices
    colStart = zeros(numel(axList),1);
    colEnd   = zeros(numel(axList),1);
    rowStart = zeros(numel(axList),1);
    rowEnd   = zeros(numel(axList),1);
    for k = 1:numel(axList)
        colStart(k) = find(abs(xStarts - lefts(k))  <= tol, 1,'first');
        colEnd(k)   = find(abs(xEnds   - rights(k)) <= tol, 1,'first');
        rowStart(k) = find(abs(yStarts - bottoms(k))<= tol, 1,'first');
        rowEnd(k)   = find(abs(yEnds   - tops(k))   <= tol, 1,'first');
    end
    colSpan = colEnd - colStart + 1; % number of cells spanned horizontally
    rowSpan = rowEnd - rowStart + 1; % number of cells spanned vertically

    % Per-band margins needed to clear ticks/labels/titles
    colMarginL = zeros(1,numCols);  % left margin of each column band
    colMarginR = zeros(1,numCols);  % right margin of each column band
    rowMarginB = zeros(1,numRows);  % bottom margin of each row band
    rowMarginT = zeros(1,numRows);  % top margin of each row band
    moveTitleRequired = false(1,numel(axList)); % advise lifting title above axis

    % measure text extents for each axes at its *current* size as a proxy
    for k = 1 : numel(axList)
        [L,R,B,T,moveTitleRequired(k)] = getTextExtents(axList(k), NVArgs, widths(k), heights(k));
        % apply left margin to the starting column only
        colMarginL(colStart(k)) = max(colMarginL(colStart(k)), L);
        % apply right margin to the ending column only
        colMarginR(colEnd(k))   = max(colMarginR(colEnd(k)),   R);
        % apply bottom margin to the starting row only
        rowMarginB(rowStart(k)) = max(rowMarginB(rowStart(k)), B);
        % apply top margin to the ending row only
        rowMarginT(rowEnd(k))   = max(rowMarginT(rowEnd(k)),   T);
    end

    % If any axes in a row needs its title moved, move all ending in that row
    if any(moveTitleRequired)
        for k = 1 : numel(axList)
            if moveTitleRequired(k)
                moveTitleRequired(rowEnd == rowEnd(k)) = true;
            end
        end
    end

    % Inner budgets and band sizes
    innerBudgetW = figW - sum(colMarginL + colMarginR);
    innerBudgetH = figH - sum(rowMarginB + rowMarginT);
    if innerBudgetW <= 0 || innerBudgetH <= 0
        error('pubPlot:TooSmall','Figure is too small for labels/ticks.');
    end

    unitInnerW = innerBudgetW / numCols;
    unitInnerH = innerBudgetH / numRows;

    bandW = colMarginL + unitInnerW + colMarginR; % width of each column band
    bandH = rowMarginB + unitInnerH + rowMarginT; % height of each row band

    % Absorb rounding so the last band ends exactly at the figure edge
    bandW(end) = figW  - sum(bandW(1:end-1));
    bandH(end) = figH - sum(bandH(1:end-1));
    bandW = max(bandW,0); bandH = max(bandH,0);

    % Cumulative band origins
    x0 = [0, cumsum(bandW(1:end-1))];
    y0 = [0, cumsum(bandH(1:end-1))];

    % Place and trim each axes
    for k = 1 : numel(axList)
        ax = axList(k);

        % full band span for this axes
        c0 = colStart(k); c1 = colEnd(k);
        r0 = rowStart(k); r1 = rowEnd(k);

        % total outer size of this spanned block
        subW = sum(colMarginL(c0:c1)) + colSpan(k)*unitInnerW + sum(colMarginR(c0:c1));
        subH = sum(rowMarginB(r0:r1)) + rowSpan(k)*unitInnerH + sum(rowMarginT(r0:r1));

        % margins for the spanned block: take edge-most bands
        margins = [colMarginL(c0), colMarginR(c1), rowMarginB(r0), rowMarginT(r1)];
        % check data
        if ~isnumeric(margins) || ~any(margins >= 0) || ~any(margins <= 1e5)
            error('The value of margins must be a nonnegative scalar.');
        end

        % absolute origin of the spanned block
        xOff = x0(c0);
        yOff = y0(r0);

        % format and trim
        if ~isempty(NVArgs.AxisOffset)
            ax = formatAxis(ax, NVArgs);
            [ax, NVArgs.AxisOffset(k,:)] = reduceWhitespace(ax, NVArgs, subW, subH, ...
                xOffset=xOff, yOffset=yOff, ...
                customOffset=NVArgs.AxisOffset(k,:), moveTitle=moveTitleRequired(k));
        else
            ax = formatAxis(ax, NVArgs);
            [ax, AxisOffset] = reduceWhitespace(ax, NVArgs, subW, subH, ...
                xOffset=xOff, yOffset=yOff, ...
                customOffset=margins, moveTitle=moveTitleRequired(k));
        end

    end

else
        % Single axes, only one plot for this Figure
        ax = formatAxis(axList(1), NVArgs);
        if ~isempty(NVArgs.AxisOffset)
            [ax, NVArgs.AxisOffset] = reduceWhitespace(ax, NVArgs, figW, figH, customOffset=NVArgs.AxisOffset);
        else
            [ax, AxisOffset] = reduceWhitespace(ax, NVArgs, figW, figH);
        end
end

%-------------------------- EXPORT FIGURE ---------------------------------
if ~isempty(NVArgs.Filename)
    Filename = NVArgs.Filename;
    FileExtensions = NVArgs.FileExtension;

    if isa(FileExtensions,'char')
        % Vector export with exact bounding box
        exportgraphics(fig, [Filename,FileExtensions], 'ContentType','vector');
    elseif isa(FileExtensions,'cell')
        for i = 1 : numel(FileExtensions)
            % Vector export with exact bounding box
            exportgraphics(fig, [Filename,FileExtensions{i}], 'ContentType','vector');
        end
    else
        warning('Issue exporting graphic.')
    end
    
    if contains(FileExtensions,'.eps')
        % Modify postscript file to enable easy opening in Illustrator
        % Lightly post-process for font embedding/renaming and artboard size

        % read all lines from original postscript file
        lines = readlines([Filename,'.eps']);

        % store updated lines in new array
        out = string();
        linesNewCounter = 1;

        % list of fonts to export
        fontList = string();
        if strcmp(NVArgs.FontName,'Arial')
            fontList = {'Arial';'Arial-Bold';'Arial-Italic';'Arial-BoldItalic'};
        end

        inFontDict = 0; % font dictionary flag
        inReencode = 0; % font reencoding flag
        replaceFontCmd = 0; % place font command

        % loop through all lines and modify
        for i = 1 : numel(lines)
            ln = lines(i);

            % Parse markers
            if strcmp(ln,'%FOPBeginFontDict')
                inFontDict = true;
                % copy BeginFontDict line
                out(end+1) = ln;
                % manually update font reencoding commands
                for j = 1 : numel(fontList)
                    out(end+1) = "%%IncludeResource: font " + fontList(j);
                end
            elseif strcmp(ln,'%FOPEndFontDict')
                inFontDict = false;
            elseif strcmp(ln,'%FOPBeginFontReencode')
                inReencode = true;
                % copy BeginFontReencode line
                out(end+1) = ln;
                % manually update font reencoding commands
                for j = 1 : size(fontList,1)
                    out(end+1) = "/"+fontList(j)+" findfont";
                    out(end+1) = "dup length dict begin";
                    out(end+1) = "  {1 index /FID ne {def} {pop pop} ifelse} forall";
                    out(end+1) = "  /Encoding WinAnsiEncoding def";
                    out(end+1) = "  currentdict";
                    out(end+1) = "end";
                    out(end+1) = "/"+fontList(j)+" exch definefont pop";
                end
                continue
            elseif strcmp(ln,'%FOPEndFontReencode')
                inReencode = false;
                continue
            end

            if contains(ln,'/Helvetica') && ~inFontDict && ~inReencode
                if ~inReencode
                    % skip command during reencode
                    replaceFontCmd = true;
                end
            end

            % Copy lines + tweak bounding box
            if ~inFontDict && ~inReencode && ~replaceFontCmd
                % set artboard size
                if contains(ln,'%%BoundingBox')
                    ln = sprintf('%%BoundingBox:     0     0   %d   %d', round(figW), round(figH));
                end
                out(end+1) = ln; %#ok<AGROW>
                continue
            end

            if replaceFontCmd
                % change depending on font type
                if contains(lines{i,1},'Bold') && (contains(lines{i,1},'Italic') || contains(lines{i,1},'Oblique'))
                    out(end+1) = replace(ln,"Helvetica-BoldOblique","Arial-BoldItalic");
                elseif contains(lines{i,1},'Italic') || contains(lines{i,1},'Oblique')
                    out(end+1) = replace(ln,"Helvetica-Oblique","Arial-Italic");
                elseif contains(lines{i,1},'Bold')
                    out(end+1) = replace(ln,"Helvetica-Bold","Arial-Bold");
                else
                    out(end+1) = replace(ln,"Helvetica","Arial");
                end
                replaceFontCmd = false;
            end

        end

        %--------------------------- REWRITE FILE -----------------------------
        fid = fopen([Filename,'.eps'],'w');
        for i = 1 : numel(out)
            fprintf(fid,'%s\n',out(i));
        end
        fclose(fid);

        if strcmpi(NVArgs.ExportMessage,'on')
            fprintf('Successfully updated and exported %s.\n',Filename)
        end

    end

end

if nargout > 1
    NVArgs.AxisOffset = AxisOffset;
    varargout{1} = NVArgs;
end

end

% ========================= END MAIN FUNCTION ============================

%% ========================= VALIDATION HELPERS ===========================

function numericPositiveScalar(x, argName)
%NUMERICPOSITIVESCALAR Validate positive scalar numeric.
    if ~isnumeric(x) || ~isscalar(x) || ~(x >= 0)
        error('The value of %s must be a nonnegative scalar.', argName);
    end
end

function journalList(s)
%JOURNALLIST Validate supported journal/publisher names.
    valid = {'Elsevier','Springer','Nature','Wiley'};
    if ~(ischar(s) || isstring(s)) || ~ismember(string(s), valid)
        error('Specified journal must be one of: %s.', strjoin(valid,', '));
    end
end

%% ========================= UTILITY FUNCTIONS ===========================
function ext = getTextSize(str, rotation, fontName, fontSize, interpreter)
%GETTEXTSIZE Return [w h] of visible text (in points).
%   Computes the glyph bounding box of STR rendered with the given TEXT
%   properties. The result is *visible text* size, not the reserved layout
%   space of a uicontrol or similar.

    % Create an invisible figure
    tmpFig = figure('Visible','off');
    ax = axes('Parent',tmpFig,'Units','points');
    
    % Create a text object
    t  = text(0,0,str, ...
        'Rotation',rotation, 'FontName',fontName, 'FontSize',fontSize, ...
        'Units','points', 'Interpreter',interpreter, 'Parent',ax);
    drawnow; % <<< important for LaTeX
    e = get(t,'Extent');  % [x y w h] in points
    ext = e(3:4);
    close(tmpFig);
    close(tempFig);  % Clean up
end

function [left, right, bottom, top] = measureTextOverhang(str, rotation, fontName, fontSize, interpreter, ha, va)
%MEASURETEXTOVERHANG Visible glyph overhangs around the text anchor.
%   [L,R,B,T] are nonnegative distances (points) from the anchor (0,0) to
%   the axis-aligned glyph bounds, given TEXT alignment and rotation.
%   This measures visible text, not bounding "box" UI padding.

    if nargin < 6 || isempty(ha), ha = 'center'; end
    if nargin < 7 || isempty(va), va = 'top';    end

    % Isolated, consistent environment
    tmpFig = figure('Visible','off','Units','points','Position',[100 100 300 200]);
    ax = axes('Parent',tmpFig,'Units','points','Position',[0 0 300 200]);

    % Place the text at the anchor point (0,0) in POINTS
    t = text(0,0,str,'Parent',ax,'Units','points', ...
        'Rotation',rotation,'HorizontalAlignment',ha,'VerticalAlignment',va, ...
        'FontName',fontName,'FontSize',fontSize,'Interpreter',interpreter);

    drawnow;   % ensure LaTeX/TeX lays out fully

    e = get(t,'Extent'); % [x y w h] in points, axis-aligned box
    %convert to overhangs relative to the anchor point (0,0)
    left   = max(0, -e(1));
    bottom = max(0, -e(2));
    right  = max(0,  e(1) + e(3));
    top    = max(0,  e(2) + e(4));
    close(tmpFig);
end

function v = roundQtr(v)
%ROUNDQTR Round to nearest 0.25 pt.
    v = round(v*4)/4;
end

function [fig, figureWidth, figureHeight] = formatFigure(fig,ppNVArgs)
%FORMATFIGURE Standardize figure object properties and size.
%   Sets figure units to points, hides toolbar if requested, and sets the
%   width based on journal column model and the requested column count.

    % change figure units to points (Illustrator default)
    fig.Units = 'points';
    
    % hide toolbar in figure
    fig.ToolBar = ppNVArgs.ToolBar;
    
    % set figure width to correct artwork sizing
    % Note: this is not necessarily the maximum width allowed, but the width I
    % have found works best in the associated Latex/Word template to fit nicely
    switch string(ppNVArgs.Journal)
        % width for [single, 1.5, double]
        case "Elsevier"
            cols = [252,394,536];
        case "Springer"
            cols = [238,365,493];
        case "Nature"
            cols = [255,382,510];
        case "Wiley"
            cols = [226,352,480];
    end
    switch string(ppNVArgs.Width)
        case "single"
            figureWidth = cols(1);
        case "1.5"
            figureWidth = cols(2);
        case "double"
            figureWidth = cols(3);
        case "custom"
            figureWidth = ppNVArgs.CustomWidth;
            if ppNVArgs.verbosity > 3
                fprintf('Custom Width requested, setting to %f',ppNVArgs.CustomWidth)
            end
    end
    
    % height of figure, defaults to 235 points. Larger height allows more
    % vertical whitespace
    figureHeight = ppNVArgs.Height;
    
    % grab current screen size and place figure at the center
    screenSize = get(0, 'ScreenSize');  % [x y width height] in pixels
    centerX = 0.75*(screenSize(3)/2); % convert to points
    centerY = 0.75*(screenSize(4)/2);
    
    % rescale and shape figure to correct size
    set(fig, 'Units', 'points');
    set(fig, 'Position', [(centerX-figureWidth/2) (centerY-figureHeight/2) ...
        figureWidth figureHeight]);  % position on screen
    
end

function ax = formatAxis(ax,ppNVArgs)
%FORMATAXIS Apply consistent axis/line/legend/text styling.
        
    % set axis properties and color
    set(ax,'FontName',ppNVArgs.FontName,'FontSize',ppNVArgs.AxisFontSize, ...
           'LineWidth',ppNVArgs.AxisLineWidth,'TickDir','out','Box','off', ...
           'XColor',ppNVArgs.AxisColor,'YColor',ppNVArgs.AxisColor,'ZColor',ppNVArgs.AxisColor);
    % change background color
    ax.Color = ppNVArgs.BoxColor;

    % set grid properties
    % turn on or off
    grid(ax, ppNVArgs.Grid);
    if ppNVArgs.Grid
        ax.GridLineWidth = ppNVArgs.GridLineWidth;
        ax.GridColor     = ppNVArgs.GridColor;
    end 
    
    % scale labels and fonts to 1.0 of text size
    ax.LabelFontSizeMultiplier = 1.0;
    ax.TitleFontSizeMultiplier = 1.0;
    
    % set all line stroke properties and marker sizes
    lines = findall(ax, 'Type', 'Line');
    for i = 1 : numel(lines)
        lines(i).LineWidth  = ppNVArgs.LineWidth;
        lines(i).MarkerSize = ppNVArgs.MarkerSize;
    end
    
    % set legend properties to font properties
    leg = findobj(gcf, 'Type', 'Legend');
    for i = 1 : length(leg)
        set(leg(i), 'FontSize', ppNVArgs.FontSize, 'Box', 'off',...
                    'Interpreter', ppNVArgs.Interpreter);
    end
    
    % set title font properties
    set(ax.Title, 'FontSize', ppNVArgs.FontSize, ...
        'FontWeight', ax.Title.FontWeight, 'Interpreter', ppNVArgs.Interpreter,...
        'Color', ppNVArgs.FontColor,'Parent',ax);
    
    % set axis label and tick font properties
    set(ax.XLabel, 'FontSize', ppNVArgs.AxisFontSize, 'FontWeight', ax.XLabel.FontWeight,...
        'Interpreter', ppNVArgs.Interpreter, 'Color', ppNVArgs.FontColor);
    set(ax.YLabel, 'FontSize', ppNVArgs.AxisFontSize, 'FontWeight', ax.YLabel.FontWeight,...
        'Interpreter', ppNVArgs.Interpreter, 'Color', ppNVArgs.FontColor);
end

function [ax, axisOffset] = reduceWhitespace(ax,ppNVArgs,subW,subH,NVArgs)
%REDUCEWHITESPACE Trim margins, redraw ticks/labels in point units, and place title.
%
%   Handles linear/log axes via mapDataToPoints, keeps labels aligned in a
%   consistent offset from the axes origin (not from tick extents), and can
%   lift titles above the axes if requested or if data overlaps the title box.

    arguments
        ax (1,1) matlab.graphics.axis.Axes
        ppNVArgs (1,1) struct
        subW (1,1) double {numericPositiveScalar}
        subH (1,1) double {numericPositiveScalar}
        NVArgs.xOffset (1,1) double {numericPositiveScalar} = 0
        NVArgs.yOffset (1,1) double {numericPositiveScalar} = 0
        NVArgs.customOffset double = []
        NVArgs.moveTitle logical = false
    end
    
    % manually reduce whitespace
    % only works for non-tiled layouts
    if isInTiledLayout(ax)
        % Defer sizing to tiled layout; still tighten padding
        tl = ancestor(ax,'matlab.graphics.layout.TiledChartLayout','toplevel');
        if isprop(tl,'TileSpacing'), tl.TileSpacing = 'compact'; end
        if isprop(tl,'Padding'),     tl.Padding     = 'compact'; end
        % Tiled Axes: CAN NOT set Position/InnerPosition/OuterPosition
        % The layout owns Position/InnerPosition/OuterPosition
        % Let the layout arrange the axes rectangle.
        drawnow;
        return
    end
        
    Spacing = ppNVArgs.SpacingOffset;
    xOffset = NVArgs.xOffset;
    yOffset = NVArgs.yOffset;
    tickLength = ppNVArgs.TickLength(1); % desired tick length
    
    % record axis limits and labels to preserve their original values
    XLim = ax.XLim;
    YLim = ax.YLim;
    XTick = ax.XTick;
    YTick = ax.YTick;
    XTickLabels = ax.XTickLabel;
    YTickLabels = ax.YTickLabel;
    XTickRotation = ax.XTickLabelRotation;
    YTickRotation = ax.YTickLabelRotation;

    if (abs(XTickRotation) > 180) || (abs(YTickRotation) > 180)
        warning('Warning: Tick Rotations > 180 or < -180 detected. Behavior may be unexpected.')
    end
    
    % Compute offsets if not supplied
    if isempty(NVArgs.customOffset)
        % determine offsets to avoid clipping text
        [leftOffset, rightOffset, bottomOffset, topOffset, NVArgs.moveTitle] = getTextExtents(ax, ppNVArgs, subW, subH);
    else
        assert(numel(NVArgs.customOffset)==4,'customOffset must be of the form [leftOffset, rightOffset, bottomOffset, topOffset].');
        leftOffset   = NVArgs.customOffset(1);
        rightOffset  = NVArgs.customOffset(2);
        bottomOffset = NVArgs.customOffset(3);
        topOffset    = NVArgs.customOffset(4);
    end

    % find x and y label extents
    
    % Flags for presence of labels
    hasXLabel = ~isempty(ax.XLabel.String);
    hasYLabel = ~isempty(ax.YLabel.String);

    % Plot area width and height (points)
    plotW = subW - (leftOffset + rightOffset);
    plotH = subH - (bottomOffset + topOffset);

    % Place axes rectangle
    ax.Position = [leftOffset + xOffset, bottomOffset + yOffset, plotW, plotH];

    % Restore limits/ticks
    ax.XLim = XLim;
    ax.YLim = YLim;
    ax.XTick = XTick;
    ax.YTick = YTick;
    
    % Manually place tick labels
    
    % Size-aware tick length (normalize by longer axis)
    % After setting axix Position, before placing any text:
    drawnow;  % ensure axis Position is final for accurate point mapping
    axisW = ax.Position(3);
    axisH = ax.Position(4);
    
    % adjust tick length: use actual axis size so the ratio is exact
    ax.TickLength = [tickLength/max([axisW,axisH]), tickLength/max([axisW,axisH])];

    % Functions to map x/y data point to position in points (linear/log)
    XDataToPt = @(x) mapDataToPoints(ax.XScale, XLim, x, axisW);
    YDataToPt = @(y) mapDataToPoints(ax.YScale, YLim, y, axisH);

    % Remove exisiting tick labels and redraw as text so we can position in
    % point space
    ax.XTickLabel = {};  ax.YTickLabel = {};

    % X ticks
    for i = 1 : numel(XTick)
        xLoc = XDataToPt(XTick(i));
        yLoc = -(tickLength + 2*Spacing); % axes-local, just below axis line
        [vertAlign,horizAlign] = tickAlignForRotation(XTickRotation,'x');

        text(xLoc, yLoc, XTickLabels{i}, ...
            'Parent', ax, 'Units','points', ...
            'Rotation', XTickRotation, ...
            'HorizontalAlignment', horizAlign, ...
            'VerticalAlignment',   vertAlign, ...
            'FontSize', ax.FontSize, 'FontName', ax.FontName, ...
            'Color', ppNVArgs.FontColor, ...
            'Clipping','off');    % allow negative y
    end

    % Y ticks
    for i = 1 : numel(YTick)
        xLoc = -(2*Spacing + tickLength); % to the left of axis line
        yLoc = YDataToPt(YTick(i));
        [vertAlign,horizAlign] = tickAlignForRotation(YTickRotation,'y');

        text(xLoc, yLoc, YTickLabels{i}, ...
            'Parent', ax, 'Units','points', ...
            'Rotation', YTickRotation, ...
            'HorizontalAlignment', horizAlign, ...
            'VerticalAlignment',   vertAlign, ...
            'FontSize', ax.FontSize, 'FontName', ax.FontName, ...
            'Color', ppNVArgs.FontColor, ...
            'Clipping','off');
    end

    % Manually place axis labels (constant offsets from origin so subplots align)

    % X label
    if hasXLabel
        text(axisW/2, -(bottomOffset - 1*Spacing), ax.XLabel.String, ...
            'Parent', ax, 'Units','points', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', ax.FontSize,'FontName',ax.FontName,...
            'Color', ppNVArgs.FontColor,'Interpreter',ppNVArgs.Interpreter);
        ax.XLabel.String = '';
    end

    % Y label
    if hasYLabel
        text(-(leftOffset + 1*Spacing), axisH/2, ax.YLabel.String, ...
            'Parent', ax, 'Units','points', ...
            'HorizontalAlignment','center', 'VerticalAlignment','top', ...
            'Rotation',90, ...
            'FontSize', ax.FontSize,'FontName',ax.FontName,...
            'Color', ppNVArgs.FontColor,'Interpreter',ppNVArgs.Interpreter);
        ax.YLabel.String = '';
    end

    % Manually place title (option to lift above axis line)
    if ~isempty(ax.Title.String)
        putAbove = ppNVArgs.MoveTitleAboveAxis || any(NVArgs.moveTitle);
        va = ternary(putAbove,'bottom','middle');
        text(axisW/2, axisH, ax.Title.String, ...
            'Parent', ax, 'Units','points', ...
            'HorizontalAlignment','center','VerticalAlignment',va, ...
            'FontSize', ax.FontSize,'FontName',ax.FontName, ...
            'Color', ppNVArgs.FontColor,'Interpreter',ppNVArgs.Interpreter);
        ax.Title.String = '';
    end

    axisOffset = [leftOffset, rightOffset, bottomOffset, topOffset];

end

function [leftOffset, rightOffset, bottomOffset, topOffset, moveTitle] = getTextExtents(ax,ppNVArgs,subW,subH)
%GETTEXTEXTENTS Compute margins required to avoid clipping ticks/labels.
%
%   Returns LEFT/RIGHT/BOTTOM/TOP offsets in points that guarantee room for
%   the longest tick label, axis labels, tick marks, and (optionally) any
%   title. Also checks whether plotted data overlaps the title glyph box,
%   suggesting the title be moved above the axes.

    Spacing = ppNVArgs.SpacingOffset;
    tickLength = ppNVArgs.TickLength(1);
    XTickRotation = ax.XTickLabelRotation;
    YTickRotation = ax.YTickLabelRotation;

    % Determine tick alignment (affects overhang)
    [xVA,xHA] = tickAlignForRotation(XTickRotation,'x');
    [yVA,yHA] = tickAlignForRotation(YTickRotation,'y');

    % For the longest X tick string, find the height for calc. extents
    if ~isempty(ax.XTick)
        Lx = cellfun(@strlength, ax.XTickLabel); [~,ix] = max(Lx);
        [~,~,b,t] = measureTextOverhang(ax.XTickLabel{ix}, XTickRotation, ppNVArgs.FontName, ppNVArgs.AxisFontSize, ppNVArgs.Interpreter, xHA, xVA);
        xTickHeight = b + t;
    else
        xTickHeight = 0;
    end

    % For the longest Y tick string, find the width for calc. extents
    if ~isempty(ax.YTick)
        Ly = cellfun(@strlength, ax.YTickLabel); [~,iy] = max(Ly);
        [l,r,~,~] = measureTextOverhang(ax.YTickLabel{iy}, YTickRotation, ppNVArgs.FontName, ppNVArgs.AxisFontSize, ppNVArgs.Interpreter, yHA, yVA);
        yTickWidth = l + r;
    else
        yTickWidth = 0;
    end
    
    % Find extent of the text box for last tick labels (the ones which may fall
    % outside of the figure, i.e. at the X and Y limits)
    if ~isempty(ax.XTick)
        [xEdgeHalfL,xEdgeHalfR,b,t] = measureTextOverhang(ax.XTickLabel{end},XTickRotation,ppNVArgs.FontName,ppNVArgs.AxisFontSize,ppNVArgs.Interpreter,xHA,xVA);
        % horizontal half-extent of the last tick label (already includes rotation)
    else
        xEdgeHalfR = 0;
    end
    if ~isempty(ax.YTick)
        [l,r,yEdgeHalfB,yEdgeHalfT] = measureTextOverhang(ax.YTickLabel{end},YTickRotation,ppNVArgs.FontName,ppNVArgs.AxisFontSize,ppNVArgs.Interpreter,yHA,yVA);
        % vertical half-extent of the last tick label (already includes rotation)
    else
        yEdgeHalfT = 0;
    end
    
    % Find X and Y axis label extents
    
    % check label limits if xLabels and yLabels contain string
    if ~isempty(ax.XLabel.String)
        set(ax.XLabel,'Units','points'); ex = get(ax.XLabel,'Extent'); xLabelSize = ex(3:4);
    else
        xLabelSize = [0 0];
    end
    if ~isempty(ax.YLabel.String)
        set(ax.YLabel,'Units','points'); ey = get(ax.YLabel,'Extent'); yLabelSize = ey(3:4);
    else
        yLabelSize = [0 0];
    end

    % add in padding for final tick  marks
    padding = max(Spacing,0.5);   % 0.75–1 pt works well visually
    
    % ---------------------------- LEFT OFFSET ----------------------------
    % distance from the left edge to the origin
    if yLabelSize(1) > 1
        leftOffset = 2*Spacing + yLabelSize(1) + yTickWidth + tickLength;
    else
        leftOffset = 3*Spacing + yLabelSize(1) + yTickWidth + tickLength;
    end

    % Ensure first X tick label fits
    % minimum left margin
    minLeft = Spacing;

    % check the first tick value
    if ~isempty(ax.XTick)
        xFirst = ax.XTick(1);
    else
        xFirst = ax.XLim(1);  % fall back to axis max
    end

    % map last tick to plot coordinates (points) from the LEFT edge
    xCenter = leftOffset + ...
              mapDataToPoints(ax.XScale, ax.XLim, xFirst, subW - leftOffset - minLeft);

    % determine how far the right edge of the label extends beyond the axis
    overhangL = xCenter - xEdgeHalfL;  

    if overhangL < 0
        % reserve extra on top of leftOffset so the label fits exactly
        leftOffset = leftOffset + abs(ceil(4*(overhangL))/4);
    end
    
    % --------------------------- RIGHT OFFSET ---------------------------
    % Ensure last X tick label fits
    % minimum right margin
    minRight = Spacing;

    % width available for the plot if only reserve minRight for now
    plotWminR = max(1, subW - leftOffset - minRight);

    % check the last tick value
    if ~isempty(ax.XTick)
        xLast = ax.XTick(end);
    else
        xLast = ax.XLim(2);  % fall back to axis max
    end

    % map last tick to plot coordinates (points) from the LEFT edge
    xCenter = leftOffset + ...
              mapDataToPoints(ax.XScale, ax.XLim, xLast, plotWminR);

    % determine how far the right edge of the label extends beyond the axis
    overhangR = (xCenter + xEdgeHalfR) - subW;  

    if overhangR > 0
        % reserve extra on top of minRight so the label fits exactly
        rightOffset = ceil(4*(minRight + overhangR + padding))/4;
    elseif ppNVArgs.AxisLineWidth > minRight
        % ensure the axes path itself doesn't touch the border
        rightOffset = ppNVArgs.AxisLineWidth + padding;
    else
        rightOffset = minRight + padding;
    end
    
    % --------------------------- BOTTOM OFFSET ---------------------------
    % distance from the bottom edge to the origin
    bottomOffset = (1*Spacing + xLabelSize(2) + xTickHeight + tickLength);

    % Ensure first Y tick label fits
    % minimum top margin
    minBottom = Spacing;

    % check the first tick value
    if ~isempty(ax.YTick)
        yFirst = ax.YTick(1);
    else
        yFirst = ax.YLim(2);
    end

    % map last tick to plot coordinates (points) from the BOTTOM edge
    yCenter = bottomOffset + ...
              mapDataToPoints(ax.YScale, ax.YLim, yFirst, subH - bottomOffset - minBottom);
    
    % determine how far the top edge of the label extends beyond the axis
    overhangB = yCenter - yEdgeHalfB;

    if overhangB < 0
        % reserve extra on top of leftOffset so the label fits exactly
        bottomOffset = bottomOffset + abs(ceil(4*(overhangB))/4);
    end

    % ---------------------------- TOP OFFSET ----------------------------
    % Ensure last Y tick label fits
    % minimum top margin
    minTop = Spacing;

    % height available for the plot if only reserve minTop for now
    plotHminT = max(1, subH - bottomOffset - minTop);

    % check the last tick value
    if ~isempty(ax.YTick)
        yLast = ax.YTick(end);
    else
        yLast = ax.YLim(2);
    end

    % map last tick to plot coordinates (points) from the BOTTOM edge
    yCenter = bottomOffset + ...
              mapDataToPoints(ax.YScale, ax.YLim, yLast, plotHminT);
    
    % determine how far the top edge of the label extends beyond the axis
    overhangT = (yCenter + yEdgeHalfT) - subH;

    if overhangT > 0
        % reserve extra on top of minRight so the label fits exactly
        topOffset = ceil(4*(minTop + overhangT  + padding))/4;
    elseif ppNVArgs.AxisLineWidth > minTop
        % ensure the axes path itself doesn't touch the border
        topOffset = ppNVArgs.AxisLineWidth + padding;
    else
        topOffset = minTop + padding;
    end

    % --------------------------- TITLE OFFSET ---------------------------
    % check for overlapping data on title
    % flag if data encroaches on title box
    moveTitle = false;

    if ~isempty(ax.Title.String)

        % Size of the inner axes (in points; origin = bottom-left of axes)
        axisW = subW  - leftOffset - rightOffset;
        axisH = subH - bottomOffset - topOffset;

        % measure visible title glyph box for a centered, middle-anchored title
        % (assumed)
        [L,R,B,T] = measureTextOverhang(ax.Title.String, ax.Title.Rotation, ...
            ppNVArgs.FontName, ppNVArgs.AxisFontSize, ...
            ppNVArgs.Interpreter, 'center','middle');

        % title bounds in points
        titleXExt = axisW/2 + [-L, R];         % [xmin xmax] in points
        titleYExt = axisH   + [-B, T];         % [ymin ymax] in points

        % Scan visible children with numeric X/Y data; convert to axes-points
        kids = ax.Children;
        for ii = 1:numel(kids)
            k = kids(ii);
            if strcmpi(k.Visible,'on') && isprop(k,'XData') && isprop(k,'YData')
                xdat = k.XData; ydat = k.YData;
                if ~isempty(xdat) && ~isempty(ydat)
                    % make vectors the same shape
                    try
                        xdat = xdat(:); ydat = ydat(:);
                    catch, continue
                    end

                    % map to physcial points (handles linear/log and nonpositive log values)
                    xpt = mapDataToPoints(ax.XScale, ax.XLim, xdat, axisW);
                    ypt = mapDataToPoints(ax.YScale, ax.YLim, ydat, axisH);

                    % ignore NaN/Inf from mapping
                    m = isfinite(xpt) & isfinite(ypt);

                    % is sample inside the title rectangle?
                    hit =  xpt(m) >= titleXExt(1) & xpt(m) <= titleXExt(2) & ...
                        ypt(m) >= titleYExt(1) & ypt(m) <= titleYExt(2);

                    if any(hit)
                        if ppNVArgs.verbosity > 2
                        fprintf(['Title overlaps plotted data. Increasing' ...
                            ' top offset and moving the title.\n']);
                        end
                        moveTitle = true;

                    end
                end
            end
        end

    end

    if moveTitle
        topOffset = max(B+T,topOffset);
    end

    % Round to quarter points for stable Illustrator alignment
    leftOffset   = roundQtr(leftOffset);
    rightOffset  = roundQtr(rightOffset);
    bottomOffset = roundQtr(bottomOffset);
    topOffset    = roundQtr(topOffset);

end

function tf = isInTiledLayout(ax)
%ISINTILEDLAYOUT True if AX belongs to a tiledlayout.
    tf = ~isempty(ancestor(ax,'matlab.graphics.layout.TiledChartLayout','toplevel'));
end

function pt = mapDataToPoints(scale, lim, val, axisLen)
%MAPDATATOPOINTS Map data to axis-local point coordinates.
%   PT = MAPDATATOPOINTS(scale,lim,val,axisLen) converts val in data units
%   into points along an axis of length axisLen, honoring linear or log scale.

    if strcmpi(scale,'log')
        % guard against nonpositive values for log axes
        v  = max(val,  realmin('double'));
        l1 = max(lim(1),realmin('double'));
        l2 = max(lim(2),realmin('double'));
        pt = axisLen * (log10(v) - log10(l1)) / (log10(l2) - log10(l1));
    else
        pt = axisLen * (val - lim(1)) / (lim(2) - lim(1));
    end
end

function [va,ha] = tickAlignForRotation(rot, axis)
%TICKALIGNFORROTATION Choose text alignments for rotated tick labels.
%   Returns VerticalAlignment VA and HorizontalAlignment HA.

    if axis == 'x'
        if rot >= 45
            if rot >= 135,  va = 'bottom'; ha = 'center';
            else           va = 'middle'; ha = 'right';
            end
        elseif rot <= -45
            if rot <= -135, va = 'bottom'; ha = 'center';
            else           va = 'middle'; ha = 'left';
            end
        else
            va = 'top'; ha = 'center';
        end
    else % 'y'
        if rot >= 45
            if rot >= 135,  ha = 'left';   va = 'middle';
            else           ha = 'center'; va = 'bottom';
            end
        elseif rot <= -45
            if rot <= -135, ha = 'left';   va = 'middle';
            else           ha = 'center'; va = 'top';
            end
        else
            ha = 'right'; va = 'middle';
        end
    end
end

function out = ternary(cond, a, b)
%TERNARY Simple ternary helper (returns A if COND, else B).
    if cond, out = a; else, out = b; end
end