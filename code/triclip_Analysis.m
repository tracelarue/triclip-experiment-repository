clc
clear all
close all


% ==========================================================================

% Get the full path of this script
scriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(scriptPath);
projectRoot = fileparts(scriptDir);

% Build data paths
dataDir = fullfile(projectRoot, 'Data');
testIdxDir = fullfile(dataDir, 'testidx');

% -------------------------- Load stable region index data for each test
load(fullfile(testIdxDir, 'test1idx.mat'));
load(fullfile(testIdxDir, 'test2idx.mat'));
load(fullfile(testIdxDir, 'test3idx.mat'));
load(fullfile(testIdxDir, 'test4idx.mat'));
load(fullfile(testIdxDir, 'test5idx.mat'));
load(fullfile(testIdxDir, 'test6idx.mat'));
load(fullfile(testIdxDir, 'test7idx.mat'));


% -------------------------- Load and preprocess data

% Define paths to data files using relative paths
folderPaths = {
    fullfile(dataDir, '05_21_25'),
    fullfile(dataDir, '06_13_25'),
    fullfile(dataDir, '07_09_25'),
    fullfile(dataDir, '07_16_25'),
    fullfile(dataDir, '07_23_25'),
    fullfile(dataDir, '08_05_25'),
    fullfile(dataDir, '08_07_25')
};

% Initialize cell arrays to store data
number_of_tests = length(folderPaths);
all_test_data = cell(number_of_tests, 1);

disp("Reading CSV files...")
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

% Load all CSV files
for test = 1:number_of_tests
    % Get file listing
    filePattern = fullfile(folderPaths{test}, '*.csv');
    csvFiles = dir(filePattern);
    
    % Initialize storage for this test
    all_test_data{test} = cell(length(csvFiles), 1);
    
    % Read each file
    for subtest = 1:length(csvFiles)
        baseFileName = csvFiles(subtest).name;
        fullFileName = fullfile(folderPaths{test}, baseFileName);
        
        % Read CSV file into table
        data = readtable(fullFileName);
        all_test_data{test}{subtest} = data;
    end
end

disp("CSV files read. Preprocessing data...")

% Unpack for easier access (up to number_of_tests)
for i = 1:number_of_tests
    eval(sprintf('test%idata = all_test_data{i};', i));
end

% -------------------------- Apply calibration to force measurements

% Calibration slopes for each pin
force_calibration = [
    -1621.4975*3.3  % pin1
    -1291.6308*3.3  % pin2
    -1368.1537*3.3  % pin3
    -1017.2041*3.3  % pin4
    -1106.5112*3.3  % pin5
    -1089.2077*3.3  % pin6
    -2948.5710*3.3  % pin7
    -1427.9047*3.3  % pin8
];

% Process each test dataset
for test = 1:number_of_tests
    for subtest = 1:length(all_test_data{test})
        % Remove the extra "." from Time column
        all_test_data{test}{subtest}.Time = strrep(all_test_data{test}{subtest}.Time, '..', '.');
        
        % Apply calibration to each pin and zero them
        for pin = 1:8
            force_x = ['Force' num2str(pin)];
            % Apply calibration
            all_test_data{test}{subtest}.(force_x) = all_test_data{test}{subtest}.(force_x) * force_calibration(pin);
            % Zero the force by subtracting the initial value
            initial_force = all_test_data{test}{subtest}.(force_x)(1);
            all_test_data{test}{subtest}.(force_x) = all_test_data{test}{subtest}.(force_x) - initial_force;
        end
        
        % Convert and zero flow rate
        all_test_data{test}{subtest}.FlowRate_ml_s_ = all_test_data{test}{subtest}.FlowRate_ml_s_ * 200.44; %% Apply flow rate calibration (200.4 ml/s/V)
        all_test_data{test}{subtest}.FlowRate_ml_s_ = all_test_data{test}{subtest}.FlowRate_ml_s_ - all_test_data{test}{subtest}.FlowRate_ml_s_(1);
    end
end

% Reassign to original variables (up to number_of_tests)
for i = 1:number_of_tests
    eval(sprintf('test%idata = all_test_data{i};', i));
end

% -------------------------- Define experimental groups for analysis

disp('Data loaded and preprocessed')

% -------------------- Calculate average data for each test

% Initialize structures for averages (up to number_of_tests)
for i = 1:number_of_tests
    eval(sprintf('test%davg = struct();', i));
end

% Calculate averages for each test (up to number_of_tests)
for test = 1:number_of_tests
    test_x = eval(['test' num2str(test) 'data']);
    test_x_index = eval(['test' num2str(test) 'idx']);
    test_x_avg = eval(['test' num2str(test) 'avg']);
    
    for subtest = 1:length(test_x)
        idx = test_x_index(subtest);
        start = idx.segmentStartIdx;
        stop = idx.segmentEndIdx;
        
        % Add test number as first field
        test_x_avg(subtest).Test = test;
        
        % Average pressure and flow
        test_x_avg(subtest).Pressure = mean(test_x{subtest}.Pressure(start:stop));
        test_x_avg(subtest).FlowRate = mean(test_x{subtest}.FlowRate_ml_s_(start:stop));
        
        % Average force for each pin
        for pin = 1:8
            force_x = ['Force' num2str(pin)];
            test_x_avg(subtest).(force_x) = mean(test_x{subtest}.(force_x)(start:stop));
        end
    end
    
    % Store back to the original variable
    eval(['test' num2str(test) 'avg = test_x_avg;']);
end

disp('Average data calculation complete')

% Data Indices based on the order of the files in the folders
%        13:15  16:18  19:21  22:24       
% Test 1 - AS,  AP,     ASAP        AP Data is bad
% Test 2 - SP,  AP,     SPAP 
% Test 3 - AS,  SP,     SPAS
% Test 4 - AS,  AP,     ASAP
% Test 5 - AS,  AP,     ASAP
% Test 6 - SP,  SPAP,   SPAS,   AP
% Test 7 - SP,  SPAP,   SPAS,   AP

healthy_data = [test1avg(1:9), test2avg(1:9), test3avg(1:9), test4avg(1:9), test5avg(1:9), test6avg(1:9), test7avg(1:9) ];
diseased_data = [test1avg(10:12), test2avg(10:12), test3avg(10:12), test4avg(10:12), test5avg(10:12), test6avg(10:12), test7avg(10:12)];
AS_data= [test1avg(13:15), test3avg(13:15), test4avg(13:15), test5avg(13:15)];
AP_data= [test2avg(16:18), test4avg(16:18), test5avg(16:18), test6avg(22:24), test7avg(22:24)];
SP_data= [test2avg(13:15), test3avg(16:18), test6avg(13:15), test7avg(13:15)];
ASAP_data= [test1avg(19:21), test4avg(19:21), test5avg(19:21)];
SPAS_data= [test3avg(19:21), test6avg(19:21), test7avg(19:21)];
SPAP_data= [test2avg(19:21), test6avg(16:18), test7avg(16:18)];

%% Flow Rate Bar Graph by Intervention Type
% Hex to RGB conversion function
hex2rgb = @(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]) / 255;

% ========================== COLOR CONTROL SECTION ==========================
% Define all colors for consistent plotting across the entire analysis
% Modify these hex codes to change colors throughout all plots

% Intervention Colors (used in all bar graphs and contour plots)
COLORS = struct();
COLORS.Diseased = hex2rgb('#EDE6CC');   % Light cream/beige
COLORS.AS       = hex2rgb('#FCCC73');   % Light golden yellow
COLORS.AP       = hex2rgb('#FC9933');   % Warm orange
COLORS.SP       = hex2rgb('#E65940');   % Red-orange
COLORS.ASAP     = hex2rgb('#CC1A1A');   % Burgundy/wine
COLORS.SPAS     = hex2rgb('#801A33');   % Wine
COLORS.SPAP     = hex2rgb('#333359');   % Dark teal/navy

% Bar graph color arrays (all interventions including diseased)
bar_colors = [ ...
    COLORS.Diseased; ...
    COLORS.AS; ...
    COLORS.AP; ...
    COLORS.SP; ...
    COLORS.ASAP; ...
    COLORS.SPAS; ...
    COLORS.SPAP ...
];

% Force plot color arrays (interventions only, no diseased)
force_colors = [ ...
    COLORS.AS; ...
    COLORS.AP; ...
    COLORS.SP; ...
    COLORS.ASAP; ...
    COLORS.SPAS; ...
    COLORS.SPAP ...
];

% Rectangle colors for contour plots (direct color references)
AS_rect_color = COLORS.AS;
AP_rect_color = COLORS.AP;
SP_rect_color = COLORS.SP;


% Define interventions and their corresponding data
interventions = {'Dis      ', 'AS     ', 'AP     ', 'SP     ', 'ASAP        ', 'SPAS        ', 'SPAP        '};
intervention_data = {diseased_data, AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

% Which tests participate in each intervention (cell array, 1 row per intervention)
test_participation = {
    [1:7];          % Diseased
    [1,3,4,5];      % AS
    [2,4,5,6,7];    % AP
    [2,3,6,7];      % SP
    [1,4,5];        % ASAP
    [3,6,7];        % SPAS
    [2,6,7];        % SPAP
};

% Preallocate arrays for means and stds
num_interv = numel(interventions);
avg_flow = zeros(1, num_interv);
std_flow = zeros(1, num_interv);

% Compute mean and std for each intervention (per-test means -> std across tests)
for i = 1:num_interv
    data = intervention_data{i};
    if isempty(data)
        avg_flow(i) = NaN;
        std_flow(i) = NaN;
        continue;
    end
    tests = unique([data.Test]);                 % unique test IDs present
    perTestMeans = zeros(1, numel(tests));
    for k = 1:numel(tests)
        t = tests(k);
        mask = [data.Test] == t;
        perTestMeans(k) = mean([data(mask).FlowRate]);  % average for this test
    end
    avg_flow(i) = mean(perTestMeans);            % mean across tests
    if numel(perTestMeans) > 1
        std_flow(i) = std(perTestMeans);         % std across tests
    else
        std_flow(i) = 0;
    end
end

% Plot setup
figure('Position', [100, 100, 1200, 700]);

%{
Old color definitions - now using centralized colors from top of script
bar_colors = [ ...
    0.2 0.2 0.2;        % Diseased - dark gray
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];
%}

x = 1:num_interv;
b = bar(x, avg_flow, 'FaceColor', 'flat');
for i = 1:num_interv
    b.CData(i,:) = bar_colors(i,:);
end
hold on;

% Error bars
errorbar(x, avg_flow, std_flow, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

% Add significance lines for comparisons to diseased state
% P-values from R analysis: AS=0.0138(*), AP=0.0187(*), SP=0.1607(ns), ASAP=0.0002(***), SPAS=0.0977(ns), SPAP=0.2717(ns)
sig_interventions = [2, 3, 5]; % AS, AP, ASAP (indices in interventions array)
sig_symbols = {'*', '*', '***'}; % Significance symbols
sig_pvals = [0.0138, 0.0187, 0.0002]; % P-values

% Calculate height for significance lines (above error bars)
max_bar_height = max(avg_flow + std_flow);
line_height_base = max_bar_height + 2; % Base height for first significance line
line_spacing = 4; % Vertical spacing between multiple significance lines

for i = 1:length(sig_interventions)
    intervention_idx = sig_interventions(i);
    diseased_idx = 1; % Diseased is always first in the array
    
    % Calculate line height (stagger multiple lines)
    line_height = line_height_base + (i-1) * line_spacing;
    
    % Draw horizontal line connecting diseased bar to intervention bar
    plot([diseased_idx, intervention_idx], [line_height, line_height], 'k-', 'LineWidth', 1);
    
    % Draw vertical ticks at each end
    tick_height = 0.8;
    plot([diseased_idx, diseased_idx], [line_height - tick_height/2, line_height + tick_height/2], 'k-', 'LineWidth', 1);
    plot([intervention_idx, intervention_idx], [line_height - tick_height/2, line_height + tick_height/2], 'k-', 'LineWidth', 1);
    
    % Add significance symbol at midpoint
    midpoint_x = (diseased_idx + intervention_idx) / 2;
    text(midpoint_x, line_height - 2, sig_symbols{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
end

% Marker/legend setup
marker_types = {'o', 's', 'd', '^', 'v', '>', '<'};
marker_colors = lines(7);
h_leg = gobjects(1, 7);

%{ 
% Plot individual test markers for each intervention
for i = 1:num_interv
    tests = test_participation{i};
    n_tests = numel(tests);
    for j = 1:n_tests
        t = tests(j);
        % Find all data points for this test in this intervention
        mask = [intervention_data{i}.Test] == t;
        test_flows = [intervention_data{i}(mask).FlowRate];
        if isempty(test_flows), continue; end
        test_mean = mean(test_flows);
        x_offset = (j - (n_tests+1)/2) * 0.08;
        scatter(i + x_offset, test_mean, 80, marker_types{t}, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', marker_colors(t,:), 'LineWidth', 1.2);
    end
end
%}
 
% Legend for test markers
for t = 1:7
    h_leg(t) = plot(NaN, NaN, marker_types{t}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', marker_colors(t,:), 'LineWidth', 1.2);
end
%legend(h_leg, arrayfun(@(x) sprintf('Test %d', x), 1:7, 'UniformOutput', false), ...
%    'Location', 'best', 'NumColumns', 2);

% Axis labels and formatting
title('', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Intervention', 'FontSize', 12);
ylabel('Flow Rate (ml/s)', 'FontSize', 12);
set(gca, 'XTick', x, 'XTickLabel', interventions, 'XTickLabelRotation',30);

% Adjust y-axis limits to accommodate significance lines
current_ylim = ylim;
new_ylim = [current_ylim(1), max(current_ylim(2), line_height_base + (length(sig_interventions)-1) * line_spacing + 4)];
ylim(new_ylim);

set(gcf, 'Color', 'white');
grid on;
hold off;

pubPlot('SpacingOffset',1,'Filename','Flow_vs_Intervention','FileExtension',{'.png','.pdf'});


%% Pressure Bar Graph by Intervention Type

% Define interventions and their corresponding data
interventions = {'Diseased              ', 'AS     ', 'AP     ', 'SP     ', 'ASAP        ', 'SPAS        ', 'SPAP        '};
intervention_data = {diseased_data, AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

% Which tests participate in each intervention (cell array, 1 row per intervention)
test_participation = {
    [1:7];      % Diseased: all tests
    [1,3,4,5];
    [2,4,5,6,7];
    [2,3,6,7];
    [1,4,5];
    [3,6,7];
    [2,6,7]
};

% Preallocate arrays for means and stds
num_interv = numel(interventions);
avg_pressure = zeros(1, num_interv);
std_pressure = zeros(1, num_interv);

% Compute mean and std for each intervention (per-test means -> std across tests)
for i = 1:num_interv
    data = intervention_data{i};
    if isempty(data)
        avg_pressure(i) = NaN;
        std_pressure(i) = NaN;
        continue;
    end
    tests = unique([data.Test]);
    perTestMeans = zeros(1, numel(tests));
    for k = 1:numel(tests)
        t = tests(k);
        mask = [data.Test] == t;
        perTestMeans(k) = mean([data(mask).Pressure]);
    end
    avg_pressure(i) = mean(perTestMeans);
    if numel(perTestMeans) > 1
        std_pressure(i) = std(perTestMeans);
    else
        std_pressure(i) = 0;
    end
end

% Plot setup
figure('Position', [100, 100, 1200, 700]);

%{
Old color definitions - now using centralized colors from top of script
bar_colors = [ ...
    0.2 0.2 0.2;        % Diseased - dark gray
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];
%}

x = 1:num_interv;
b = bar(x, avg_pressure, 'FaceColor', 'flat');
for i = 1:num_interv
    b.CData(i,:) = bar_colors(i,:);
end
hold on;

% Error bars
errorbar(x, avg_pressure, std_pressure, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

% Marker/legend setup
marker_types = {'o', 's', 'd', '^', 'v', '>', '<'};
marker_colors = lines(7);
h_leg = gobjects(1, 7);

%{
% Plot individual test markers for each intervention
for i = 1:num_interv
    tests = test_participation{i};
    n_tests = numel(tests);
    for j = 1:n_tests
        t = tests(j);
        % Find all data points for this test in this intervention
        mask = [intervention_data{i}.Test] == t;
        test_pressures = [intervention_data{i}(mask).Pressure];
        if isempty(test_pressures), continue; end
        test_mean = mean(test_pressures);
        x_offset = (j - (n_tests+1)/2) * 0.08;
        scatter(i + x_offset, test_mean, 80, marker_types{t}, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', marker_colors(t,:), 'LineWidth', 1.2);
    end
end
%}

% Legend for test markers
for t = 1:7
    h_leg(t) = plot(NaN, NaN, marker_types{t}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', marker_colors(t,:), 'LineWidth', 1.2);
end
%legend(h_leg, arrayfun(@(x) sprintf('Test %d', x), 1:7, 'UniformOutput', false), ...
%    'Location', 'best', 'NumColumns', 2);

% Axis labels and formatting
% title('Average Pressures by Intervention Type with Individual Test Values', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Intervention Type', 'FontSize', 12);
ylabel('Pressure (mmHg)', 'FontSize', 12);
set(gca, 'XTick', x, 'XTickLabel', interventions, 'XTickLabelRotation',30);
set(gcf, 'Color', 'white');
grid on;
hold off;

pubPlot('Filename','Pressure_vs_Intervention','FileExtension',{'.png','.pdf'});

%% Force Difference Bar Graphs by Intervention Type (6 Subplots)

% List of intervention types to plot (excluding 'Diseased')
interventions = {'AS', 'AP', 'SP', 'ASAP', 'SPAS', 'SPAP'};
force_data = {AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};
force_test_participation = {
    [1,3,4,5];      % AS
    [2,4,5,6,7];    % AP
    [2,3,6,7];      % SP
    [1,4,5];        % ASAP
    [3,6,7];        % SPAS
    [2,6,7]         % SPAP
};

%{
Old force_colors definition - now using centralized colors from top of script
force_colors = [ ...
    hex2rgb('#FCCC73');     % AS - light golden yellow  
    hex2rgb('#FC9933');     % AP - warm orange
    hex2rgb('#E65940');     % SP - red-orange
    hex2rgb('#CC1A1A');     % ASAP - burgundy/wine
    hex2rgb('#801A33');     % SPAS - wine
    hex2rgb('#333359');     % SPAP - dark teal/navy
];
%}
%{
Old commented force_colors = [
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];
%}
number_of_pins = 8;
marker_types = {'o', 's', 'd', '^', 'v', '>', '<'};
marker_colors = lines(7);

figure('Position', [100, 100, 1400, 700]);

for interv_idx = 1:length(interventions)
    subplot(2, 3, interv_idx);
    
    curr_data = force_data{interv_idx};
    curr_color = force_colors(interv_idx, :);
    participating_tests = force_test_participation{interv_idx};
    
    % Calculate force differences for each pin (per-test means -> std across tests)
    number_of_pins = 8;
    avg_force_diff = zeros(1, number_of_pins);
    std_force_diff = zeros(1, number_of_pins);

    % Determine which participating tests actually have data for this intervention
    tests_present = intersect(unique([curr_data.Test]), participating_tests);

    % Build per-test x pin matrix of differences (NaN where missing)
    perTestDiffs = NaN(numel(tests_present), number_of_pins);
    for k = 1:numel(tests_present)
        t = tests_present(k);
        interv_mask = [curr_data.Test] == t;
        diseased_mask = [diseased_data.Test] == t;
        if ~any(interv_mask) || ~any(diseased_mask)
            continue;
        end
        for pin = 1:number_of_pins
            force_col = ['Force' num2str(pin)];
            interv_mean = mean([curr_data(interv_mask).(force_col)]);
            diseased_mean = mean([diseased_data(diseased_mask).(force_col)]);
            perTestDiffs(k, pin) = interv_mean - diseased_mean;
        end
    end

    % Compute mean and std across tests (ignore NaNs)
    for pin = 1:number_of_pins
        vals = perTestDiffs(:, pin);
        vals = vals(~isnan(vals));
        if isempty(vals)
            avg_force_diff(pin) = NaN;
            std_force_diff(pin) = NaN;
        else
            avg_force_diff(pin) = mean(vals);          % mean of per-test means
            if numel(vals) > 1
                std_force_diff(pin) = std(vals);      % std across per-test means
            else
                std_force_diff(pin) = NaN;            % not enough samples to define std
            end
        end
    end

    % Plot bars
    x = 1:number_of_pins;
    b = bar(x, avg_force_diff, 'FaceColor', curr_color, 'FaceAlpha', 1);
    hold on;

    % Error bars (skip NaNs)
    err = std_force_diff;
    err(isnan(err)) = 0;
    errorbar(x, avg_force_diff, err, 'k', 'LineWidth', 1.2, 'LineStyle', 'none');

    %{
    % Plot individual test markers using the perTestDiffs matrix
    for k = 1:numel(tests_present)
        test_num = tests_present(k);
        test_idx = find(participating_tests == test_num, 1);
        if isempty(test_idx), continue; end
        x_offset = (test_idx - (length(participating_tests)+1)/2) * 0.08;
        for pin = 1:number_of_pins
            val = perTestDiffs(k, pin);
            if isnan(val), continue; end
            scatter(pin + x_offset, val, 60, marker_types{test_num}, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', marker_colors(test_num,:), 'LineWidth', 1.1);
        end
    end
    %}
    
    % Formatting
    title(interventions{interv_idx}, 'FontWeight', 'bold', 'FontSize', 12);
    
    % Only add x and y labels to subplot 4 (ASAP)
    if interv_idx == 4
        xlabel('Pin Number');
        ylabel('\Delta F (N)');
    end
    
    set(gca, 'XTick', 1:number_of_pins);
    grid on;
    set(gca, 'FontSize', 10);
    hold off;
end

pubPlot('Width','double','Height',400,'Filename','ForceDiff_Pin','FileExtension',{'.png','.pdf'});

%% Force Contour Spline Visualization (6 Subplots)

% Manual position settings for the 8 force measurement pins + 1 shape control point 
% Arranged in circular pattern with Pin 1 at 12 o'clock, going counterclockwise
force_measurement_pin_positions = [
     0.0,  1.0;     % Pin 1 - 12 o'clock
     0.5,  0.75;    % Pin 2 - 1:30 position
     0.8,  0.3;     % Pin 3 - 3 o'clock
     0.85, -0.3;    % Pin 4 - 4:30 position
     0.58, -0.8;    % Pin 5 - 6 o'clock
     0.0,  -1.0;    % Pin 6 - 7:30 position
    -0.7,  -0.5;    % Pin 7 - 9 o'clock
    -0.7,   0.5;    % Pin 8 - 10:30 position
    -0.55,  0.85    % Point 9 - 11 o'clock (shape control point, no force measurement)
];

% Calculate polar angles from Cartesian coordinates for spline interpolation
polar_angles_all_points = atan2(force_measurement_pin_positions(:,2), force_measurement_pin_positions(:,1))';
reference_radius = 1.0; % Reference radius for scaling calculations

% Create figure with 6 subplots for different intervention types
figure('Position', [100, 100, 1400, 700]);

for intervention_index = 1:length(interventions)
    subplot(2, 3, intervention_index);
    
    current_intervention_data = force_data{intervention_index};
    tests_participating_in_intervention = force_test_participation{intervention_index};
    
    % Calculate average force differences for this intervention (only for 8 force measurement pins)
    average_force_differences = zeros(1, 8); % Only 8 force measurement pins
    for pin_number = 1:8
        force_column_name = ['Force' num2str(pin_number)];
        force_differences_across_tests = [];
        
        for test_number = tests_participating_in_intervention
            % Get intervention and diseased state forces for this test
            intervention_data_mask = [current_intervention_data.Test] == test_number;
            diseased_data_mask = [diseased_data.Test] == test_number;
            
            if any(intervention_data_mask) && any(diseased_data_mask)
                intervention_mean_force = mean([current_intervention_data(intervention_data_mask).(force_column_name)]);
                diseased_mean_force = mean([diseased_data(diseased_data_mask).(force_column_name)]);
                force_differences_across_tests(end+1) = intervention_mean_force - diseased_mean_force;
            end
        end
        
        if ~isempty(force_differences_across_tests)
            average_force_differences(pin_number) = mean(force_differences_across_tests);
        end
    end
    
    % Create baseline spline curve that passes through all pin positions (including shape control point)
    % Use parametric interpolation (x,y) over angle to avoid radial oscillations
    % Sort pins by angle to ensure monotonic parameterization for smooth interpolation
    normalized_angles = mod(polar_angles_all_points, 2*pi);
    [sorted_angles, angle_sort_indices] = sort(normalized_angles);
    sorted_pin_positions = force_measurement_pin_positions(angle_sort_indices, :);

    % Build extended angle vector for periodic interpolation (wraps around circle)
    extended_angles_for_interpolation = [sorted_angles - 2*pi, sorted_angles, sorted_angles + 2*pi];

    % Repeat point coordinates to match extended angles for smooth periodic interpolation
    extended_x_coordinates = [sorted_pin_positions(:,1); sorted_pin_positions(:,1); sorted_pin_positions(:,1)];
    extended_y_coordinates = [sorted_pin_positions(:,2); sorted_pin_positions(:,2); sorted_pin_positions(:,2)];

    % Interpolate x and y coordinates as functions of angle using spline for smooth closed curve
    fine_angle_resolution = linspace(0, 2*pi, 600);
    baseline_contour_x = interp1(extended_angles_for_interpolation, extended_x_coordinates, fine_angle_resolution, 'spline');
    baseline_contour_y = interp1(extended_angles_for_interpolation, extended_y_coordinates, fine_angle_resolution, 'spline');

    % Ensure curve closure for proper circular contour
    baseline_contour_x(end) = baseline_contour_x(1);
    baseline_contour_y(end) = baseline_contour_y(1);

    % Calculate contour displacement based on force measurements (only for 8 force measurement pins)
    contour_displacement_scale = 0.3; % Scale factor for contour displacement magnitude
    maximum_absolute_force = max(abs(average_force_differences));
    if maximum_absolute_force > 0
        normalized_force_magnitudes = average_force_differences / maximum_absolute_force;
    else
        normalized_force_magnitudes = zeros(1, 8);
    end

    % For force-based displacement, only use the first 8 points (force measurement pins, not shape control)
    force_pin_positions_only = force_measurement_pin_positions(1:8, :);
    force_pin_polar_angles = atan2(force_pin_positions_only(:,2), force_pin_positions_only(:,1))';
    
    % Sort force measurement pins by polar angle for consistent interpolation
    normalized_force_pin_angles = mod(force_pin_polar_angles, 2*pi);
    [sorted_force_pin_angles, force_pin_sort_indices] = sort(normalized_force_pin_angles);
    sorted_force_pin_positions = force_pin_positions_only(force_pin_sort_indices, :);

    % Compute displaced pin coordinates based on force measurements
    force_pin_vectors_from_origin = sorted_force_pin_positions; % vectors from origin to each pin
    force_pin_distances_from_origin = sqrt(sum(force_pin_vectors_from_origin.^2, 2));
    % Avoid division by zero for pins at origin
    unit_vectors_to_pins = force_pin_vectors_from_origin ./ max(force_pin_distances_from_origin, eps);

    % Ensure normalized forces are in same sorted order as force pin positions
    normalized_forces_sorted_order = normalized_force_magnitudes(force_pin_sort_indices);

    % Create displaced pin positions based on force magnitudes (radial displacement)
    displaced_force_pin_positions = sorted_force_pin_positions + (unit_vectors_to_pins .* (normalized_forces_sorted_order(:) * contour_displacement_scale));

    % Build extended angle vector for periodic interpolation of force-displaced contour
    extended_force_pin_angles = [sorted_force_pin_angles - 2*pi, sorted_force_pin_angles, sorted_force_pin_angles + 2*pi];

    % Extend displaced pin coordinates for periodic interpolation
    extended_displaced_x_coordinates = [displaced_force_pin_positions(:,1); displaced_force_pin_positions(:,1); displaced_force_pin_positions(:,1)];
    extended_displaced_y_coordinates = [displaced_force_pin_positions(:,2); displaced_force_pin_positions(:,2); displaced_force_pin_positions(:,2)];

    % Interpolate force-displaced contour smoothly using force pin angles
    interpolation_angle_resolution = linspace(0, 2*pi, 600);
    force_displaced_contour_x = interp1(extended_force_pin_angles, extended_displaced_x_coordinates, interpolation_angle_resolution, 'makima');
    force_displaced_contour_y = interp1(extended_force_pin_angles, extended_displaced_y_coordinates, interpolation_angle_resolution, 'makima');

    % Close force-displaced contour for proper circular shape
    force_displaced_contour_x(end) = force_displaced_contour_x(1);
    force_displaced_contour_y(end) = force_displaced_contour_y(1);

    % Apply light smoothing to remove interpolation artifacts and ensure smooth curve
    force_displaced_contour_x = smoothdata(force_displaced_contour_x, 'loess', 30);
    force_displaced_contour_y = smoothdata(force_displaced_contour_y, 'loess', 30);

    % Re-close contour after smoothing to maintain circular shape
    force_displaced_contour_x(end) = force_displaced_contour_x(1);
    force_displaced_contour_y(end) = force_displaced_contour_y(1);
    
    % Plot baseline spline curve (diseased state reference)
    plot(baseline_contour_x, baseline_contour_y, 'k-', 'LineWidth', 2);
    hold on;
    
    % Plot force-displaced contour curve (intervention state)
    %plot(force_displaced_contour_x, force_displaced_contour_y, '-.', 'Color', force_colors(intervention_index,:), 'LineWidth', 2);
    
    % Add valve leaflet splines
    % First spline: from between pins 8 and 1 to pin 4
    % Find midpoint between pins 8 and 1 on the baseline contour
    pin1_pos = force_measurement_pin_positions(1, :);  % Pin 1
    pin8_pos = force_measurement_pin_positions(8, :);  % Pin 8
    pin4_pos = force_measurement_pin_positions(4, :);  % Pin 4
    pin6_pos = force_measurement_pin_positions(6, :);  % Pin 6
    pin7_pos = force_measurement_pin_positions(7, :);  % Pin 7
    
    % Calculate midpoint between pins 8 and 1 (anterior commissure)
    midpoint_8_1 = (pin8_pos + pin1_pos) / 2;
    % Normalize to be on the contour circle
    midpoint_8_1 = midpoint_8_1 / norm(midpoint_8_1);
    
    % Calculate midpoint between pins 7 and 6 (posterior commissure)
    midpoint_7_6 = (pin7_pos + pin6_pos) / 2;
    % Normalize to be on the contour circle
    midpoint_7_6 = midpoint_7_6 / norm(midpoint_7_6);
    
    % First leaflet spline: anterior commissure to pin 4 (anterior leaflet)
    spline1_x = [midpoint_8_1(1), 0, pin4_pos(1)-0.05];
    spline1_y = [midpoint_8_1(2), 0, pin4_pos(2)];
    spline1_t = linspace(0, 1, 50);
    spline1_interp_x = interp1([0, 0.5, 1], spline1_x, spline1_t, 'spline');
    spline1_interp_y = interp1([0, 0.5, 1], spline1_y, spline1_t, 'spline');
    
    % Second leaflet spline: posterior commissure to center point of first spline (posterior leaflet)
    center_first_spline = [spline1_interp_x(30), spline1_interp_y(30)]; % midpoint of first spline
    
    % Add a control point to create a curve in the second spline
    % Position the control point to create a realistic curved leaflet
    control_point = [-0.2, -0.4]; % Offset downward for curve
    
    spline2_x = [-0.4, control_point(1), center_first_spline(1)];
    spline2_y = [-0.8, control_point(2), center_first_spline(2)];
    spline2_t = linspace(0, 1, 25);
    spline2_interp_x = interp1([0, 0.5, 1], spline2_x, spline2_t, 'spline');
    spline2_interp_y = interp1([0, 0.5, 1], spline2_y, spline2_t, 'spline');
    
    % Plot the valve leaflet splines
    plot(spline1_interp_x, spline1_interp_y, 'k-', 'LineWidth', 1.5);
    plot(spline2_interp_x, spline2_interp_y, 'k-', 'LineWidth', 1.5);
    
    % ========== RECTANGLE CONFIGURATION - EDIT ALL VALUES HERE ==========
    % Shared Rectangle Size
    rect_width = 0.4;               % Rectangle width (used for all rectangles)
    rect_height = 0.15;             % Rectangle height (used for all rectangles)
    
    % AS Rectangle Parameters (appears on AS, ASAP, SPAS subplots)
    AS_rect_center_x = -0.11;        % Rectangle center X coordinate
    AS_rect_center_y = 0.13;         % Rectangle center Y coordinate
    AS_rect_angle = 37;             % Rotation angle in degrees (positive = counterclockwise)
    % AS_rect_color defined in centralized color section at top of script

    % SP Rectangle Parameters (appears on SP, SPAP, SPAS subplots)
    SP_rect_center_x = -0.2;        % Rectangle center X coordinate
    SP_rect_center_y = -0.4;        % Rectangle center Y coordinate
    SP_rect_angle = -40;            % Rotation angle in degrees (negative = clockwise)
    % SP_rect_color defined in centralized color section at top of script

    % AP Rectangle Parameters (appears on AP, ASAP, SPAP subplots)
    AP_rect_center_x = 0.4;         % Rectangle center X coordinate
    AP_rect_center_y = -0.23;        % Rectangle center Y coordinate
    AP_rect_angle = 72;             % Rotation angle in degrees (positive = counterclockwise)
    % AP_rect_color defined in centralized color section at top of script
    
    % Shared Rectangle Appearance
    rect_edge_color = 'k';          % Edge color (black)
    rect_line_width = 1.5;          % Edge line width
    rect_face_alpha = 1.0;          % Transparency (0 = transparent, 1 = opaque)
    % ====================================================================
    
    % Add rectangle for AS-related configurations (AS, ASAP, SPAS)
    if any(strcmp(interventions{intervention_index}, {'AS', 'ASAP', 'SPAS'}))
        % Create rectangle vertices centered at origin
        half_width = rect_width / 2;
        half_height = rect_height / 2;
        rect_x_orig = [-half_width, half_width, half_width, -half_width, -half_width];
        rect_y_orig = [-half_height, -half_height, half_height, half_height, -half_height];
        
        % Apply rotation transformation
        angle_rad = AS_rect_angle * pi / 180;  % Convert to radians
        cos_angle = cos(angle_rad);
        sin_angle = sin(angle_rad);
        
        % Rotate each vertex
        rect_x_rotated = rect_x_orig * cos_angle - rect_y_orig * sin_angle;
        rect_y_rotated = rect_x_orig * sin_angle + rect_y_orig * cos_angle;
        
        % Translate to final position
        rect_x = rect_x_rotated + AS_rect_center_x;
        rect_y = rect_y_rotated + AS_rect_center_y;
        
        % Plot rectangle
        fill(rect_x, rect_y, AS_rect_color, 'EdgeColor', rect_edge_color, 'LineWidth', rect_line_width, 'FaceAlpha', rect_face_alpha);
    end
    
    % Add rectangle for SP-related configurations (SP, SPAP, SPAS)
    if any(strcmp(interventions{intervention_index}, {'SP', 'SPAP', 'SPAS'}))
        % Create rectangle vertices centered at origin
        half_width_sp = rect_width / 2;
        half_height_sp = rect_height / 2;
        rect_x_orig_sp = [-half_width_sp, half_width_sp, half_width_sp, -half_width_sp, -half_width_sp];
        rect_y_orig_sp = [-half_height_sp, -half_height_sp, half_height_sp, half_height_sp, -half_height_sp];
        
        % Apply rotation transformation
        angle_rad_sp = SP_rect_angle * pi / 180;  % Convert to radians
        cos_angle_sp = cos(angle_rad_sp);
        sin_angle_sp = sin(angle_rad_sp);
        
        % Rotate each vertex
        rect_x_rotated_sp = rect_x_orig_sp * cos_angle_sp - rect_y_orig_sp * sin_angle_sp;
        rect_y_rotated_sp = rect_x_orig_sp * sin_angle_sp + rect_y_orig_sp * cos_angle_sp;
        
        % Translate to final position
        rect_x_sp = rect_x_rotated_sp + SP_rect_center_x;
        rect_y_sp = rect_y_rotated_sp + SP_rect_center_y;
        
        % Plot rectangle
        fill(rect_x_sp, rect_y_sp, SP_rect_color, 'EdgeColor', rect_edge_color, 'LineWidth', rect_line_width, 'FaceAlpha', rect_face_alpha);
    end
    
    % Add rectangle for AP-related configurations (AP, ASAP, SPAP)
    if any(strcmp(interventions{intervention_index}, {'AP', 'ASAP', 'SPAP'}))
        % Create rectangle vertices centered at origin
        half_width_ap = rect_width / 2;
        half_height_ap = rect_height / 2;
        rect_x_orig_ap = [-half_width_ap, half_width_ap, half_width_ap, -half_width_ap, -half_width_ap];
        rect_y_orig_ap = [-half_height_ap, -half_height_ap, half_height_ap, half_height_ap, -half_height_ap];
        
        % Apply rotation transformation
        angle_rad_ap = AP_rect_angle * pi / 180;  % Convert to radians
        cos_angle_ap = cos(angle_rad_ap);
        sin_angle_ap = sin(angle_rad_ap);
        
        % Rotate each vertex
        rect_x_rotated_ap = rect_x_orig_ap * cos_angle_ap - rect_y_orig_ap * sin_angle_ap;
        rect_y_rotated_ap = rect_x_orig_ap * sin_angle_ap + rect_y_orig_ap * cos_angle_ap;
        
        % Translate to final position
        rect_x_ap = rect_x_rotated_ap + AP_rect_center_x;
        rect_y_ap = rect_y_rotated_ap + AP_rect_center_y;
        
        % Plot rectangle
        fill(rect_x_ap, rect_y_ap, AP_rect_color, 'EdgeColor', rect_edge_color, 'LineWidth', rect_line_width, 'FaceAlpha', rect_face_alpha);
    end
    
    % Add force magnitude bars using thick lines (only for force measurement pins 1-8)
    bar_length_scale_factor = 0.4; % Scale factor for bar length visualization
    maximum_absolute_force_for_bars = max(abs(average_force_differences));
    
    if maximum_absolute_force_for_bars > 0
        for pin_number = 1:8
            pin_x_position = force_measurement_pin_positions(pin_number, 1);
            pin_y_position = force_measurement_pin_positions(pin_number, 2);
            force_magnitude_at_pin = average_force_differences(pin_number);
            
            % Find the closest point on the spline to get normal direction
            pin_distances = sqrt((baseline_contour_x - pin_x_position).^2 + (baseline_contour_y - pin_y_position).^2);
            [~, closest_idx] = min(pin_distances);
            
            % Calculate tangent vector at the closest point using numerical derivative
            if closest_idx == 1
                % Use forward difference at start
                tangent_x = baseline_contour_x(2) - baseline_contour_x(1);
                tangent_y = baseline_contour_y(2) - baseline_contour_y(1);
            elseif closest_idx == length(baseline_contour_x)
                % Use backward difference at end
                tangent_x = baseline_contour_x(end) - baseline_contour_x(end-1);
                tangent_y = baseline_contour_y(end) - baseline_contour_y(end-1);
            else
                % Use central difference in middle
                tangent_x = baseline_contour_x(closest_idx+1) - baseline_contour_x(closest_idx-1);
                tangent_y = baseline_contour_y(closest_idx+1) - baseline_contour_y(closest_idx-1);
            end
            
            % Normalize tangent vector
            tangent_magnitude = sqrt(tangent_x^2 + tangent_y^2);
            if tangent_magnitude > eps
                tangent_unit_x = tangent_x / tangent_magnitude;
                tangent_unit_y = tangent_y / tangent_magnitude;
            else
                tangent_unit_x = 1;
                tangent_unit_y = 0;
            end
            
            % Calculate normal vector (perpendicular to tangent, pointing outward)
            % Rotate tangent 90 degrees clockwise to get outward normal
            normal_unit_x = tangent_unit_y;
            normal_unit_y = -tangent_unit_x;
            
            % Scale bar length by normalized force magnitude
            normalized_bar_length = (force_magnitude_at_pin / maximum_absolute_force_for_bars) * bar_length_scale_factor;
            bar_x_component = normal_unit_x * normalized_bar_length;
            bar_y_component = normal_unit_y * normalized_bar_length;
            
            % Only plot bar if its length exceeds minimum threshold for visibility
            minimum_bar_length_threshold = 0.1;
            if abs(normalized_bar_length) > minimum_bar_length_threshold
                % Calculate end point of the bar
                end_x = pin_x_position + bar_x_component;
                end_y = pin_y_position + bar_y_component;
                
                % Plot thick line from pin to end point - positive forces point outward normal to spline, negative forces point inward
                plot([pin_x_position, end_x], [pin_y_position, end_y], 'r-', 'LineWidth', 8);
            end
        end
    end
    
    % Plot pin positions with identifying numbers and markers
    for pin_number = 1:9
        pin_x_position = force_measurement_pin_positions(pin_number, 1);
        pin_y_position = force_measurement_pin_positions(pin_number, 2);
        
        if pin_number <= 8
            % Force measurement pins - black circular markers with pin numbers
            plot(pin_x_position, pin_y_position, '.', 'MarkerSize', 15, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            
            % Add pin number label positioned outside the pin location
            text(pin_x_position * 1.3, pin_y_position * 1.3, num2str(pin_number), 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
        end
    end
    
    % Format subplot appearance and labels
    title(interventions{intervention_index}, 'FontWeight', 'bold', 'FontSize', 12);
    axis equal;
    axis off;
    axis([-1.5 1.5 -1.5 1.5]);
    grid off;
    set(gca, 'FontSize', 10);
    set(gca, 'XTickLabels',{})
    set(gca, 'YTickLabels',{})
    hold off;
end

% Add overall title and legend

pubPlot('Width','double','Height',400,'Filename','ForceDiff_Contour','FileExtension',{'.png','.pdf'});


%{
%% Export Force Differences Relative to Diseased Data

% Initialize empty structure array for collecting force difference rows
diff_rows = struct('Intervention', {}, 'Test', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {}, 'ForceDifference', {}, 'PressureDifference', {}, 'FlowDifference', {});

% Get diseased data reference for each test
diseased_reference = containers.Map('KeyType', 'int32', 'ValueType', 'any');
for i = 1:length(diseased_data)
    test_num = diseased_data(i).Test;
    if ~isKey(diseased_reference, test_num)
        diseased_reference(test_num) = [];
    end
    diseased_reference(test_num) = [diseased_reference(test_num), diseased_data(i)];
end

% Process non-diseased interventions (skip 'Diseased' at index 1)
intervention_names = {'Diseased','AS', 'AP', 'SP', 'ASAP', 'SPAS', 'SPAP'};
intervention_datasets = {diseased_data, AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

for i = 1:numel(intervention_names)
    data = intervention_datasets{i};
    if isempty(data)
        continue;
    end
    
    tests = unique([data.Test]);
    for tt = tests
        % Check if we have diseased reference data for this test
        if ~isKey(diseased_reference, tt)
            warning('No diseased reference data found for test %d, skipping intervention %s', tt, intervention_names{i});
            continue;
        end
        
        % Get intervention data for this test
        interv_mask = [data.Test] == tt;
        interv_data = data(interv_mask);
        
        % Get diseased reference data for this test
        diseased_ref = diseased_reference(tt);
        
        if isempty(interv_data) || isempty(diseased_ref)
            continue;
        end
        
        % Calculate means for intervention data
        interv_pressure_mean = mean([interv_data.Pressure]);
        interv_flow_mean = mean([interv_data.FlowRate]);
        
        % Calculate means for diseased reference data
        diseased_pressure_mean = mean([diseased_ref.Pressure]);
        diseased_flow_mean = mean([diseased_ref.FlowRate]);
        
        % Calculate pressure and flow differences
        pressure_diff = interv_pressure_mean - diseased_pressure_mean;
        flow_diff = interv_flow_mean - diseased_flow_mean;
        
        % Calculate force differences for each pin
        for pin = 1:8
            force_col = ['Force' num2str(pin)];
            
            % Calculate mean forces
            interv_force_mean = mean([interv_data.(force_col)]);
            diseased_force_mean = mean([diseased_ref.(force_col)]);
            
            % Calculate force difference (intervention - diseased)
            force_diff = interv_force_mean - diseased_force_mean;
            
            diff_rows(end+1) = struct( ...
                'Intervention', intervention_names{i}, ...
                'Test', tt, ...
                'Pressure', interv_pressure_mean, ...
                'FlowRate', interv_flow_mean, ...
                'Pin', pin, ...
                'Force', interv_force_mean, ...
                'ForceDifference', force_diff, ...
                'PressureDifference', pressure_diff, ...
                'FlowDifference', flow_diff ...
            );
        end
    end
end

% Convert to table and write to Excel
if ~isempty(diff_rows)
    T_diff = struct2table(diff_rows);
    outFile_diff = fullfile(projectRoot, 'TriClipXT_Statistics.xlsx');
    try
        writetable(T_diff, outFile_diff, 'Sheet', 'ForceDifferences', 'WriteMode', 'overwrite');
        fprintf('Exported %d force difference rows to %s\n', height(T_diff), outFile_diff);
    catch ME
        warning('MATLAB:ExcelWrite', 'Failed to write force difference Excel file: %s. Attempting to save as CSV instead.', ME.message);
        csvFile_diff = fullfile(projectRoot, 'TriClipXT_ForceDifferences_vs_Diseased.csv');
        writetable(T_diff, csvFile_diff);
        fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
    end
else
    warning('No force difference rows collected for export; no file written.');
end
%}


