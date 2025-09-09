clc
clear all
close all

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

% Define interventions and their corresponding data
interventions = {'Diseased', 'AS', 'AP', 'SP', 'ASAP', 'SPAS', 'SPAP'};
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
bar_colors = [ ...
    0.2 0.2 0.2;        % Diseased - dark gray
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];

x = 1:num_interv;
b = bar(x, avg_flow, 'FaceColor', 'flat');
for i = 1:num_interv
    b.CData(i,:) = bar_colors(i,:);
end
hold on;

% Error bars
errorbar(x, avg_flow, std_flow, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

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
set(gcf, 'Color', 'white');
grid on;
hold off;

pubPlot('Filename','Flow_vs_Intervention','FileExtension',{'.png','.pdf'});


%% Pressure Bar Graph by Intervention Type

% Define interventions and their corresponding data
interventions = {'Diseased', 'AS', 'AP', 'SP', 'ASAP', 'SPAS', 'SPAP'};
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
bar_colors = [ ...
    0.2 0.2 0.2;        % Diseased - dark gray
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];

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

force_colors = [
    1.0 0.4 0.0;        % AS - #FF6700 orange
    0.95 0.16 0.31;     % AP - #F14F50 red
    0.71 0.04 0.45;     % SP - #B50A72 magenta
    0.39 0.22 0.58;     % ASAP - #643894 purple
    0.32 0.38 0.93;     % SPAS - #5260EE blue
    0.29 0.45 0.72;     % SPAP - #4A73B8 teal blue
];

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
    b = bar(x, avg_force_diff, 'FaceColor', curr_color, 'FaceAlpha', 0.8);
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
    xlabel('Pin Number');
    ylabel('Force Difference (N)');
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
    plot(force_displaced_contour_x, force_displaced_contour_y, '-.', 'Color', force_colors(intervention_index,:), 'LineWidth', 2);
    
    % Add force magnitude arrows using quiver (only for force measurement pins 1-8)
    arrow_length_scale_factor = 0.25; % Scale factor for arrow length visualization
    maximum_absolute_force_for_arrows = max(abs(average_force_differences));
    
    if maximum_absolute_force_for_arrows > 0
        for pin_number = 1:8
            pin_x_position = force_measurement_pin_positions(pin_number, 1);
            pin_y_position = force_measurement_pin_positions(pin_number, 2);
            force_magnitude_at_pin = average_force_differences(pin_number);
            
            % Calculate radial direction vector (outward from center for positive forces)
            pin_position_vector = [pin_x_position, pin_y_position];
            distance_from_origin = norm(pin_position_vector);
            if distance_from_origin > 0
                radial_unit_vector = pin_position_vector / distance_from_origin;
            else
                radial_unit_vector = [1, 0]; % fallback direction for pins at origin
            end
            
            % Scale arrow length by normalized force magnitude
            normalized_arrow_length = (force_magnitude_at_pin / maximum_absolute_force_for_arrows) * arrow_length_scale_factor;
            arrow_x_component = radial_unit_vector(1) * normalized_arrow_length;
            arrow_y_component = radial_unit_vector(2) * normalized_arrow_length;
            
            % Only plot arrow if its length exceeds minimum threshold for visibility
            minimum_arrow_length_threshold = 0.1;
            if abs(normalized_arrow_length) > minimum_arrow_length_threshold
                % Plot arrow - positive forces point outward, negative forces point inward
                quiver(pin_x_position, pin_y_position, arrow_x_component, arrow_y_component, 0, 'Color', 'k', ...
                    'LineWidth', 2, 'MaxHeadSize', 2);
            end
        end
    end
    
    % Plot pin positions with identifying numbers and markers
    for pin_number = 1:9
        pin_x_position = force_measurement_pin_positions(pin_number, 1);
        pin_y_position = force_measurement_pin_positions(pin_number, 2);
        
        if pin_number <= 8
            % Force measurement pins - red circular markers with pin numbers
            plot(pin_x_position, pin_y_position, '.', 'MarkerSize', 15, 'MarkerFaceColor', 'r', ...
                'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
            
            % Add pin number label positioned outside the pin location
            text(pin_x_position * 1.3, pin_y_position * 1.3, num2str(pin_number), 'FontSize', 10, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
        end
    end
    
    % Format subplot appearance and labels
    title(interventions{intervention_index}, 'FontWeight', 'bold', 'FontSize', 12);
    axis equal;
    axis([-1.5 1.5 -1.5 1.5]);
    grid on;
    set(gca, 'FontSize', 10);
    set(gca, 'XTickLabels',{})
    set(gca, 'YTickLabels',{})
    hold off;
end

% Add overall title and legend

pubPlot('Width','double','Height',400,'Filename','ForceDiff_Contour','FileExtension',{'.png','.pdf'});

%% Export averaged intervention data to Excel
%{
% Build intervention list and corresponding data cell array (reuse existing variables)
interventions = {'Diseased','AS','AP','SP','ASAP','SPAS','SPAP'};
intervention_data = {diseased_data, AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

% Initialize empty structure array for collecting rows
rows = struct('Intervention', {}, 'Test', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {});

for i = 1:numel(interventions)
    data = intervention_data{i};
    if isempty(data)
        continue;
    end
    tests = unique([data.Test]);
    for tt = tests
        mask = [data.Test] == tt;
        if ~any(mask)
            continue;
        end
        % Calculate means for this test
        pressure_mean = mean([data(mask).Pressure]);
        flow_mean = mean([data(mask).FlowRate]);
        
        % Add one row per pin for this test
        for pin = 1:8
            force_col = ['Force' num2str(pin)];
            force_mean = mean([data(mask).(force_col)]);
            rows(end+1) = struct( ...
                'Intervention', interventions{i}, ...
                'Test', tt, ...
                'Pressure', pressure_mean, ...
                'FlowRate', flow_mean, ...
                'Pin', pin, ...
                'Force', force_mean ...
            );
        end
    end
end

% Convert to table and write to Excel
if ~isempty(rows)
    T = struct2table(rows);
    outFile = fullfile('C:\Users\trace\OneDrive\Documents\Rausch Lab\TriClip Experiment', 'TriClipXT_Averages_by_Intervention.xlsx');
    try
        writetable(T, outFile, 'Sheet', 'Averages', 'WriteMode', 'overwrite');
        fprintf('Exported %d rows to %s\n', height(T), outFile);
    catch ME
        warning('MATLAB:ExcelWrite', 'Failed to write Excel file: %s. Attempting to save as CSV instead.', ME.message);
        csvFile = fullfile('C:\Users\trace\OneDrive\Documents\Rausch Lab\TriClip Experiment', 'TriClipXT_Averages_by_Intervention.csv');
        writetable(T, csvFile);
        fprintf('Exported %d rows to %s\n', height(T), csvFile);
    end
else
    warning('No rows collected for export; no file written.');
end


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


