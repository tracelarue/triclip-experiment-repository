clc
clear all
close all


%% ==========================================================================

% Get the full path of this script
scriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(scriptPath);
projectRoot = fileparts(scriptDir);

% Build data paths
dataDir = fullfile(projectRoot, 'data');
testIdxDir = fullfile(dataDir, 'testidx');
annotationsDir = fullfile(dataDir, 'videoAnnotations');

% -------------------------- Load stable region index data for each test
load(fullfile(testIdxDir, 'test1idx.mat'));
load(fullfile(testIdxDir, 'test2idx.mat'));
load(fullfile(testIdxDir, 'test3idx.mat'));
load(fullfile(testIdxDir, 'test4idx.mat'));
load(fullfile(testIdxDir, 'test5idx.mat'));
load(fullfile(testIdxDir, 'test6idx.mat'));
load(fullfile(testIdxDir, 'test7idx.mat'));
load(fullfile(testIdxDir, 'test8idx.mat'));
load(fullfile(testIdxDir, 'test9idx.mat'));
load(fullfile(testIdxDir, 'test10idx.mat'));
load(fullfile(testIdxDir, 'test11idx.mat'));
load(fullfile(testIdxDir, 'test12idx.mat'));
load(fullfile(testIdxDir, 'test13idx.mat'));

% -------------------------- Load video annotation data for each test
load(fullfile(annotationsDir, 'test1Annotations.mat'));
load(fullfile(annotationsDir, 'test2Annotations.mat'));
load(fullfile(annotationsDir, 'test3Annotations.mat'));
load(fullfile(annotationsDir, 'test4Annotations.mat'));
load(fullfile(annotationsDir, 'test5Annotations.mat'));
load(fullfile(annotationsDir, 'test6Annotations.mat'));
load(fullfile(annotationsDir, 'test7Annotations.mat'));
load(fullfile(annotationsDir, 'test8Annotations.mat'));
load(fullfile(annotationsDir, 'test9Annotations.mat'));
load(fullfile(annotationsDir, 'test10Annotations.mat'));
load(fullfile(annotationsDir, 'test11Annotations.mat'));
load(fullfile(annotationsDir, 'test12Annotations.mat'));
load(fullfile(annotationsDir, 'test13Annotations.mat'));
annotationValveNames = ['2025_05_21';'2025_06_13';'2025_07_09';...
    '2025_07_16';'2025_07_23';'2025_08_05';'2025_08_07';...
    '2026_01_14';'2026_01_21';'2026_01_23';'2026_01_30';...
    '2026_02_04';'2026_02_06'];

% -------------------------- Load morphology data for each test
load(fullfile(dataDir,'TriClipMorphology.mat'));

% -------------------------- Load and preprocess data

% Define paths to data files using relative paths
folderPaths = {
    fullfile(dataDir, '05_21_25'),
    fullfile(dataDir, '06_13_25'),
    fullfile(dataDir, '07_09_25'),
    fullfile(dataDir, '07_16_25'),
    fullfile(dataDir, '07_23_25'),
    fullfile(dataDir, '08_05_25'),
    fullfile(dataDir, '08_07_25'),
    fullfile(dataDir, '01_14_26'),
    fullfile(dataDir, '01_21_26'),
    fullfile(dataDir, '01_23_26'),
    fullfile(dataDir, '01_30_26'),
    fullfile(dataDir, '02_04_26'),
    fullfile(dataDir, '02_06_26')
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

force_calibration_updated = [
    -1581.340707*3.3  % pin1
    -1356.716945*3.3  % pin2
    -1436.661363*3.3  % pin3
    -943.844378*3.3  % pin4
    -939.596181*3.3  % pin5
    -803.385345*3.3  % pin6
    -2958.500987*3.3  % pin7
    -1410.673529*3.3  % pin8
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
            if test < 8
                all_test_data{test}{subtest}.(force_x) = all_test_data{test}{subtest}.(force_x) * force_calibration(pin);
            else
                all_test_data{test}{subtest}.(force_x) = all_test_data{test}{subtest}.(force_x) * force_calibration_updated(pin);
            end
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

% NOTES:
% Test 1 -  AP Data is bad, damaged septal leaflet on AS removal
%           exclude 16-18 (AP) for noise
% Test 2 -  
% Test 3 -  ripped small chunk of anterior leaflet, 
%           consistent noise/vibration 10-12, 16-18
% Test 4 -  lost some anterior leaflet, restarted and re-ran disease
% Test 5 -  exclude 15 for creep in #7, little bit of noise 10-12
% Test 6 -  10-12 disease doesn't plateau nicely
%           exclude 1-4 due to noise
% Test 7 -  unreliable healhty flow data, lost some septal leafelt on AS removal
%           22-24 doesn't plateau that nice, small noise 10-12, 19-21
% Test 8 -  AS didn't work 
%           exclude 19-21, 17, 04-06 for noise. 10-12 and 16-18 also look a little sus
% Test 9 - 
% Test 10 - exclude 15, 17. 16-21 look like noise, very small values
% Test 11 - 
% Test 12 - 19 has some noise on pin #7
% Test 13 - 

% Data Indices based on the order of the files in the folders
%         13:15 16:18    19:21  22:24       
% Test 1 -  AS,  AP,     ASAP
% Test 2 -  SP,  AP,     SPAP
% Test 3 -  AS,  SP,     SPAS
% Test 4 -  AS,  AP,     ASAP
% Test 5 -  AS,  AP,     ASAP
% Test 6 -  SP,  SPAP,   SPAS,  AP
% Test 7 -  SP,  SPAP,   SPAS,  AP 
% Test 8 -  SP,  SPAS,   AS 
% Test 9 -  SP,  SPAP,   AP,    ASAP,   AS
% Test 10 - AS,  ASAP,   AP,    SPAP,   SP
% Test 11 - SP,  SPAS,   AS
% Test 12 - AS,  ASAP,   AP,    SPAP,   SP
% Test 13 - AS,  SPAS,   SP,    SPAP,   AP

healthy_data = [test1avg(1:9), test2avg(1:9), test3avg(1:9), test4avg(1:9), test5avg(1:9), test6avg(5:9), test7avg(1:9), test8avg([1:3,7:9]), test9avg(1:9), test10avg([1:3,7:9]), test11avg(1:9), test12avg(1:9), test13avg(1:9)];
diseased_data = [test1avg(10:12), test2avg(10:12), test3avg(10:12), test4avg(10:12), test5avg(10:12), test6avg(10:12), test7avg(10:12), test8avg(10:12), test9avg(10:12), test10avg(10:12), test11avg(10:12), test12avg(10:12), test13avg(10:12)];
AS_data= [test1avg(13:15), test3avg(13:15), test4avg(13:15), test5avg(13:14), test9avg(25:27), test10avg(13:14), test11avg(19:21), test12avg(13:15), test13avg(13:15)];
AP_data= [test2avg(16:18), test4avg(16:18), test5avg(16:18), test6avg(22:24), test7avg(22:24), test9avg(19:21), test10avg(19:21), test12avg(19:21), test13avg(25:27)];
SP_data= [test2avg(13:15), test3avg(16:18), test6avg(13:15), test7avg(13:15), test8avg(13:15), test9avg(13:15), test10avg(25:27), test11avg(13:15), test12avg(25:27), test13avg(19:21)];
ASAP_data= [test1avg(19:21), test4avg(19:21), test5avg(19:21), test9avg(22:24), test10avg([16,18]), test12avg(16:18)];
SPAS_data= [test3avg(19:21), test6avg(19:21), test7avg(19:21), test8avg([16,18]), test11avg(16:18), test13avg(16:18)];
SPAP_data= [test2avg(19:21), test6avg(16:18), test7avg(16:18), test9avg(16:18), test10avg(22:24), test12avg(22:24), test13avg(22:24)];


test_order = {
    {'AS','ASAP','AP'};
    {'AP','SPAP','SP'};
    {'SP','SPAS','AS'};
    {'AS','ASAP','AP'};
    {'AS','ASAP','AP'};
    {'SP','SPAS','SPAP','AP'};
    {'SP','SPAS','SPAP','AP'};
    {'SP','SPAS','AS'};
    {'SP','SPAP','AP','ASAP','AS'};
    {'AS','ASAP','AP','SPAP','SP'};
    {'SP','SPAS','AS'};
    {'AS','ASAP','AP','SPAP','SP'};
    {'AS','SPAS','SP','SPAP','AP'};
    };

% treatment procedure
% 1: AS -> ASAP -> AP
% 2: SP -> SPAS -> AS
% 3: AP -> SPAP -> SP
% 4: SP -> SPAP -> AP
% 5: AP -> ASAP -> AS
% 6: SP -> SPAS -> SPAP -> AP
% 7: SP -> SPAP -> AP -> ASAP -> AS
% 8: AS -> ASAP -> AP -> SPAP -> SP
% 9: AS -> SPAS -> SP -> SPAP -> AP
test_procedure = {
    [1]; % Test 1
    [3]; % Test 2
    [2]; % Test 3
    [1]; % Test 4
    [1]; % Test 5
    [6]; % Test 6
    [6]; % Test 7
    [2]; % Test 8
    [7]; %[4, 5]; % Test 9
    [8];%[1, 3]; % Test 10
    [2]  % Test 11
    [8]  % Test 12
    [9]  % Test 13
    };



% Which tests participate in each intervention (cell array, 1 row per intervention)
test_participation = {
    [1:13];          % Diseased
    [1,3,4,5,9,10,11,12,13];   % AS
    [2,4,5,6,7,9,10,12,13];      % AP
    [2,3,6,7,8,9,10,11,12,13];   % SP
    [1,4,5,9,10,12];             % ASAP
    [3,6,7,8,11,13];             % SPAS
    [2,6,7,9,10,12,13];          % SPAP
};

% % Tests with flow difference < 5
% % AS: 8, 10
% % AP: 2, 9
% % ASAP: 4, 10 
% AS_data= [test1avg(13:15), test3avg(13:15), test4avg(13:15), test5avg(13:15), test9avg(25:27), test10avg(13:14), test11avg(19:21)];
% AP_data= [test4avg(16:18), test5avg(16:18), test6avg(22:24), test7avg(22:24), test10avg(19:21)];
% SP_data= [test2avg(13:15), test3avg(16:18), test6avg(13:15), test7avg(13:15), test8avg(13:15), test9avg(13:15), test10avg(25:27), test11avg(13:15)];
% ASAP_data= [test1avg(19:21), test5avg(19:21), test9avg(22:24)];
% SPAS_data= [test3avg(19:21), test6avg(19:21), test7avg(19:21), test8avg([16,18]), test11avg(16:18)];
% SPAP_data= [test2avg(19:21), test6avg(16:18), test7avg(16:18), test9avg(16:18), test10avg(22:24)];
% 
% % Which tests participate in each intervention (cell array, 1 row per intervention)
% test_participation = {
%     [1:10];          % Diseased
%     [1,3,4,5,9];      % AS
%     [4,5,6,7,10];      % AP
%     [2,3,6,7,8,9,10];      % SP
%     [1,5,9];          % ASAP
%     [3,6,7,8];             % SPAS
%     [2,6,7,9,10];          % SPAP
% };

% % Tests with flow difference < 10
% % AS: 3, 8, 10
% % AP: 2, 4, 9, 10
% % SP: 3, 6, 11
% % ASAP: 4, 9, 10 
% % SPAS: 11
% AS_data= [test1avg(13:15), test3avg(13:15), test4avg(13:15), test5avg(13:15), test9avg(25:27), test10avg(13:14), test11avg(19:21)];
% AP_data= [test5avg(16:18), test6avg(22:24), test7avg(22:24)];
% SP_data= [test2avg(13:15), test7avg(13:15), test8avg(13:15), test9avg(13:15), test10avg(25:27)];
% ASAP_data= [test1avg(19:21), test5avg(19:21)];
% SPAS_data= [test3avg(19:21), test6avg(19:21), test7avg(19:21), test8avg([16,18])];
% SPAP_data= [test2avg(19:21), test6avg(16:18), test7avg(16:18), test9avg(16:18), test10avg(22:24)];
% 
% test_participation = {
%     [1:10];          % Diseased
%     [1,4,5,9];      % AS
%     [5,6,7];      % AP
%     [2,7,8,9,10];      % SP
%     [1,5];          % ASAP
%     [3,6,7,8];             % SPAS
%     [2,6,7,9,10];          % SPAP
% };


force_test_participation = {test_participation{2:end}};

% Hex to RGB conversion function
hex2rgb = @(hex) sscanf(hex(2:end), '%2x%2x%2x', [1 3]) / 255;

% ========================== COLOR CONTROL SECTION ==========================
% Define all colors for consistent plotting across the entire analysis
% Modify these hex codes to change colors throughout all plots

% Intervention Colors (used in all bar graphs and contour plots)
COLORS = struct();
% COLORS.Diseased = hex2rgb('#EDE6CC');   % Light cream/beige
% COLORS.AS       = hex2rgb('#FCCC73');   % Light golden yellow
% COLORS.AP       = hex2rgb('#FC9933');   % Warm orange
% COLORS.SP       = hex2rgb('#E65940');   % Red-orange
% COLORS.ASAP     = hex2rgb('#CC1A1A');   % Burgundy/wine
% COLORS.SPAS     = hex2rgb('#801A33');   % Wine
% COLORS.SPAP     = hex2rgb('#333359');   % Dark teal/navy
% 
% % Rainbow Neon Fusion
% COLORS.Diseased = hex2rgb('#333359');   % Light cream/beige
% COLORS.AS       = hex2rgb('#8338EC');   % Light golden yellow
% COLORS.AP       = hex2rgb('#3A86FF');   % Warm orange
% COLORS.SP       = hex2rgb('#619EFF');   % Red-orange
% COLORS.ASAP     = hex2rgb('#FFBE0B');   % Burgundy/wine
% COLORS.SPAS     = hex2rgb('#FB5607');   % Wine
% COLORS.SPAP     = hex2rgb('#FF006E');   % Dark teal/navy
% 
% % Tropical Sunset Mix
% COLORS.Diseased = hex2rgb('#333359');   % Light cream/beige
% COLORS.AS       = hex2rgb('#219ebc');   % Light golden yellow
% COLORS.AP       = hex2rgb('#A8E8F9');   % Warm orange
% COLORS.SP       = hex2rgb('#39B89A');   % Red-orange
% COLORS.ASAP     = hex2rgb('#79D3BE');   % Burgundy/wine
% COLORS.SPAS     = hex2rgb('#FF5883');   % Wine
% COLORS.SPAP     = hex2rgb('#FF91AD');   % Dark teal/navy

% Rainbow Neon Fusion
COLORS.Diseased = hex2rgb('#333359');   % Dark teal/navy
COLORS.AS       = hex2rgb('#741FEA');
COLORS.AP       = hex2rgb('#4455EE');
COLORS.SP       = hex2rgb('#85B4FF');
COLORS.ASAP     = hex2rgb('#FFBE0B');
COLORS.SPAS     = hex2rgb('#FB5607');
COLORS.SPAP     = hex2rgb('#FF006E');   

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
interventions = {'Dis', 'AS', 'AP', 'SP ', 'ASAP', 'SPAS', 'SPAP'};
intervention_data = {diseased_data, AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

%% Pressure Bar Graph by Intervention Type

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
        std_pressure(i) = std(perTestMeans)./sqrt(numel(tests)); % standard error across tests
    else
        std_pressure(i) = 0;
    end
end

% Plot setup
figure('Position', [100, 100, 1200, 700]);
subplot(1,2,1)

x = 1:num_interv;
b = bar(x, avg_pressure, 'FaceColor', 'flat');
for i = 1:num_interv
    b.CData(i,:) = bar_colors(i,:);
end
hold on;

% Error bars
errorbar(x, avg_pressure, std_pressure, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

% Add significance lines for comparisons to diseased state
sig_interventions = [2]; % AS
sig_symbols = {'***'}; % Significance symbols

% Calculate height for significance lines (above error bars)
max_bar_height = max(avg_pressure + std_pressure);
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
    text(midpoint_x, line_height - 3, sig_symbols{i}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
end

% Axis labels and formatting
xlabel('Intervention', 'FontSize', 12);
ylabel('Pressure (mmHg)', 'FontSize', 12);
ylim([0,40])
set(gca, 'XTick', x, 'XTickLabel', interventions, 'XTickLabelRotation',30);
set(gcf, 'Color', 'white');
grid on;
hold off;

% pubPlot('Filename','Intervention','FileExtension',{'.png','.eps'});
% pubPlot('Width','double','Height',300,'Filename','Intervention','FileExtension',{'.png','.eps'});

%% Flow Rate Bar Graph by Intervention Type

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
        std_flow(i) = std(perTestMeans)./sqrt(numel(tests)); % standard error across tests
    else
        std_flow(i) = 0;
    end
end

% Plot setup
% figure('Position', [100, 100, 1200, 700]);
subplot(1,2,2)

x = 1:num_interv;
b = bar(x, avg_flow, 'FaceColor', 'flat');
for i = 1:num_interv
    b.CData(i,:) = bar_colors(i,:);
end
hold on;

% Error bars
errorbar(x, avg_flow, std_flow, 'k', 'LineWidth', 1.5, 'LineStyle', 'none');

% Add significance lines for comparisons to diseased state
sig_interventions = [];
sig_symbols = {}; % Significance symbols

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

% Axis labels and formatting
xlabel('Intervention', 'FontSize', 12);
ylabel('Flow Rate (ml/s)', 'FontSize', 12);
set(gca, 'XTick', x, 'XTickLabel', interventions, 'XTickLabelRotation',30);

% Adjust y-axis limits to accommodate significance lines
ylim([0,60])

set(gcf, 'Color', 'white');
grid on;
hold off;

% pubPlot('SpacingOffset',1,'Filename','Flow_vs_Intervention','FileExtension',{'.png','.eps'});
pubPlot('Width','double','Height',300,'Filename','Intervention','FileExtension',{'.png','.eps'});

%% Force Difference Bar Graphs by Intervention Type (6 Subplots)

% List of intervention types to plot (excluding 'Diseased')
interventions = {'AS', 'AP', 'SP', 'ASAP', 'SPAS', 'SPAP'};
force_data = {AS_data, AP_data, SP_data, ASAP_data, SPAS_data, SPAP_data};

sig_interventions = {...
    '','***','***','','','*','','***';... % AS
    '*','','','','','','','***';... % AP
    '***','','','','***','*','','***';... % SP
    '','***','***','','','*','','';... % ASAP
    '','***','***','','***','***','','***';... % SPAS
    '***','**','','','***','***','','***';... % SPAP
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
    % err = std_force_diff; % standard deviations
    err = std_force_diff./sqrt(numel(tests_present)); % standard errors
    err(isnan(err)) = 0;
    errorbar(x, avg_force_diff, err, 'k', 'LineWidth', 1.2, 'LineStyle', 'none');

    % Plot significance markers
    for pin = 1 : number_of_pins
        if ~isempty(sig_interventions{interv_idx,pin})
            text(pin,avg_force_diff(pin)+sign(avg_force_diff(pin))*(err(pin)+0.02),sig_interventions{interv_idx,pin},'HorizontalAlignment','center')
        end
    end
    
    % Formatting
    title([interventions{interv_idx},', n=',num2str(length(tests_present))], 'FontWeight', 'bold', 'FontSize', 12);
    
    % Only add x and y labels to subplot 4 (ASAP)
    if interv_idx == 5
        xlabel('Pin Number');
        % ylabel('\Delta F (N)');
    end
    ylabel('\Delta F (N)');
    
    set(gca, 'XTick', 1:number_of_pins);
    grid on;
    set(gca, 'FontSize', 10);
    ylim([-0.5,0.20])
    yticks(-0.5:0.1:0.2)
    hold off;
end

pubPlot('Width','double','Height',400,'Filename','ForceDiff_Pin','FileExtension',{'.png','.eps'});

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

pubPlot('Width','double','Height',400,'Filename','ForceDiff_Contour','FileExtension',{'.png','.eps'});



%% Export Force Differences Relative to Diseased Data

% Initialize empty structure array for collecting force difference rows
diff_rows = struct('Heart', {}, 'Intervention', {}, 'NumClips', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {}, ...
    'ForceDifference', {}, 'PressureDifference', {}, 'FlowDifference', {}, 'ClipOrder', {}, 'Treatment', {}, ...
    'AnnularArea', {}, 'AnnularPerimeter', {}, 'SLDiameter', {}, 'APDiameter', {}, 'CoaptationGapArea', {});
diff_rows_full = struct('Heart', {}, 'Intervention', {}, 'NumClips', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {}, ...
    'ForceDifference', {}, 'PressureDifference', {}, 'FlowDifference', {}, 'ClipOrder', {}, 'Treatment', {}, 'ReplicateTestNum', {}, ...
    'AnnularArea', {}, 'AnnularPerimeter', {}, 'SLDiameter', {}, 'APDiameter', {}, 'CoaptationGapArea', {});
intervention_diff_rows = struct('Heart', {}, 'Intervention', {}, 'NumClips', {}, 'Pressure', {}, 'FlowRate', {}, ...
    'AcrossAxisForceDiff', {},'SLForce', {}, 'APForce', {}, 'ClipDistanceMax', {}, 'ClipDistanceMin', {}, 'ClipDistanceMean', {},...
     'ClipAngle', {}, 'PeakForceAngle', {}, 'ForceClipAngleDifference', {}, 'ClipOrder', {}, 'Treatment', {}, ...
    'AnnularArea', {}, 'AnnularPerimeter', {}, 'SLDiameter', {}, 'APDiameter', {}, 'CoaptationGapArea', {});
morph_full_diff_rows = struct('Heart', {}, 'Intervention', {}, 'NumClips', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {}, ...
    'ForceDifference', {}, 'PressureDifference', {}, 'FlowDifference', {}, 'ClipOrder', {}, 'Treatment', {}, ...
    'AnnularArea', {}, 'AnnularPerimeter', {}, 'SLDiameter', {}, 'APDiameter', {}, 'CoaptationGapArea', {}, ...
    'TotalLeafletArea', {}, 'TotalLeafletPerimeter', {});
morph_leaflet_diff_rows = struct('Heart', {}, 'Intervention', {}, 'NumClips', {}, 'Pressure', {}, 'FlowRate', {}, 'Pin', {}, 'Force', {}, ...
    'ForceDifference', {}, 'PressureDifference', {}, 'FlowDifference', {}, 'ClipOrder', {}, 'Treatment', {}, ...
    'AnnularArea', {}, 'AnnularPerimeter', {}, 'SLDiameter', {}, 'APDiameter', {}, 'CoaptationGapArea', {}, ...
    'Leaflet', {}, 'LeafletArea', {}, 'LeafletPerimeter', {}, 'LeafletHeight', {}, 'LeafletWidth', {});

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

        % number of clips
        if strcmp(intervention_names{i},'AS') | strcmp(intervention_names{i},'AP') | strcmp(intervention_names{i},'SP')
            numClips = 1;
        elseif strcmp(intervention_names{i},'ASAP') | strcmp(intervention_names{i},'SPAS') | strcmp(intervention_names{i},'SPAP')
            numClips = 2;
        else
            numClips = 0;
        end

        % clip order and treatment
        if numClips == 0
            ClipOrder = 'None';
            Treatment = 0;
        elseif numClips == 2
            ClipOrder = 'Multiple';
            Treatment = test_procedure{tt};
        else
            ClipOrder = find(strcmp(test_order{tt},intervention_names{i}));
            if ClipOrder == 1
                ClipOrder = 'First';
            elseif ClipOrder == 3
                ClipOrder = 'Second';
            elseif ClipOrder == 4
                ClipOrder = 'Second';
            elseif ClipOrder == 5
                ClipOrder = 'Third';
            end
            Treatment = test_procedure{tt};
        end

        % annotation data
        annotationData = eval(['test' num2str(tt) 'Annotations']);
        annotationTestDex = 0;

        if strcmp(intervention_names{i},'Diseased')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,'10-') || contains(annotationData(j).Test,'11-') || contains(annotationData(j).Test,'12-')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'ASAP')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,'ASAP')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'SPAS')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,'SPAS')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'SPAP')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,'SPAP')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'AP')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,' AP ')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'AS')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,' AS ')
                    annotationTestDex = j;
                end
            end
        elseif strcmp(intervention_names{i},'SP')
            for j = 1 : numel(annotationData)
                if contains(annotationData(j).Test,' SP ')
                    annotationTestDex = j;
                end
            end
        end

        a = annotationData(annotationTestDex);
        V = a.AnnulusOutline.Vertices;

        AnnularArea = a.Area;
        AnnularPerimeter = a.Perimeter;
        SL_Diameter = a.SL_Diameter;
        AP_Diameter = a.AP_Diameter;
        Coaptation_Gap_Area = a.Gap_Area;

        % find nearest pin to clip line
        closest_pins = [];
        closest_pins_clip_2 = [];
        if ~isempty(annotationData(annotationTestDex).ClipAxis_Idx)
            closest_pins_1 = knnsearch(V(a.Pin_Idx,:),...
                V(a.ClipAxis_Idx(1),:),'K',2);
            closest_pins_2 = knnsearch(V(a.Pin_Idx,:),...
                V(a.ClipAxis_Idx(2),:),'K',2);
            closest_pins = [closest_pins_1, closest_pins_2];
            closest_pins = unique(closest_pins);
        end
        if ~isempty(a.ClipAxis2_Idx)
            closest_pins_3 = knnsearch(V(a.Pin_Idx,:),...
                V(a.ClipAxis2_Idx(1),:),'K',2);
            closest_pins_4 = knnsearch(V(a.Pin_Idx,:),...
                V(a.ClipAxis2_Idx(2),:),'K',2);
            closest_pins_clip_2 = [closest_pins_3, closest_pins_4];
            closest_pins_clip_2 = unique(closest_pins_clip_2);
        end
        
        % assume pins measure radial forces towards centroid
        pin_vecs = -V(a.Pin_Idx,:)./...
            vecnorm(V(a.Pin_Idx,:),2,2);

        % find distances from clip to annulus
        clip_distance_max = 0;
        clip_distance_min = 0;
        clip_distance_mean = 0;
        clip_2_distance_max = 0;
        clip_2_distance_min = 0;
        clip_2_distance_mean = 0;
        if ~isempty(a.ClipAxis_Idx)
            clip_distances = vecnorm([V(a.ClipAxis_Idx(1),:) - a.ClipCenter; ...
                V(a.ClipAxis_Idx(2),:) - a.ClipCenter],2,2);
            clip_distances = a.mmPerPixel.*clip_distances;
            clip_distance_max = max(clip_distances);
            clip_distance_min = min(clip_distances);
            clip_distance_mean = (max(clip_distances) + min(clip_distances))/2;
        end
        if ~isempty(a.ClipAxis2_Idx)
            clip2_distances = vecnorm([V(a.ClipAxis2_Idx(1),:) - a.ClipCenter2; ...
                V(a.ClipAxis2_Idx(2),:) - a.ClipCenter2],2,2);
            clip2_distances = a.mmPerPixel.*clip2_distances;
            clip_2_distance_max = max(clip2_distances);
            clip_2_distance_min = min(clip2_distances);
            clip_2_distance_mean = (max(clip2_distances) + min(clip2_distances))/2;
        end


        % morphology data
        if tt >= 8
            MorphFullDex = tt - 7;
            MorphLeafletDex = [3*(tt-8)+1:3*(tt-8)+3];
        end

        
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
        interv_force_mean = zeros(1,8);
        diseased_force_mean = zeros(1,8);
        force_diff = zeros(1,8);
        for pin = 1:8
            force_col = ['Force' num2str(pin)];
            
            % Calculate mean forces
            interv_force_mean(pin) = mean([interv_data.(force_col)]);
            diseased_force_mean(pin) = mean([diseased_ref.(force_col)]);
            
            % Calculate force difference (intervention - diseased)
            force_diff(pin) = interv_force_mean(pin) - diseased_force_mean(pin);

            % switch "Diseased" to "Control"
            if strcmp(intervention_names{i},'Diseased')
                treatment_name = 'Control';
            else
                treatment_name = intervention_names{i};
            end

            % do not average techincal replicates
            for replicate = 1 : length(interv_data)
                diff_rows_full(end+1) = struct( ...
                    'Heart', tt, ...
                    'Intervention', treatment_name, ...
                    'NumClips', numClips,...
                    'Pressure', interv_data(replicate).Pressure, ...
                    'FlowRate', interv_data(replicate).FlowRate, ...
                    'Pin', pin, ...
                    'Force', interv_data(replicate).(force_col), ...
                    'ForceDifference', interv_data(replicate).(force_col) - diseased_force_mean(pin), ...
                    'PressureDifference', interv_data(replicate).Pressure - diseased_pressure_mean, ...
                    'FlowDifference', interv_data(replicate).FlowRate - diseased_flow_mean, ...
                    'ClipOrder', ClipOrder, ...
                    'Treatment', Treatment, ...
                    'ReplicateTestNum', replicate, ...
                    'AnnularArea', AnnularArea, ...
                    'AnnularPerimeter', AnnularPerimeter, ...
                    'SLDiameter', SL_Diameter, ...
                    'APDiameter', AP_Diameter, ...
                    'CoaptationGapArea', Coaptation_Gap_Area ...
                );
            end
            
            % average technical replicates
            diff_rows(end+1) = struct( ...
                'Heart', tt, ...
                'Intervention', treatment_name, ...
                'NumClips', numClips,...
                'Pressure', interv_pressure_mean, ...
                'FlowRate', interv_flow_mean, ...
                'Pin', pin, ...
                'Force', interv_force_mean(pin), ...
                'ForceDifference', force_diff(pin), ...
                'PressureDifference', pressure_diff, ...
                'FlowDifference', flow_diff, ...
                'ClipOrder', ClipOrder, ...
                'Treatment', Treatment, ...
                'AnnularArea', AnnularArea, ...
                'AnnularPerimeter', AnnularPerimeter, ...
                'SLDiameter', SL_Diameter, ...
                'APDiameter', AP_Diameter, ...
                'CoaptationGapArea', Coaptation_Gap_Area ...
            );

            % morphological data if present
            if tt > 8
                morph_full_diff_rows(end+1) = struct( ...
                    'Heart', tt, ...
                    'Intervention', treatment_name, ...
                    'NumClips', numClips,...
                    'Pressure', interv_pressure_mean, ...
                    'FlowRate', interv_flow_mean, ...
                    'Pin', pin, ...
                    'Force', interv_force_mean(pin), ...
                    'ForceDifference', force_diff(pin), ...
                    'PressureDifference', pressure_diff, ...
                    'FlowDifference', flow_diff, ...
                    'ClipOrder', ClipOrder, ...
                    'Treatment', Treatment, ...
                    'AnnularArea', AnnularArea, ...
                    'AnnularPerimeter', AnnularPerimeter, ...
                    'SLDiameter', SL_Diameter, ...
                    'APDiameter', AP_Diameter, ...
                    'CoaptationGapArea', Coaptation_Gap_Area, ...
                    'TotalLeafletArea', MorphologyFull(MorphFullDex).LeafletArea, ...
                    'TotalLeafletPerimeter', MorphologyFull(MorphFullDex).LeafletPerimeter ...
                );

                for l = 1 : 3
                    morph_leaflet_diff_rows(end+1) = struct( ...
                        'Heart', tt, ...
                        'Intervention', treatment_name, ...
                        'NumClips', numClips,...
                        'Pressure', interv_pressure_mean, ...
                        'FlowRate', interv_flow_mean, ...
                        'Pin', pin, ...
                        'Force', interv_force_mean(pin), ...
                        'ForceDifference', force_diff(pin), ...
                        'PressureDifference', pressure_diff, ...
                        'FlowDifference', flow_diff, ...
                        'ClipOrder', ClipOrder, ...
                        'Treatment', Treatment, ...
                        'AnnularArea', AnnularArea, ...
                        'AnnularPerimeter', AnnularPerimeter, ...
                        'SLDiameter', SL_Diameter, ...
                        'APDiameter', AP_Diameter, ...
                        'CoaptationGapArea', Coaptation_Gap_Area, ...
                        'Leaflet', MorphologyLeaflet(MorphLeafletDex(l)).LeafletID, ...
                        'LeafletArea', MorphologyLeaflet(MorphLeafletDex(l)).LeafletArea, ...
                        'LeafletPerimeter', MorphologyLeaflet(MorphLeafletDex(l)).LeafletPerimeter, ...
                        'LeafletHeight', MorphologyLeaflet(MorphLeafletDex(l)).LeafletHeight, ...
                        'LeafletWidth', MorphologyLeaflet(MorphLeafletDex(l)).LeafletWidth ...
                    );
                end
            end

        end

        % Calculate force difference across the clip axis
        % Note: Watch scale on this. Summing different numbers of "closest
        % pins" can scale force F reported to 2F, 4F, etc.

        % AcrossAxisForceDiff = 0;
        % % AcrossAxisForceDiff = AcrossAxisForceDiff + mean([force_diff(closest_pins)]);
        % % AcrossAxisForceDiff2 = AcrossAxisForceDiff + mean([force_diff(closest_pins_clip_2)]);
        % AcrossAxisForceDiff = AcrossAxisForceDiff + sum([force_diff(closest_pins)]);
        % AcrossAxisForceDiff2 = AcrossAxisForceDiff + sum([force_diff(closest_pins_clip_2)]);
        % if isnan(AcrossAxisForceDiff)
        %     AcrossAxisForceDiff = 0;
        % end

        % Interpolate force differences at the two annulus/clip-axis crossings.
        AcrossAxisForceDiff  = NaN;
        AcrossAxisForceDiff2 = NaN;
        crossingForces1 = [];
        crossingForces2 = [];

        if ~isempty(a.ClipAxis_Idx)
            crossingForces1 = interpolateAnnularForce( ...
                V, a.Pin_Idx, force_diff, a.ClipAxis_Idx);
            AcrossAxisForceDiff = sum(crossingForces1);
        end

        if ~isempty(a.ClipAxis2_Idx)
            crossingForces2 = interpolateAnnularForce( ...
                V, a.Pin_Idx, force_diff, a.ClipAxis2_Idx);
            AcrossAxisForceDiff2 = sum(crossingForces2);
        end

        % Compute radial force vectors
        pin_force_vecs = force_diff'.*pin_vecs;
        
        % sum force, which would be like total force vector? Don't think it
        % actually means anything
        sum_rad_force = sum(pin_force_vecs);

        % SL & AP force: 
        SL_force = 0;
        AP_force = 0;
        for pin = 1:8
            if pin_vecs(pin,1) < 0 && pin_vecs(pin,2) < 0 % Q1
                SL_force = SL_force - pin_force_vecs(pin,1);
                AP_force = AP_force - pin_force_vecs(pin,2);
            elseif pin_vecs(pin,1) >= 0 && pin_vecs(pin,2) < 0 % Q2
                SL_force = SL_force + pin_force_vecs(pin,1);
                AP_force = AP_force - pin_force_vecs(pin,2);
            elseif pin_vecs(pin,1) >= 0 && pin_vecs(pin,2) >= 0 % Q3
                SL_force = SL_force + pin_force_vecs(pin,1);
                AP_force = AP_force + pin_force_vecs(pin,2);
            elseif pin_vecs(pin,1) < 0 && pin_vecs(pin,2) >= 0 % Q4
                SL_force = SL_force - pin_force_vecs(pin,1);
                AP_force = AP_force + pin_force_vecs(pin,2);
            end
        end

        % Clip Orientation
        % Clip-axis angles relative to AP (+y), modulo 180 degrees.
        clipIdx = {a.ClipAxis_Idx, a.ClipAxis2_Idx};
        clipAngles = nan(1,2);

        for c = 1:2
            if numel(clipIdx{c}) == 2
                d = V(clipIdx{c}(2),:) - V(clipIdx{c}(1),:);
                clipAngles(c) = mod(atan2d(d(1), d(2)), 180);
            end
        end

        % Combined orientation: axial circular mean of available clips.
        validAngles = clipAngles(isfinite(clipAngles));
        ClipAngle = NaN;
        if ~isempty(validAngles)
            z = mean(exp(2i * deg2rad(validAngles)));
            if abs(z) > 1e-8  % Perpendicular axes have no unique mean orientation.
                ClipAngle = mod(rad2deg(angle(z))/2, 180);
            end
        end

        % Most negative pin force CHANGE; use interv_force_mean for measured force.
        F = force_diff(:);
        F(~isfinite(F)) = Inf;
        [peakNegativeForce, peakPin] = min(F);

        PeakForceAngle = NaN;
        if isfinite(peakNegativeForce)
            r = V(a.Pin_Idx(peakPin),:);  % Centroid-to-pin radial vector.
            PeakForceAngle = mod(atan2d(r(1), r(2)), 180);
        else
            peakNegativeForce = NaN;
            peakPin = NaN;
        end

        % Smallest angular separation between force and clip axes: 0–90 degrees.
        ForceClipAngleDifference = abs( ...
            mod(PeakForceAngle - ClipAngle + 90, 180) - 90);

        % store in struct
        intervention_diff_rows(end+1) = struct( ...
                'Heart', tt, ...
                'Intervention', treatment_name, ...
                'NumClips', numClips,...
                'Pressure', interv_pressure_mean, ...
                'FlowRate', interv_flow_mean, ...
                'AcrossAxisForceDiff', AcrossAxisForceDiff, ...
                'SLForce', SL_force, ...
                'APForce', AP_force, ...
                'ClipDistanceMax', clip_distance_max, ...
                'ClipDistanceMin', clip_distance_min, ...
                'ClipDistanceMean', clip_distance_mean, ...
                'ClipAngle',ClipAngle, ...
                'PeakForceAngle',PeakForceAngle, ...
                'ForceClipAngleDifference',ForceClipAngleDifference, ...
                'ClipOrder', ClipOrder, ...
                'Treatment', Treatment, ...
                'AnnularArea', AnnularArea, ...
                'AnnularPerimeter', AnnularPerimeter, ...
                'SLDiameter', SL_Diameter, ...
                'APDiameter', AP_Diameter, ...
                'CoaptationGapArea', Coaptation_Gap_Area ...
            );

        if ~isempty(annotationData(annotationTestDex).ClipAxis2_Idx)
            % add another row
            intervention_diff_rows(end+1) = struct( ...
                'Heart', tt, ...
                'Intervention', treatment_name, ...
                'NumClips', numClips,...
                'Pressure', interv_pressure_mean, ...
                'FlowRate', interv_flow_mean, ...
                'AcrossAxisForceDiff', AcrossAxisForceDiff2, ...
                'SLForce', SL_force, ...
                'APForce', AP_force, ...
                'ClipDistanceMax', clip_2_distance_max, ...
                'ClipDistanceMin', clip_2_distance_min, ...
                'ClipDistanceMean', clip_2_distance_mean, ...
                'ClipAngle',ClipAngle, ...
                'PeakForceAngle',PeakForceAngle, ...
                'ForceClipAngleDifference',ForceClipAngleDifference, ...
                'ClipOrder', ClipOrder, ...
                'Treatment', Treatment, ...
                'AnnularArea', AnnularArea, ...
                'AnnularPerimeter', AnnularPerimeter, ...
                'SLDiameter', SL_Diameter, ...
                'APDiameter', AP_Diameter, ...
                'CoaptationGapArea', Coaptation_Gap_Area ...
            );
        end

    end
end

% Convert to table and write to Excel
if ~isempty(diff_rows_full)
    T_diff = struct2table(diff_rows_full);
    csvFile_diff = fullfile(projectRoot, 'TriClipXT_Full_Statistics.csv');
    writetable(T_diff, csvFile_diff, 'WriteMode', 'overwrite');
    fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
else
    warning('No force difference rows collected for export; no file written.');
end

if ~isempty(diff_rows)
    T_diff = struct2table(diff_rows);
    csvFile_diff = fullfile(projectRoot, 'TriClipXT_Statistics.csv');
    writetable(T_diff, csvFile_diff, 'WriteMode', 'overwrite');
    fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
else
    warning('No force difference rows collected for export; no file written.');
end

if ~isempty(intervention_diff_rows)
    T_diff = struct2table(intervention_diff_rows);
    csvFile_diff = fullfile(projectRoot, 'TriClipXT_Intervention_Statistics.csv');
    writetable(T_diff, csvFile_diff, 'WriteMode', 'overwrite');
    fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
else
    warning('No force difference rows collected for export; no file written.');
end

if ~isempty(morph_full_diff_rows)
    T_diff = struct2table(morph_full_diff_rows);
    csvFile_diff = fullfile(projectRoot, 'TriClipXT_Full_Morphology_Statistics.csv');
    writetable(T_diff, csvFile_diff, 'WriteMode', 'overwrite');
    fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
else
    warning('No force difference rows collected for export; no file written.');
end

if ~isempty(morph_leaflet_diff_rows)
    T_diff = struct2table(morph_leaflet_diff_rows);
    csvFile_diff = fullfile(projectRoot, 'TriClipXT_Leaflet_Morphology_Statistics.csv');
    writetable(T_diff, csvFile_diff, 'WriteMode', 'overwrite');
    fprintf('Exported %d force difference rows to %s\n', height(T_diff), csvFile_diff);
else
    warning('No force difference rows collected for export; no file written.');
end


%% Plot and save all 10 
testDataCurrent = test13data;
dataNum = 13;
run_dex = [10:21];

% startdex = [];
% for i = 10:length(testDataCurrent)
% startdex(i) = find(test12data{i}.Pressure(:) > 1.2*mean(test12data{i}.Pressure(1:10)),1)
% end
% enddex = [];
% for i = 10:length(testDataCurrent)
% temp = find(test12data{i}.Pressure(:) > 1.2*mean(test12data{i}.Pressure(end-10:end)));
% enddex(i) = temp(end);
% end

if dataNum == 2
    % Test 2 - SP,  AP,     SPAP
    start_time_dex = [42,37,48,...
        58,58,51,...
        47,42,46,...
        40,52,55];
    end_time_dex = [327,328,323,...
        348,355,219,...
        342,333,206,...
        373,406,221];
    end_time_dex = [test2idx(run_dex).segmentEndIdx];
    start_time_dex = start_time_dex - 30;
    data_Colors = [COLORS.Diseased; COLORS.SP; COLORS.AP; COLORS.SPAP];
elseif dataNum == 3
    % Test 3 - AS,  SP,     SPAS
    start_time_dex = [72,56,92,...
        40,46,57,...
        39,47,106,...
        42,43,73];
    end_time_dex = [340,328,374,...
        348,337,293,...
        345,358,330,...
        343,344,315];
    end_time_dex = [test3idx(run_dex).segmentEndIdx];
    start_time_dex = start_time_dex - 30;
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.SP; COLORS.SPAS];
elseif dataNum == 4
    % Test 4 - AS,  AP,     ASAP
    % disease 10:12, AS 13:15, AP 16:18, ASAP  19:21
    start_time_dex = [40,34,80,...
        40,37,69,...
        50,43,45,...
        41,35,97];
    end_time_dex = [342,348,318,...
        350,356,283,...
        357,365,314,...
        362,367,332];
    end_time_dex = [test4idx(run_dex).segmentEndIdx];
    start_time_dex = start_time_dex - 30;
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.AP; COLORS.ASAP];
elseif dataNum == 5
    % Test 5 - AS,  AP,     ASAP
    % disease 10:12, AS 13:15, AP 16:18, ASAP  19:21
    start_time_dex = [52,44,51,...
        47,51,42,...
        46,45,52,...
        52,57,42];
    end_time_dex = [377,366,204,...
        357,350,212,...
        392,421,272,...
        476,450,275];
    start_time_dex = start_time_dex - 30;
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.AP; COLORS.ASAP];
elseif dataNum == 8
    % Test 8 -  SP,  SPAS,   AS 
    % disease 10:12, SP 13:15, SPAS 16:18, AS  19:21
    start_time_dex = [46,43,73,...
        58,50,65,...
        47,36,41,...
        45,46,67];
    end_time_dex = [348,359,478,...
        408,350,463,...
        396,397,410,...
        417,368,425];
    start_time_dex = start_time_dex - 30;
    data_Colors = [COLORS.Diseased; COLORS.SP; COLORS.SPAS; COLORS.AS];
elseif dataNum == 9
    % Test 9 -  SP,  SPAP,   AP,    ASAP,   AS
    run_dex = [10:27];
    start_time_dex = [35,44,50,...
        32,49,35,...
        45,44,46,...
        44,41,72,...
        29,40,51,...
        47,45,41];
    end_time_dex = [404,355,280,...
        386,400,323,...
        455,386,334,...
        404,348,363,...
        367,372,293,...
        379,389,313];
    start_time_dex = start_time_dex - 20;
    data_Colors = [COLORS.Diseased; COLORS.SP; COLORS.SPAP; COLORS.AP; COLORS.ASAP; COLORS.AS];
elseif dataNum == 10
    % Test 10 - AS,  ASAP,   AP,    SPAP,   SP 
    run_dex = [10:27];
    start_time_dex = [36,36,46,...
        41,40,73,...
        38,38,87,...
        28,42,80,...
        27,38,46,...
        39,42,52];
    end_time_dex = [319,347,289,...
        320,323,323,...
        317,330,322,...
        319,344,308,...
        382,335,353,...
        349,290,328];
    start_time_dex = start_time_dex - 20;
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.ASAP; COLORS.AP; COLORS.SPAP; COLORS.SP];
elseif dataNum == 11
    % Test 8 -  SP,  SPAS,   AS
    % disease 10:12, SP 13:15, SPAS 16:18, AS  19:21
    start_time_dex = [49,38, 70,...
        42, 42, 46,...
        53, 23, 42,...
        43, 46, 45];
    end_time_dex = [340, 342, 292,...
        322, 325, 361,...
        332, 340, 326,...
        324, 332, 304];
    start_time_dex = start_time_dex - 20;
    data_Colors = [COLORS.Diseased; COLORS.SP; COLORS.SPAS; COLORS.AS];
elseif dataNum == 12
    % Test 12 - AS,  ASAP,   AP,    SPAP,   SP 
    run_dex = [10:27];
    start_time_dex = [38, 36, 47,...
        33, 46, 54,...
    	41, 41, 58,...
        44,	37,	186,...
        45,	40,	54,...
    	43,	41,	78];
    end_time_dex = [300, 311, 281,...
    	354, 317, 288,...
        316, 319, 291,...
    	315, 328, 448,...
    	354, 323, 296,...
        357, 351, 384];
    start_time_dex = start_time_dex - 20;
    end_time_dex = end_time_dex + 10;
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.ASAP; COLORS.AP; COLORS.SPAP; COLORS.SP];
elseif dataNum == 13
    % Test 13 - AS,  SPAS,   SP,    SPAP,   AP 
    run_dex = [10:27];
    start_time_dex = [41, 48, 42,...
        40, 37, 72,...
        39, 36, 69,...
    	35, 49, 68,...
        51, 39, 68,...
        43, 39, 62];
    end_time_dex = [606, 433, 538,...
        365, 373, 379,...
        435, 417, 397,...
        315, 342, 361,...
        458, 381, 333,...
        357, 335, 335];
    start_time_dex = start_time_dex - 20;
    end_time_dex = end_time_dex + 10;
    % for i = 1 : length(run_dex)
    %     start_time_dex(i) = test13idx(9+i).segmentStartIdx;
    %     end_time_dex(i) = test13idx(9+i).segmentEndIdx;
    % end
    data_Colors = [COLORS.Diseased; COLORS.AS; COLORS.SPAS; COLORS.SP; COLORS.SPAP; COLORS.AP];
end

% resample data
num_pts = 1000;
time_norm = linspace(0,1,num_pts);
resampled_data = zeros(11,num_pts,3);
resampled_data(1,:,:) = repmat(time_norm,1,1,size(resampled_data,3));
% 9 x num_pts x numel(run_dex)
% Row 1: Normalized Time
% Row 2: Pressure
% Row 3: Flow
% Row 4-11: Pin Forces 1-8 in order
% 3rd dimension is run
for i = 1 : numel(run_dex)
    % pressure
    time_count = ([1:length(testDataCurrent{run_dex(i)}.Pressure(start_time_dex(i):end_time_dex(i)))]-1)./(length(testDataCurrent{run_dex(i)}.Pressure(start_time_dex(i):end_time_dex(i)))-1);
    pressure = testDataCurrent{run_dex(i)}.Pressure(start_time_dex(i):end_time_dex(i))-testDataCurrent{run_dex(i)}.Pressure(1);
    flow = testDataCurrent{run_dex(i)}.FlowRate_ml_s_(start_time_dex(i):end_time_dex(i))-testDataCurrent{run_dex(i)}.FlowRate_ml_s_(1);
    resampled_data(2,:,i) = interp1(time_count,pressure,time_norm);
    resampled_data(3,:,i) = interp1(time_count,flow,time_norm);
    for j = 1:8
        forceCol = ['Force' num2str(j)];
        resampled_data(j+3,:,i) = interp1(time_count,testDataCurrent{run_dex(i)}.(forceCol)(start_time_dex(i):end_time_dex(i))-testDataCurrent{run_dex(i)}.(forceCol)(1),time_norm);
    end
end

resampled_avg = zeros(size(resampled_data));
resampled_std = zeros(size(resampled_data));
% take averages and stds at each point
for i = 1 : num_pts
    for j = 1 : numel(run_dex)/3
        for k = 2 : size(resampled_data,1)
            resampled_avg(k,i,j) = mean(resampled_data(k,i,3*j-2:3*j));
            resampled_std(k,i,j) =  std(resampled_data(k,i,3*j-2:3*j));
        end
    end
end

% subplot_dex = [1,7,3,4,5,6,9,10,11,12];
subplot_dex = [1,6,2,3,4,5,7,8,9,10];
figure()
for j = 1 : numel(run_dex)/3
    for k = 2 : size(resampled_data,1)
        subplot(2,5,subplot_dex(k-1))
        hold on
        boundedline(time_norm, resampled_avg(k,:,j), resampled_std(k,:,j),'cmap',data_Colors(j,:),'alpha');
        xticks([0,1])
        xticklabels({'0','1'})
        xtickangle(0)
        xlabel('Normalized Time (-)')
        if k-1 == 1
            ylabel('Pressure (mmHg)')
            if dataNum == 3
                ylim([0,40])
                yticks(0:10:40)
            else
                ylim([0,40])
                yticks(0:10:40)
            end
        elseif k-1 == 2
            ylabel('Flow Rate (mL/s)')
            if dataNum == 3
                ylim([0,60])
                yticks(0:10:60)
            elseif dataNum == 8 || 12
                ylim([0,70])
                yticks(0:10:70)
            else
                ylim([0,60])
                yticks(0:10:60)
            end
        else
            if dataNum == 3
                ylim([-0.5,0.2])
                yticks(-0.5:0.1:0.2)
            elseif dataNum == 8
                if k-1 == 5
                    legend('','Diseased','','SP','','SPAS','','AS','')
                end
                ylim([-1,0.5])
                yticks(-1:0.25:0.5)
            elseif dataNum == 9
                if k-1 == 5
                    legend('','Diseased','','SP','','SPAP','','AP','','ASAP','','AS')
                end
                ylim([-1,0.5])
                yticks(-1:0.25:0.5)
            elseif dataNum == 10
                if k-1 == 5
                    legend('','Diseased','','AS','','ASAP','','AP','','SPAP','','SP')
                end
                ylim([-1,0.5])
                yticks(-1:0.25:0.5)
            elseif dataNum == 11
                if k-1 == 5
                    legend('','Diseased','','SP','','SPAS','','AS','')
                end
                ylim([-.75,0.25])
                yticks(-0.75:0.25:0.25)
            elseif dataNum == 12
                if k-1 == 5
                    legend('','Diseased','','AS','','ASAP','','AP','','SPAP','','SP')
                end
                ylim([-1,0.5])
                yticks(-1:0.25:0.5)
            elseif dataNum == 13
                if k-1 == 5
                    legend('','Diseased','','AS','','SPAS','','SP','','SPAP','','AP')
                end
                ylim([-1,0.5])
                yticks(-1:0.25:0.5)
            else
                ylim([-0.5,0.2])
                yticks(-0.5:0.1:0.2)
            end
            if k-1 == 3 || k-1 == 7
                ylabel('\Delta F (N)');
            end
        end
    end
end
% subplot(2,6,2)
% plot(time_norm,time_norm)
% axis off
% subplot(2,6,8)
% plot(time_norm,time_norm)
% axis off

% if 0
%     pubPlot('Width','double','Height',400,'Filename',['Raw_Data_Test',num2str(dataNum)],'FileExtension',{'.png','.eps'});
% else
%     pubPlot('Width','double','Height',400);
% end


%% Plot from long form
p = plotInternventionLongForm(intervention_diff_rows,'AcrossAxisForceDiff');
p = plotInternventionLongForm(intervention_diff_rows,'SLForce');
p = plotInternventionLongForm(intervention_diff_rows,'APForce');
p = plotInternventionLongForm(intervention_diff_rows,'SumRadialForceSL');
p = plotInternventionLongForm(intervention_diff_rows,'SumRadialForceAP');

%% HELPER FUNCTIONS
function Fcross = interpolateAnnularForce(V, pinIdx, Fpin, crossingIdx)
% Periodic linear interpolation using distance along a closed annulus.
% V:           ordered annular boundary vertices, N-by-2
% pinIdx:      vertex index for each pin, in force-channel order
% Fpin:        force or force difference for each pin
% crossingIdx: vertex indices of the two clip-axis crossings

    pinIdx = pinIdx(:);
    Fpin = Fpin(:);
    crossingIdx = crossingIdx(:);

    assert(size(V,2) == 2 && all(isfinite(V(:))), ...
        'Expected one continuous annular boundary without NaN separators.');
    assert(numel(pinIdx) == numel(Fpin), ...
        'Each pin must have a corresponding force value.');
    assert(numel(crossingIdx) == 2, ...
        'Expected exactly two clip-axis crossings.');

    % Cumulative arc length at each vertex, including closing edge.
    edgeLength = vecnorm(diff([V; V(1,:)], 1, 1), 2, 2);
    s = [0; cumsum(edgeLength(1:end-1))];
    perimeterLength = sum(edgeLength);

    assert(perimeterLength > 0, 'Annular perimeter must be positive.');

    % Sort pins by position along the boundary, preserving force pairing.
    sPin = mod(s(pinIdx), perimeterLength);
    [sPin, order] = sort(sPin);
    FpinSorted = Fpin(order);

    assert(numel(sPin) >= 2 && all(diff(sPin) > 0), ...
        'Pins must occupy distinct positions along the annulus.');

    % Extend periodically so interpolation also works across the seam
    % between the last and first boundary vertices.
    sExtended = [sPin(end) - perimeterLength; ...
                 sPin; ...
                 sPin(1) + perimeterLength];
    FExtended = [FpinSorted(end); FpinSorted; FpinSorted(1)];

    sCross = mod(s(crossingIdx), perimeterLength);
    Fcross = interp1(sExtended, FExtended, sCross, 'linear');

    % % plot for checking
    % figure()
    % hold on
    % axis equal
    % 
    % plot([V(:,1); V(1,1)], [V(:,2); V(1,2)], 'k-');
    % 
    % P = V(pinIdx,:);
    % C = V(crossingIdx,:);
    % 
    % plot(P(:,1), P(:,2), 'co', 'MarkerFaceColor', 'c');
    % plot(C(:,1), C(:,2), 'mo--', 'LineWidth', 1.5);
    % 
    % for i = 1:numel(pinIdx)
    %     text(P(i,1), P(i,2), ...
    %         sprintf('  P%d: %.3f N', i, Fpin(i)), ...
    %         'Color', [0 0.5 0.5]);
    % end
    % 
    % for i = 1:numel(crossingIdx)
    %     text(C(i,1), C(i,2), ...
    %         sprintf('  Crossing %d: %.3f N', i, Fcross(i)), ...
    %         'Color', 'm');
    % end
    % 
    % % Radial directions assume the annulus is centered at the origin.
    % pin_vecs   = P ./ vecnorm(P, 2, 2);
    % cross_vecs = C ./ vecnorm(C, 2, 2);
    % 
    % vec_scale = 1e3;
    % quiver(P(:,1), P(:,2), ...
    %     vec_scale * Fpin(:) .* pin_vecs(:,1), ...
    %     vec_scale * Fpin(:) .* pin_vecs(:,2), ...
    %     0, 'Color', [0 0.6 0.6]);
    % 
    % quiver(C(:,1), C(:,2), ...
    %     vec_scale * Fcross(:) .* cross_vecs(:,1), ...
    %     vec_scale * Fcross(:) .* cross_vecs(:,2), ...
    %     0, 'Color', 'm');

end

function fig = plotInternventionLongForm(intervention_diff_rows,DataName)
    fig = figure('Name',[sprintf(DataName),' vs Interventions']);
    hold on

    bar([0:6],[...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''Control'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AS'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''ASAP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAS'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAP'')).',sprintf(DataName),']'])),...
        ])
    
    scatter(0*ones(sum(strcmp({intervention_diff_rows.Intervention},'Control')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''Control'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(1*ones(sum(strcmp({intervention_diff_rows.Intervention},'AS')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AS'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(2*ones(sum(strcmp({intervention_diff_rows.Intervention},'AP')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AP'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(3*ones(sum(strcmp({intervention_diff_rows.Intervention},'SP')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SP'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(4*ones(sum(strcmp({intervention_diff_rows.Intervention},'ASAP')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''ASAP'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(5*ones(sum(strcmp({intervention_diff_rows.Intervention},'SPAS')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAS'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    scatter(6*ones(sum(strcmp({intervention_diff_rows.Intervention},'SPAP')),1),eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAP'')).',sprintf(DataName),']']),'MarkerEdgeColor','k')
    
    errorbar([0:6],[...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''Control'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AS'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''ASAP'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAS'')).',sprintf(DataName),']'])),...
        mean(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAP'')).',sprintf(DataName),']'])),...
        ],[...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''Control'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'Control'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AS'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'AS'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''AP'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'AP'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SP'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'SP'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''ASAP'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'ASAP'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAS'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'SPAS'))),...
        std(eval(['[intervention_diff_rows(strcmp({intervention_diff_rows.Intervention},''SPAP'')).',sprintf(DataName),']']))./sqrt(sum(strcmp({intervention_diff_rows.Intervention},'SPAP'))),...
        ], 'k', 'LineWidth', 1.2, 'LineStyle', 'none')
    
    xlim([0.5,6.5])
    xticks([0:6])
    xticklabels({'Dis', 'AS', 'AP', 'SP ', 'ASAP', 'SPAS', 'SPAP'})
    ylabel(DataName)
end