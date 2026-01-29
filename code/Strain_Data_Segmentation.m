%% Strain Gauge Analysis
% Clear workspace and command window
clc
clear all
close all

askUser = true;

%Before you begin, esnure that all excel files are saved as .csv, and that
%there are no duplicate column headers from the experiment starting and
%stopping. Use CTRL+F on each file and search for "pressure" as an easy
%way to find duplicates. For some reason, you may also need to to retype
%the "Time" column header in each file. 

% After running this file, you will get a 2D array containing the rlevant
% indices for each test. Rename the result as you wish and save as a .mat
% file for 

% Define your folder path containing CSV data files
% folderPath = 'C:\Users\trace\OneDrive\Documents\Rausch Lab\TriClip Experiment\Data\08_07_25\CSV Files';

% Get the full path of this script
scriptPath = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(scriptPath);
projectRoot = fileparts(scriptDir);

% Build data paths
dataDir = fullfile(projectRoot, 'data');

% ask user to select folder 
if askUser
    folderPath = uigetdir(dataDir);
else
    folderPath = fullfile(dataDir,"01_14_26/");
end

% Create pattern to identify all CSV files
filePattern = fullfile(folderPath, '*.csv');
% Get information about all matching CSV files
csvFiles = dir(filePattern);


% Initialize cell arrays to store data and filenames
allData = cell(length(csvFiles), 1);
names = cell(length(csvFiles), 1);

% Loop through each CSV file
for k = 1:length(allData)
    baseFileName = csvFiles(k).name;
    fullFileName = fullfile(folderPath, baseFileName);
    fprintf('Now reading %s\n', fullFileName);
    
    % Read CSV file into table
    data = readtable(fullFileName);
    names{k} = baseFileName;
    allData{k} = data;
end

% Create column vector of filenames and clean up
names_list = names.';
clear data


% Clean up time format - remove double periods
for i = 1:length(allData)
    disp(names_list(i))
    allData{i}.Time = strrep(allData{i}.Time, '..', '.');
end


% Calibration slopes for each pin
fcal = [
    -1621.4975  % pin1
    -1291.6308  % pin2
    -1368.1537  % pin3
    -1017.2041  % pin4
    -1106.5112  % pin5
    -1089.2077  % pin6
    -2948.5710  % pin7
    -1427.9047  % pin8
];

% Apply calibration slopes to each pin
for i = 1:length(allData)
    for j = 1:8
        allData{i}.(['Force' num2str(j)]) = allData{i}.(['Force' num2str(j)]) * fcal(j);
    end
    allData{i}.FlowRate_ml_s_ = allData{i}.FlowRate_ml_s_ * 200.44; % Convert to ml/s
    % Zero flow rate
    allData{i}.FlowRate_ml_s_ = allData{i}.FlowRate_ml_s_ - allData{i}.FlowRate_ml_s_(1);
end



% Create a new cell array for smoothed data (keeps original data intact)
smoothedData = cell(size(allData));

% Copy and smooth the data
for i = 1:length(allData)
    % Start with a copy of the original data
    smoothedData{i} = allData{i};
    
    % Smooth force data
    for j = 1:8
        smoothedData{i}.(['Force' num2str(j)]) = smoothdata(allData{i}.(['Force' num2str(j)]), 'gaussian', 20);
    end
    
    % Smooth pressure and flow rate
    smoothedData{i}.Pressure = smoothdata(allData{i}.Pressure, 'gaussian', 20);
    smoothedData{i}.FlowRate_ml_s_ = smoothdata(allData{i}.FlowRate_ml_s_, 'gaussian', 20);
end


% Initialize arrays to store region data for later use in force plots
stableRegionIdx = struct('maxPin', {}, 'maxIndex', {}, 'minIndex', {}, 'segmentStartIdx', {}, 'segmentEndIdx', {});

% Plot the derivatives of the smooth data - Minimalist version
for i=1:length(smoothedData)
    % Create simple integer time points
    timePoints = 1:height(smoothedData{i}.Time);

    %figure()
    %hold on
    
    % Colors for each pin
    pinColors = lines(8);
        
    % Find pin with greatest absolute derivative
    globalAbsMaxVal = 0;
    globalMaxPin = 0;
    allMaxIndices = zeros(8, 1);
    allMinIndices = zeros(8, 1);
    
    % Calculate and plot derivatives, find pin with max absolute derivative
    for j = 1:8
        zeroed_force = smoothedData{i}.(['Force' num2str(j)]) - smoothedData{i}.(['Force' num2str(j)])(1);
        derivative = gradient(zeroed_force);
        % Plot the derivative line
        %plot(timePoints, derivative, 'Color', pinColors(j,:));

        % Find max and min derivative values
        [maxVal, maxIdx] = max(derivative);
        [minVal, minIdx] = min(derivative);
        
        % Store indices
        allMaxIndices(j) = maxIdx;
        allMinIndices(j) = minIdx;
        
        % Calculate absolute maximum derivative
        absMax = max(abs(maxVal), abs(minVal));
        
        % Update if this pin has the greatest absolute derivative
        if absMax > globalAbsMaxVal
            globalAbsMaxVal = absMax;
            globalMaxPin = j;
        end
    end
        % Get indices of max and min points for selected pin
    maxIndex = allMaxIndices(globalMaxPin);
    minIndex = allMinIndices(globalMaxPin);
    
    % Identify segment between max and min with offset
    startOffsetPercent = 0.3; 
    endOffsetPercent = 0.1;
    startIdx = min(maxIndex, minIndex);
    endIdx = max(maxIndex, minIndex);
    segmentLength = endIdx - startIdx;
    segmentStartIdx = startIdx + round(segmentLength * startOffsetPercent);
    segmentEndIdx = endIdx - round(segmentLength * endOffsetPercent);

    % Save region data for use in force plots
    stableRegionIdx(i).maxPin = globalMaxPin;
    stableRegionIdx(i).maxIndex = maxIndex;
    stableRegionIdx(i).minIndex = minIndex;
    stableRegionIdx(i).segmentStartIdx = segmentStartIdx;
    stableRegionIdx(i).segmentEndIdx = segmentEndIdx;

    % Mark the region of interest
    %if segmentStartIdx < segmentEndIdx
        % Shade the region
        %x = [timePoints(segmentStartIdx), timePoints(segmentEndIdx), timePoints(segmentEndIdx), timePoints(segmentStartIdx)];
        %yLimits = ylim();
        %y = [yLimits(1), yLimits(1), yLimits(2), yLimits(2)];
        %patch(x, y, pinColors(globalMaxPin,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        % Add region label
        %text(mean([timePoints(segmentStartIdx), timePoints(segmentEndIdx)]), yLimits(2)*0.9, ...
        %    ['Region: ' num2str(segmentStartIdx) ' to ' num2str(segmentEndIdx)], ...
        %    'Color', pinColors(globalMaxPin,:), 'HorizontalAlignment', 'center');
    %end
    
    % Mark max and min points of selected pin
    %xline(timePoints(maxIndex), '--', 'LineWidth', 1.5, 'Color', pinColors(globalMaxPin,:));
    %xline(timePoints(minIndex), ':', 'LineWidth', 1.5, 'Color', pinColors(globalMaxPin,:));
    
    % Set plot attributes
    %title(['Force Derivatives - ' names{i}])
    %xlabel('Time Point')
    %ylabel('Force Derivative')
    %legend('Pin 1', 'Pin 2', 'Pin 3', 'Pin 4', 'Pin 5', 'Pin 6', 'Pin 7', 'Pin 8', 'Location', 'best')
end

disp('Done')

%% Plot each segment for validation
clc
close all


averageForces = zeros(length(allData), 8);
averageFlowRates = zeros(length(allData), 1);
averagePressures = zeros(length(allData), 1);

% Loop through each file and calculate average forces
for i = 1:length(allData)
        % Get the stable region indices
        segmentStartIdx = stableRegionIdx(i).segmentStartIdx;
        segmentEndIdx = stableRegionIdx(i).segmentEndIdx;
        
        % Calculate average forces for each pin
        for j = 1:8
            averageForces(i, j) = mean(allData{i}.(['Force' num2str(j)])(segmentStartIdx:segmentEndIdx));
        end
        
        % Calculate average flow rate and pressure
        averageFlowRates(i) = mean(allData{i}.FlowRate_ml_s_(segmentStartIdx:segmentEndIdx));
        averagePressures(i) = mean(allData{i}.Pressure(segmentStartIdx:segmentEndIdx));
end

% Plot
for i=1:length(allData)
    % Create simple integer time points
    timePoints = 1:height(allData{i}.Time);

    figure()
    hold on
    
    % Colors for each pin
    pinColors = lines(8);
    
    % Plot each force channel
    for j = 1:8
        zeroed_force = allData{i}.(['Force' num2str(j)]) - allData{i}.(['Force' num2str(j)])(1);
        plot(timePoints, zeroed_force, 'Color', pinColors(j,:));
    end
    % Get the pin and region indices
    maxPin = stableRegionIdx(i).maxPin;
    maxIndex = stableRegionIdx(i).maxIndex;
    minIndex = stableRegionIdx(i).minIndex;
    segmentStartIdx = stableRegionIdx(i).segmentStartIdx;
    segmentEndIdx = stableRegionIdx(i).segmentEndIdx;
    
    % Mark the region of interest
    if segmentStartIdx < segmentEndIdx
        % Shade the region
        x = [timePoints(segmentStartIdx), timePoints(segmentEndIdx), timePoints(segmentEndIdx), timePoints(segmentStartIdx)];
        yLimits = ylim();
        y = [yLimits(1), yLimits(1), yLimits(2), yLimits(2)];
        patch(x, y, pinColors(maxPin,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        % Add region label
        text(mean([timePoints(segmentStartIdx), timePoints(segmentEndIdx)]), yLimits(2)*0.9, ...
            ['Region: ' num2str(segmentStartIdx) ' to ' num2str(segmentEndIdx)], ...
            'Color', pinColors(maxPin,:), 'HorizontalAlignment', 'center');
    end
        
    % Mark max and min points of selected pin
    xline(timePoints(maxIndex), '--', 'LineWidth', 1.5, 'Color', pinColors(maxPin,:));
    xline(timePoints(minIndex), ':', 'LineWidth', 1.5, 'Color', pinColors(maxPin,:));
    title(['Force Data - ' names{i}])
    xlabel('Time Point')
    ylabel('Force')
    legend('Pin 1', 'Pin 2', 'Pin 3', 'Pin 4', 'Pin 5', 'Pin 6', 'Pin 7', 'Pin 8', 'Location', 'best')

    % subplot(2,1,2)
    % hold on
    % for j = 1:8
    %     bar(j, mean(averageForces(i, j))-allData{i}.(['Force' num2str(j)])(1), 'FaceColor', pinColors(j,:));
    % end
    % % Set x-ticks and labels
    % xticks(1:8);
    % xticklabels({'Pin 1', 'Pin 2', 'Pin 3', 'Pin 4', 'Pin 5', 'Pin 6', 'Pin 7', 'Pin 8'});
    % % Set title and labels
    % title('Average Forces for Each Pin');
    % xlabel('Pin Number');
    % ylabel('Average Force (N)');
    % 
    % hold off
    
end


disp('Done')