# TriClip Experiment Analysis

## Overview

This repository contains the data analysis pipeline for the TriClip experiment conducted at the Rausch Lab. The experiment investigates the biomechanical effects of different TriClip configurations on tricuspid valve annulus force distributions, pressure, and regurgitant flow rate. 

## Experiment Description

The TriClip experiment examines how various clip configurations affect tricuspid valve function by measuring:
- **Force distributions** across 8 measurement points around the valve annulus
- **Pressure** in the ventricle
- **Flow rates** through the valve

### Test Configurations

we tested the following states:
- **Healthy**: Baseline normal valve function
- **Diseased**: Pathological valve state (reference condition)
- **Intervention configurations**:
  - **AS**: Anterior-septal clip placement
  - **AP**: Anterior-posterior clip placement  
  - **SP**: Septal-posterior clip placement
  - **ASAP**: Combined anterior-septal and anterior-posterior clips
  - **SPAS**: Combined septal-posterior and anterior-septal clips
  - **SPAP**: Combined septal-posterior and anterior-posterior clips

## Repository Structure

```
triclip-experiment-repository/
├── code/                           # Analysis scripts
│   ├── triclip_Analysis.m         # Main MATLAB analysis pipeline
│   ├── pubPlot.m                  # Publication-quality plotting function
│   └── triclip_satistics.R        # R statistical analysis scripts
├── data/                          # Processed data files
├── statistics/                    # Statistical analysis outputs
└── README.md                      # This file
```


## Data Analysis Pipeline

### 1. Data Loading and Preprocessing
- Loads CSV files containing time-series measurements
- Applies calibration factors to force measurements
- Zeros initial values for relative measurements
- Identifies stable measurement regions for analysis

### 2. Force Calibration
The system uses 8 strain gauge force sensors with individual calibration factors:
```matlab
force_calibration = [
    -1621.4975*3.3  % Pin 1
    -1291.6308*3.3  % Pin 2
    -1368.1537*3.3  % Pin 3
    -1017.2041*3.3  % Pin 4
    -1106.5112*3.3  % Pin 5
    -1089.2077*3.3  % Pin 6
    -2948.5710*3.3  % Pin 7
    -1427.9047*3.3  % Pin 8
];
```

### 3. Flow Rate Calibration
Flow rate measurements are calibrated using:
- Calibration factor: 200.44 ml/s/V
- Zero-offset correction applied

### 4. Analysis Outputs

#### Visualizations Generated:
1. **Flow Rate Bar Graphs**: Average flow rates by intervention type with error bars
2. **Pressure Bar Graphs**: Average pressures by intervention type with error bars
3. **Force Difference Bar Graphs**: Force changes relative to diseased state (6 subplots)
4. **Force Contour Visualizations**: 2D contour plots showing spatial force distributions with:
   - Baseline contour (diseased state reference)
   - Force-displaced contour (intervention effects)
   - Quiver arrows indicating force magnitude and direction

#### Statistical Analysis:
- Mean and standard deviation calculations across tests
- Force difference calculations (intervention - diseased)
- Per-test statistical summaries
- Cross-intervention comparisons

## Usage

### Prerequisites
- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox
- R (for statistical analysis scripts)

### Running the Analysis

1. **Main Analysis Pipeline**:
   ```matlab
   % Navigate to the code directory
   cd('path/to/triclip-experiment-repository/code')
   
   % Run the main analysis
   triclip_Analysis
   ```

2. **Key Analysis Sections**:
   - Data loading and preprocessing (automatic)
   - Flow rate analysis and visualization
   - Pressure analysis and visualization  
   - Force difference analysis (bar graphs)
   - Force contour visualization (spline plots)

### Customization

#### Adding New Tests:
1. Place new CSV data in appropriately named date folders
2. Create corresponding test index files in `testidx/` directory
3. Update folder paths in `triclip_Analysis.m`

#### Modifying Visualizations:
- Adjust `force_colors` array for different intervention colors
- Modify `force_measurement_pin_positions` for different sensor layouts
- Change `arrow_length_scale_factor` for arrow visibility

#### Publication Plots:
Use `pubPlot()` function for publication-quality formatting:
```matlab
pubPlot('Width', 'double', 'SpacingOffset', 2)
```

## Data Export

The analysis pipeline can export results to Excel files:
- `TriClipXT_Statistics.csv`: Force differences relative to diseased state
- `TriClipXT_Flow_Statistics.csv`: Average measurements by intervention


## Contributing

When modifying the analysis pipeline:
1. Maintain descriptive variable names
2. Add detailed comments for complex calculations
3. Use consistent naming conventions
4. Test with multiple datasets before committing
5. Update this README for significant changes

## Contact

Rausch Lab  
[Institution/Department Information]

## License

[Add appropriate license information]