% MATLAB test script to test adapted VERDICT fitting on in vivo data

% Specify INNOVATE patient number 
PatNum = "BAR_003";

% Specify model type (noVASCex3000)
ModelNum = "Model 999";

% Define output folder
OutputFolder = string(fileread("..\..\output_folder.txt"));

% Call VERDICT function
VERDICT(PatNum, ModelNum, solver = 'histreg', parent_folder = OutputFolder)