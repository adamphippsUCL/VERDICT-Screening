% MATLAB script to run VERDICT processing on specified patients

% Read list of patents with downloaded data
PatTable = readtable("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE Data\DownloadedPatNums.xlsx", ReadVariableNames=false);

PatNums = string( PatTable{:,1} );


% Define output folder
OutputFolder = string(fileread("..\..\output_folder.txt"));


%%% ======== Model types

% %% RDI Model
% 
% % Model type
% ModelType = "RDIex903000";
% solver = 'lsqnonnegTikhonov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver, parent_folder=OutputFolder);
% end
% 
% pause(60)
% 
% 
% % Model 1: Original model 
% 
% % Model type
% ModelType = "Original VERDICT";
% solver = 'lsqnonnegTikhonov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver, parent_folder=OutputFolder);
% end
% 
% pause(60)


% Model: No VASC ex 90, 3000

% Model type
ModelType = "noVASCex903000";
solver = 'lsqnonnegTikhonov';

% Run VERDICT processing
for indx = 1:size(PatNums,1)
    PatNum = PatNums(indx);
    disp(["----------->" PatNums(indx)])
    VERDICT(PatNum, ModelType, solver = solver);
end

pause(60)
