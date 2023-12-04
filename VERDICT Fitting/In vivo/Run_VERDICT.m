% MATLAB script to run VERDICT processing on specified patients

% Read list of patents with downloaded data
PatTable = readtable("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\INNOVATE Data\DownloadedPatNums.xlsx", ReadVariableNames=false);

PatNums = string( PatTable{:,1} );


% Define output folder
OutputFolder = string(fileread("..\..\output_folder.txt"));


%%% ======== Model types

%% RDI Model

% Model type
ModelType = "RDI";
solver = 'lsqnonnegTikhonov';

% Run VERDICT processing
for indx = 1:size(PatNums,1)
    PatNum = PatNums(indx);
    disp(["----------->" PatNums(indx)])
    VERDICT(PatNum, ModelType, solver = solver, parent_folder=OutputFolder);
end

pause(60)


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

% % Model 2: No VASC model
% 
% % Model type
% ModelType = "Model 2";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 3: No VASC Reduced Rs [0.1,5.1,10.1,15.1]
% 
% 
% % Model type
% ModelType = "Model 3";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 4: No VASC Reduced Rs [3,6,9,12]
% 
% % Model type
% ModelType = "Model 4";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 5: No VASC Reduced Rs [6,9,12,15]
% 
% % Model type
% ModelType = "Model 5";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 6: No VASC Reduced Rs [0.1,3.1,6.1,9.1]
% 
% % Model type
% ModelType = "Model 6";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 7: No VASC, No 3000, Reduced Rs [0.1,7.6,15.1]
% 
% % % Model type
% % ModelType = "Model 7";
% % solver = 'lsqnonneg';
% % 
% % % Run VERDICT processing
% % for indx = 1:size(PatNums,1)
% %     PatNum = PatNums(indx);
% %     disp(["----------->" PatNums(indx)])
% %     VERDICT(PatNum, ModelType, solver = solver);
% % end
% % 
% % pause(60)

% %% Model 8: No VASC, No 3000, Reduced Rs [3.75, 7.5, 11.25]
% 
% % Model type
% ModelType = "Model 8";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% %% Model 9: No VASC, No 3000, Reduced Rs [7.5, 11.25, 15]
% 
% % Model type
% ModelType = "Model 9";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% %% Model 10: No VASC, No 3000, Reduced Rs [0.1, 3.85, 7.6]
% 
% % Model type
% ModelType = "Model 10";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% %% Model 11:  No 3000
% 
% % Model type
% ModelType = "Model 11";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% %% Model 12: No 3000, No 2000
% % Model type
% ModelType = "Model 12";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 13: No Vasc, SVD optimal Rs 
% 
% % Model type
% ModelType = "Model 13";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 14: No VASC, No 3000, SVD optimal Rs
% 
% % Model type
% ModelType = "Model 14";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 

% %% Model 15: No VASC, No 3000
% 
% % Model type
% ModelType = "Model 15";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% %% Model 16: No 3000, Rs = [0.1, 7.6, 15.1]
% 
% % Model type
% ModelType = "Model 16";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 17: No 2000, No 3000, Rs = [3.75, 11.25]
% 
% % Model type
% ModelType = "Model 17";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 18: Rs = [0.1, 5.1, 10.1, 15.1]
% 
% % Model type
% ModelType = "Model 18";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

% %% Model 19: No VASC
% 
% % Model type
% ModelType = "Model 19";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)
% 
% %% Model 20: No VASC
% 
% % Model type
% ModelType = "Model 20";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)
% 
% %% Model 21: No VASC
% 
% % Model type
% ModelType = "Model 21";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for indx = 1:size(PatNums,1)
%     PatNum = PatNums(indx);
%     disp(["----------->" PatNums(indx)])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% pause(60)

