
function VERDICT(pat_num, model_type, opts)

% MATLAB function to carry out VERDICT processing 

% INNOVATE patient number and VERDICT model type specified as inputs

% Model type can be: 'Original' or 'No VASC'

arguments
    pat_num % INNOVATE patient number
    model_type % Type of VERDICT model to fit
    opts.INNOVATE_path = "D:\UCL PhD Imaging Data\INNOVATE VERDICT\" % Path to INNOVATE imaging folder
    opts.parent_folder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Submission\Outputs" % Parent folder to save outputs
    opts.solver = 'lsqnonnegTikhonov'
end

% Define DICOM folder path
DICOM_path = join( [opts.INNOVATE_path pat_num "\scans"], "");
if exist(DICOM_path, "dir")
    disp('')
else
    DICOM_path = join( [opts.INNOVATE_path pat_num], "");
end

% Define output folder
output_path = join([opts.parent_folder "/VERDICT outputs/" pat_num "/" model_type], "");

if ~exist(output_path, "dir")
   mkdir(output_path)
end


%% Define options for model type

% Model 0: ADC excluding b2000 and b3000
if strcmp(model_type, 'Model 0')

    ncompart = 1;
    fitting_excludebvals = [2000,3000];
    Rs = [10];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver,...
        calcADC = true,...
        ADCoverride = true...
        );

% Model 0.5: ADC excluding b3000
elseif strcmp(model_type, 'Model 0.5')

    ncompart = 1;
    fitting_excludebvals = [3000];
    Rs = [10];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver,...
        calcADC = true,...
        ADCoverride = true...
        );
    

% Model 1: Original VERDICT
elseif strcmp(model_type, 'Model 1')
    
    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        solver = opts.solver... 
        );

% Model 2: No VASC, original Rs
elseif strcmp(model_type, 'Model 2')

    ncompart = 1;
    fitting_excludebvals = [90];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        solver = opts.solver...    
        );


% Model 3: No VASC, Reduced Rs = [0.1, 5.1, 10.1, 15.1]  
elseif strcmp(model_type, 'Model 3')

    ncompart = 1;
    fitting_excludebvals = [90];
    Rs = linspace(0.1,15.1,4);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 4: No VASC, Reduced Rs = [0.1, 6.1, 12.1, 18.1]  (SVD optimal Rs )
elseif strcmp(model_type, 'Model 4')

    ncompart = 1;
    fitting_excludebvals = [90];
    Rs = linspace(0.1,18.1,4);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 5: No VASC, Reduced Rs = [6,9,12,15]  
elseif strcmp(model_type, 'Model 5')

    ncompart = 1;
    fitting_excludebvals = [90];
    Rs = linspace(6,15,4);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 6: No VASC, Reduced Rs = [0.1, 3.1, 6.1, 9.1]
elseif strcmp(model_type, 'Model 6')

    ncompart = 1;
    fitting_excludebvals = [90];
    Rs = linspace(0.1,9.1,4);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...
        );

% Model 7: No VASC, No 3000, Reduced Rs = [0.1, 7.6, 15.1]  
elseif strcmp(model_type, 'Model 7')

    ncompart = 1;
    fitting_excludebvals = [90,3000];
    Rs = linspace(0.1,15.1,3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 8: No VASC, No 3000, Reduced Rs = [3.75, 7.5, 11.25]  
elseif strcmp(model_type, 'Model 8')

    ncompart = 1;
    fitting_excludebvals = [90,3000];
    Rs = linspace(3.75, 11.25, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );



% Model 9: No VASC, No 3000, Reduced Rs = [7.5, 11.25, 15]  
elseif strcmp(model_type, 'Model 9')

    ncompart = 1;
    fitting_excludebvals = [90,3000];
    Rs = linspace(7.5, 15, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 10: No VASC, No 3000, Reduced Rs = [0.1, 3.85, 7.6]  
elseif strcmp(model_type, 'Model 10')

    ncompart = 1;
    fitting_excludebvals = [90, 3000];
    Rs = linspace(0.1, 7.6, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 11: No 3000
elseif strcmp(model_type, 'Model 11')

    ncompart = 2;
    fitting_excludebvals = [3000];
    Rs = linspace(0.1, 15.1, 17);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );

% Model 12: No 3000, No 2000
elseif strcmp(model_type, 'Model 12')

    ncompart = 2;
    fitting_excludebvals = [2000, 3000];
    Rs = linspace(0.1, 15.1, 17);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 13: No Vasc, Rs = [3,6,9,12]
elseif strcmp(model_type, 'Model 13')

    ncompart = 1;
    fitting_excludebvals = [90];
    Rs = linspace(3, 12, 4);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );



% Model 14: No Vasc, No 3000, SVD optimal Rs 
elseif strcmp(model_type, 'Model 14')

    ncompart = 1;
    fitting_excludebvals = [90, 3000];
    Rs = linspace(0.1, 18.1, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 15: No VASC, No 3000
elseif strcmp(model_type, 'Model 15')

    ncompart = 1;
    fitting_excludebvals = [90, 3000];
    Rs = linspace(0.1, 15.1, 17);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 16: No 3000, Rs = [0.1, 7.6, 15.1]
elseif strcmp(model_type, 'Model 16')

    ncompart = 2;
    fitting_excludebvals = [3000];
    Rs = linspace(0.1, 15.1, 3);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );

% Model 17: No 2000, No 3000, Rs = [3.75, 11.25]
elseif strcmp(model_type, 'Model 17')

    ncompart = 2;
    fitting_excludebvals = [2000, 3000];
    Rs = linspace(3.75, 11.25, 2);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 18:  Rs = [0.1, 5.1, 10.1, 15.1]
elseif strcmp(model_type, 'Model 18')

    ncompart = 2;
%     fitting_excludebvals = [2000, 3000];
    Rs = linspace(3.75, 11.25, 2);

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        Rs = Rs,...
        solver = opts.solver...    
        );


% Model 19: No VASC 
elseif strcmp(model_type, 'Model 19')

    ncompart = 1;
%     fitting_excludebvals = [2000, 3000];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        solver = opts.solver...    
        );



% Model 20:  No VASC, No 3000
elseif strcmp(model_type, 'Model 20')

    ncompart = 1;
    fitting_excludebvals = [3000];

    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        solver = opts.solver...    
        );


% Model 21:  No VASC, No 2000, No 3000
elseif strcmp(model_type, 'Model 21')

    ncompart = 1;
    fitting_excludebvals = [2000, 3000];


    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        solver = opts.solver...    
        );


% Model 999:  RDI fitting, No VASC, No 3000
elseif strcmp(model_type, 'Model 999')

    ncompart = 1;
    fitting_excludebvals = [90, 3000];

    fitting = 'RDI';


    % Run VERDICT processing code
    [scheme, Y, fIC, fEES, fVASC, R, rmse] = verdict_Adam( ...
        convertStringsToChars(DICOM_path), ...
        convertStringsToChars(output_path), ...
        ncompart = ncompart, ...
        fitting_excludebvals = fitting_excludebvals,...
        solver = opts.solver,...    
        fitting = fitting...
        );
else
    disp('Incorrect model specification')
    exit()
end


% Save results
save([convertStringsToChars(output_path) '/scheme.mat'], 'scheme')
save([convertStringsToChars(output_path) '/Y.mat'], 'Y')   
save([convertStringsToChars(output_path) '/fIC.mat'], 'fIC')
% save([convertStringsToChars(output_path) '/fEES.mat'], 'fEES')
% save([convertStringsToChars(output_path) '/fVASC.mat'], 'fVASC')
% save([convertStringsToChars(output_path) '/R.mat'], 'R')
% save([convertStringsToChars(output_path) '/rmse.mat'], 'rmse')

end
