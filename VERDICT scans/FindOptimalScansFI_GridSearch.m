% fmincon not really working so going to do a grid search instead...

clear all;

%% Define model parameters

% First, define model type
modeltype = 'RDI';

% Define model parameter index (to optimise scan for)
paramIndx = 2; 

% Define number of paramaters (number of volume fractions)
Nparam = 4;


% Define min and max volume fractions 
fmin = 0;
fmax = 1;
nf = 11;

% Construct paramater grid
fs = linspace(fmin, fmax, nf);

[varargout{1:Nparam}] = ndgrid(fs);

paramGrid = cat(Nparam+1, varargout{1:Nparam});

paramSpacings = ((fmax-fmin)/nf)*ones(Nparam,1);


%% Define scan parameter ranges

% delta
deltaMin = 15;
deltaMax = 50;
ndelta = 20;
deltas = linspace(deltaMin, deltaMax, ndelta);

% Delta
DeltaMin = 25;
DeltaMax = 60;
nDelta = 20;
Deltas = linspace(DeltaMin, DeltaMax, nDelta);

% bval
bMin = 2500;
bMax = 3000;
nbval = 5;
bvals = linspace(bMin, bMax, nbval);

% Gradient strength
Gmax = 60;


% Define grid
[deltaVals, DeltaVals, bVals] = ndgrid(deltas, Deltas, bvals);
ScanParamGrid = cat(4, deltaVals, DeltaVals, bVals);
FlatScanParamGrid = reshape(ScanParamGrid, ndelta*nDelta*nbval, 3);


%% Iterate over grid

FIs = zeros(length(FlatScanParamGrid), 1);

for scanparamIndx = 1:length(FlatScanParamGrid)

    disp('scan param indx')
    disp(scanparamIndx)

    scan_params = FlatScanParamGrid(scanparamIndx,:);

    delta = scan_params(1);
    Delta = scan_params(2);
    bval = scan_params(3);
% 
%     disp(delta)
%     disp(Delta)
%     disp(bval)

    % Test delta<Delta
    if Delta <= delta
        disp('failed delta')
        continue

    % Test gradient strength
    elseif stejskal(delta, Delta, bval = bval) > Gmax
        disp('failed G')
        continue
    else
        disp('ok')
    end


    % == Evaluate Fisher Information
    FI = FisherInfoGaussian(paramIndx, paramGrid, paramSpacings, scan_params, modeltype = modeltype);
    FIs(scanparamIndx,1) = FI;







end




