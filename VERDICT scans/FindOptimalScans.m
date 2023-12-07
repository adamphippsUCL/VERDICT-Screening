% Matlab script to find scan parameters which minimise condition number

% Run function minimiseCN multiple times and find case with lowest CN


% Number of scans
nscan = 3;

% Fix b values?
bvals = [];



%% Iterate optimisation function

% Number of optimisation iterations
Niter = 10;

% Loop
scan_params = [];
CNs = [];

for iterIndx = 1:Niter

    [params, CN] = minimiseCN(nscan, bvals = bvals);

    scan_params = [scan_params, params];
    CNs = [CNs, CN];

end


%% Find best solution

[minCN, indx] = min(CNs);

best_params = scan_params(:,indx);
