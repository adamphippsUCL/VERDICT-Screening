function CN = conditionNumber(ScanParamVector, params)


arguments

    ScanParamVector % vector [3*Nscan, 1] : (delta1, Delta1, b1, ...)
    % Excludes b=0 scans!!!

    % Fitting paramaters 

    % (for verdict fitting)
    params.fitting_type = 'RDI' % ('verdict for fitting to induvidual radii, RDI for distribution fitting')
    params.ncompart = 1 % number of compartments 
    params.Rmin = 0.1
    params.Rmax = 15.1
    params.nR = 17 % Number of radii in original VERDICT fitting
    params.DIC = 2;
    params.DEES = 2;
    params.DVASC = 8;

    % (for RDI fitting)
    params.muRs = [6, 7.5, 9]
    params.sigmaRs = [2,2,2]
    params.DEESs = [2]

    % Other
    params.NoisyMatrix = false;
    params.NoiseSigma = 0.05;
    params.regulariser = 'Tikhonov';
    params.lambda = sqrt(0.001);
end



%%  First construct scheme from scan vector
nscan = (1/3)*length(ScanParamVector);
scans = transpose( ( reshape( (ScanParamVector), 3, nscan) ) );

scheme = BuildScheme(scans);

% Add b=0 scans to scheme
for i = 1:nscan
    scheme(end+1).delta = 1;
    scheme(end).DELTA = 2;
    scheme(end).bval = 0;
    scheme(end).G = 0;
end

%% Second construct matrix

A = constructMatrix(scheme, ...
    fitting_type=params.fitting_type, ...
    ncompart = params.ncompart, ...
    Rmin = params.Rmin, ...
    Rmax = params.Rmax, ...
    nR = params.nR, ...
    DIC = params.DIC, ...
    DEES = params.DEES, ...
    DVASC = params.DVASC, ...
    muRs = params.muRs, ...
    sigmaRs = params.sigmaRs, ...
    DEESs = params.DEESs, ...
    NoisyMatrix = params.NoisyMatrix, ...
    NoiseSigma = params.NoiseSigma, ...
    regulariser = params.regulariser, ...
    lambda = params.lambda);



%% Find condition number
CN = cond(A);








end