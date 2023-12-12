% Matlab script to compare fitting ability of different schemes

%% Define schemes to evaluate

% == Original scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
%     V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
%     V0; V5
    ];

original_scheme = BuildScheme(Vs);


% == New scheme

V0 = [1,2,0];
V1 = [18.3, 31.0, 1800];
V2 = [18.4, 33.2, 1980];
V3 = [27.8, 35.6, 1942];

% Addition of b=0
Vs = [...
    V0; V1;...
    V0; V2;...
    V0; V3;...
    ];

new_scheme = BuildScheme(Vs);



%% Set up tissue parameters

% == Define possible ranges of paramaters

fICs = [0.1, 0.9];
fVASCs = [0, 0.1];

DEESs = [1.5,2.5];

Rs = linspace(0.1, 15.1, 50);

muRs = [6,9];
sigmaRs = [1,3];


% == Specify tissue params

fIC = fICs(1) + ( fICs(2) - fICs(1) )*rand();
simfICs(iterIndx) = fIC;

fVASC = fVASCs(1) + ( fVASCs(2) - fVASCs(1) )*rand();
simfVASCs(iterIndx) = fVASC;

fEES = 1-fIC-fVASC;

muR = muRs(1) + ( muRs(2) - muRs(1) )*rand();
simmuRs(iterIndx) = muR;

sigmaR = sigmaRs(1) + ( sigmaRs(2) - sigmaRs(1) )*rand();

fRs = normpdf(Rs, muR, sigmaR);

DEES = DEESs(1) + ( DEESs(2) - DEESs(1) )*rand();
simDEESs(iterIndx) = DEES;


%% Evaluate fitting performance of different schemes

% Noise level
NoiseSigma = 0.05;
NoisyMatrix = false;

% Repetitions
Nrep = 2000;

% == Fitting

% Model name
ModelName = 'RDI';


% original scheme
[original_bias, original_variance] = EvaluatefICfitting( ...
    original_scheme, ...
    ModelName, ...
    Nrep = Nrep,...
    fIC = fIC, ...
    fEES = fEES, ...
    DEES = DEES, ...
    NoiseSigma = NoiseSigma, ...
    fRs=fRs,...
    Rs = Rs,...
    NoisyMatrix=NoisyMatrix);

% new scheme
[new_bias, new_variance] = EvaluatefICfitting( ...
    new_scheme, ...
    ModelName, ...
    Nrep = Nrep,...
    fIC = fIC, ...
    fEES = fEES, ...
    DEES = DEES, ...
    NoiseSigma = NoiseSigma, ...
    fRs=fRs,...
    Rs = Rs,...
    NoisyMatrix=NoisyMatrix);

