% MATLAB script to conduct comparison of different models on simulated data

% 1. Randomise tissue paramaters
% 2. Evaluate fitting performance
% 3. repeat and find distribution of bias variances
% 
% Additional (can I show how bias variance varies with each tissue paramater?)


%% Define VERDICT scheme

V0 = [1,2,0];
V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
V5 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
    V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    V0; V5
    ];

scheme = BuildScheme(Vs);


%% Define noise level

NoiseSigma = 0.05;
NoisyMatrix = false;

%% Define simulation settings

% Number of parameter sets
Niter = 400;

% Number of repetitions at each parameter set
Nrep = 200;

%% Define possible ranges of paramaters

fICs = [0.1, 0.9];
fVASCs = [0, 0.1];

DEESs = [1.5,2.5];

Rs = linspace(0.1, 15.1, 50);

muRs = [4,12];
sigmaRs = [1,3];


%% Set up different models

% == Original VERDICT (all scans)

OrigAll_bias = zeros(Niter, 1);
OrigAll_var = zeros(Niter, 1);

% == Original VERDICT (ex. 90, 3000)

OrigEx903000_bias = zeros(Niter, 1);
OrigEx903000_var = zeros(Niter, 1);

% == no VASC (all scans)

noVASCAll_bias = zeros(Niter, 1);
noVASCAll_var = zeros(Niter,1);

% == no VASC (ex. 90, 3000)

noVASCEx903000_bias = zeros(Niter, 1);
noVASCEx903000_var = zeros(Niter, 1);

% == RDI (all scans)

RDIAll_bias = zeros(Niter, 1);
RDIAll_var = zeros(Niter,1);

% == RDI (ex. 90, 3000)

RDIEx903000_bias = zeros(Niter, 1);
RDIEx903000_var = zeros(Niter, 1);

% == RDI (two scans: 1500, 2000)

RDITwoScan_bias = zeros(Niter, 1);
RDITwoScan_var = zeros(Niter, 1);


%% Iterate over many tissue paramaters


simfICs = zeros(Niter, 1);
simfVASCs = zeros(Niter, 1);
simmuRs = zeros(Niter, 1);
simDEESs = zeros(Niter, 1);

for iterIndx = 1:Niter

    % Specify tissue params

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

    % == Original VERDICT (all scans)

    [bias, variance] = EvaluatefICfitting( ...
        scheme, ...
        'Original VERDICT', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    OrigAll_bias(iterIndx) = bias;
    OrigAll_var(iterIndx) = variance;



    % == no VASC VERDICT (all scans)

    [bias, variance] = EvaluatefICfitting( ...
        scheme, ...
        'no VASC VERDICT', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    noVASCAll_bias(iterIndx) = bias;
    noVASCAll_var(iterIndx) = variance;



    % == RDI (all scans)

    [bias, variance] = EvaluatefICfitting( ...
        scheme, ...
        'RDI', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    RDIAll_bias(iterIndx) = bias;
    RDIAll_var(iterIndx) = variance;


    % == Original VERDICT (ex. 90, 3000)

    [bias, variance] = EvaluatefICfitting( ...
        scheme(3:end-2), ...
        'Original VERDICT', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    OrigEx903000_bias(iterIndx) = bias;
    OrigEx903000_var(iterIndx) = variance;

 
     % == no VASC VERDICT (ex. 90, 3000)

    [bias, variance] = EvaluatefICfitting( ...
        scheme(3:end-2), ...
        'no VASC VERDICT', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    noVASCEx903000_bias(iterIndx) = bias;
    noVASCEx903000_var(iterIndx) = variance;



    % == RDI (ex. 90, 3000)

    [bias, variance] = EvaluatefICfitting( ...
        scheme(3:end-2), ...
        'RDI', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    RDIEx903000_bias(iterIndx) = bias;
    RDIEx903000_var(iterIndx) = variance;



     % == RDI (Two scan: 1500, 2000)

    [bias, variance] = EvaluatefICfitting( ...
        scheme([false false false false true true true true false false]), ...
        'RDI', ...
        Nrep = Nrep,...
        fIC = fIC, ...
        fEES = fEES, ...
        DEES = DEES, ...
        NoiseSigma = NoiseSigma, ...
        fRs=fRs,...
        Rs = Rs,...
        NoisyMatrix=NoisyMatrix);

    RDITwoScan_bias(iterIndx) = bias;
    RDITwoScan_var(iterIndx) = variance;


end

% 
% figure;
% scatter(simmuRs, OrigAll_bias)
% ylim([-0.3,0.3])
% xlabel('muR')
% ylabel('bias')
% title('Original All')
% 
% figure;
% scatter(simmuRs, RDIEx903000_bias)
% ylim([-0.3,0.3])
% xlabel('muR')
% ylabel('bias')
% title('RDI ex. 90, 3000')


% figure;
% scatter(simfICs, OrigAll_bias)
% ylim([-0.3,0.3])
% xlabel('fIC')
% ylabel('bias')
% title('Original All')
% 
% figure;
% scatter(simfIC, RDIEx903000_bias)
% ylim([-0.3,0.3])
% xlabel('fIC')
% ylabel('bias')
% title('RDI ex. 90, 3000')


figure;
scatter(simDEESs, OrigAll_bias)
ylim([-0.3,0.3])
xlabel('DEES')
ylabel('bias')
title('Original All')

figure;
scatter(simDEESs, RDIEx903000_bias)
ylim([-0.3,0.3])
xlabel('DEES')
ylabel('bias')
title('RDI ex. 90, 3000')


figure;
scatter(ones(Niter,1), OrigAll_bias)
hold on
scatter(2*ones(Niter, 1), noVASCAll_bias)
hold on
scatter(3*ones(Niter, 1), RDIAll_bias)
hold on
scatter(4*ones(Niter,1), OrigEx903000_bias)
hold on
scatter(5*ones(Niter, 1), noVASCEx903000_bias)
hold on
scatter(6*ones(Niter, 1), RDIEx903000_bias)
hold on
scatter(7*ones(Niter, 1), RDITwoScan_bias)
% Overlay boxplots
boxplot([OrigAll_bias, noVASCAll_bias, RDIAll_bias, OrigEx903000_bias, noVASCEx903000_bias, RDIEx903000_bias, RDITwoScan_bias], [1,2, 3, 4, 5, 6, 7])
% Add x ticks
xticks([1 2 3 4 5 6 7])
xticklabels({'Orig All', 'no VASC All', 'RDI All', 'Orig ex. 90, 3000', 'no VASC ex. 90, 3000', 'RDI ex. 90, 3000', 'RDI ex. 90, 500, 3000'})
ylabel('fIC bias')
grid on

figure;
scatter(ones(Niter,1), OrigAll_var)
hold on
scatter(2*ones(Niter, 1), noVASCAll_var)
hold on
scatter(3*ones(Niter, 1), RDIAll_var)
hold on
scatter(4*ones(Niter,1), OrigEx903000_var)
hold on
scatter(5*ones(Niter, 1), noVASCEx903000_var)
hold on
scatter(6*ones(Niter, 1), RDIEx903000_var)
hold on
scatter(7*ones(Niter, 1), RDITwoScan_var)
% Overlay boxplots
boxplot([OrigAll_var, noVASCAll_var, RDIAll_var, OrigEx903000_var, noVASCEx903000_var, RDIEx903000_var, RDITwoScan_var], [1, 2, 3, 4, 5, 6, 7])
% Add x ticks
xticks([1 2 3 4 5 6 7])
xticklabels({'Orig All', 'no VASC All', 'RDI All', 'Orig ex. 90, 3000', 'no VASC ex. 90, 3000', 'RDI ex. 90, 3000', 'RDI ex. 90, 500, 3000'})
ylabel('fIC variance')
grid on


