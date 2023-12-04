% MATLAB script to test RDI_fit function

clear;

%% Define VERDICT scheme

V0 = [1,2,0];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];

% Addition of b=0
Vs = [...
    V0; V2;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    ];


%% Define tissue parameters

% Volume fractions
fIC = 0.4;
fVASC = 0.0;
fEES = 1-fIC-fVASC;

% Spheres distribution
tissueRmin = 0.1;
tissueRmax = 15.1;
tissuenR = 100; % Number of radii

% Normal distributon
Rs = linspace(tissueRmin, tissueRmax, tissuenR);

muR = 10;
sigmaR = 2;

fRs = normpdf(Rs, muR, sigmaR);


% Noise level
NoiseSigma = 0.01;


%% Simulate signal

% simulate signal
[Signals, Vscheme] = SimulateSignal(Vs, fIC, fVASC, Rs, fRs);

% Add noise
SignalsNoisy = abs( AddNoise(Signals, NoiseSigma = NoiseSigma) );
SignalsNoisy(1:2:end) = 1;


%% Apply fitting

% Construct Y matrix
Y = zeros([1,1,length(SignalsNoisy)]);
Y(1,1,:) = SignalsNoisy;

% Noisy Matrix bool
NoisyMatrix = false;

% RDI fit
[fIC_fit, fEES_fit, rmse] = RDI_fit( ...
    Vscheme, ...
    Y,...
    NoisyMatrix=NoisyMatrix,...
    NoiseSigma = NoiseSigma...
    );
