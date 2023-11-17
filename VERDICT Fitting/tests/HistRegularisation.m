% Matlab script to test implementation of histology regularisation


%% Define VERDICT scheme

V0 = [1,2,0];
% V1 = [3.9, 23.8, 90];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];
% V5 = [18.9, 38.8, 3000];

% Addition of b=0
Vs = [...
%     V0; V1;...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    V0; V5
    ];


%% Define tissue parameters

% Volume fractions
fIC = 0.4;
fVASC = 0.1;
fEES = 1-fIC-fVASC;

% Spheres distribution
tissueRmin = 0.1;
tissueRmax = 15;
tissuenR = 100; % Number of radii

Rs = linspace(tissueRmin, tissueRmax, tissuenR);

muR = 11;
sigmaR = 1;

fRs = normpdf(Rs, muR, sigmaR);

% Noise level
NoiseSigma = 0.01 ;


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

% Rs used in fitting
fitRs = linspace(0.1, 15.1, 25);

% Number of comparments
ncompart = 1;

% Noisy Matrix bool
NoisyMatrix = false;


[fIC_fit, fEES_fit, fVASC_fit, R, rmse, A, t, opt, x] = verdict_fit( ...
    Vscheme, ...
    Y, ...
    ncompart = ncompart, ...
    Rs = fitRs, ...
    solver = 'histreg',...
    NoisyMatrix = NoisyMatrix,...
    NoiseSigma = NoiseSigma);

fRs_fit = x(1:length(fitRs));

figure;
plot(Rs, fRs, '-*')
hold on
plot(fitRs, fRs_fit, '-*')

