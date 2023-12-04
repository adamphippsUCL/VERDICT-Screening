% MATLAB code to evaluate fIC fitting performance of RDIover multiple noise instances


clear;

%% Define VERDICT scheme

V0 = [1,2,0];
V2 = [11.4, 31.3, 500];
V3 = [23.9, 43.8, 1500];
V4 = [12.4, 32.3, 2000];

% Addition of b=0
Vs = [...
    V0; V2;...
    V0; V3;...
    V0; V4;...
    ];


%% Define tissue parameters

% Volume fractions
fIC = 0.6;
fVASC = 0.0;
fEES = 1-fIC-fVASC;
% Spheres distribution
tissueRmin = 0.1;
tissueRmax = 15.1;
tissuenR = 100; % Number of radii

% Normal distributon
Rs = linspace(tissueRmin, tissueRmax, tissuenR);

muR = 8;
sigmaR = 10;

fRs = normpdf(Rs, muR, sigmaR);
% 
% % Random uniform distribution
% Rmin = 5;
% Rmax = 10;
% nR = 100;
% Rs = Rmin + (Rmax-Rmin)*rand(1,100);
% fRs = (1/nR)*ones(1,nR);
% 
% muR = sum(fRs.*Rs);

%% Define noise

% Noise level (b=0 SNR)
NoiseSigma = 0.05;


% Number of instances
Ninst = 1000;


%% Define fitting settings

% Noisy matrix
NoisyMatrix = false;

% muR and sigmaR for RDI fitting
RDI_muR = 8;
RDI_sigmaR = 4;

RDI_fRs = normpdf(Rs, RDI_muR, RDI_sigmaR);


%% Iterate over noise instances

% Results
RDI_fICs = zeros(Ninst,1);
noVASC_fICs = zeros(Ninst,1);


for indx = 1:Ninst

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


    % No VASC fitting
    [noVASC_fIC, noVASC_fEES, noVASC_fVASC_fit, noVASC_R, noVASC_rmse] = verdict_fit( ...
        ...
        Vscheme, ...
        Y, ...
        ncompart = 1, ...
        solver = 'lsqnonnegTikhonov',...
        NoisyMatrix = NoisyMatrix,...
        NoiseSigma = NoiseSigma ...
        );
    
    % Record key results
    noVASC_fICs(indx) = noVASC_fIC;

    % RDI fitting
    [RDI_fIC, RDI_fEES, RDI_rmse] = RDI_fit( ...
        Vscheme, ...
        Y,...
        muR = RDI_muR, ...
        sigmaR = RDI_sigmaR ...
        );

    % Record key results
    RDI_fICs(indx) = RDI_fIC;

end


%% Evaluation of fitting performance

% == fIC fitting

noVASC_fIC_errors = noVASC_fICs - fIC;
RDI_fIC_errors = RDI_fICs - fIC;

fig1 = figure;
% Induvidual noVASC points
scatter(ones(size(noVASC_fIC_errors)), ...
    noVASC_fIC_errors, Marker ='*', ...
    MarkerEdgeAlpha = 0.1, ...
    MarkerFaceAlpha = 0.1,...
    DisplayName = 'noVASC')
hold on
% Induvidual RDI points
scatter(2*ones(size(RDI_fIC_errors)), ...
    RDI_fIC_errors, Marker ='*', ...
    MarkerEdgeAlpha = 0.1, ...
    MarkerFaceAlpha = 0.1, ...
    DisplayName = 'RDI')
hold on
% Overlay boxplots 
boxplot([noVASC_fIC_errors, RDI_fIC_errors],[1,2])
legend
title('fIC fitting error')
ylim([-0.6,0.6])
grid on

saveas(fig1, 'Figures/fIC_errors.png')


% == Plot simulated and fitting radii distributions


fig2 = figure;
plot(Rs, fRs, 'p', DisplayName = ['Simulated (mu = ' num2str(muR) ', sigma = ' num2str(sigmaR) ')'])
hold on
plot(Rs, RDI_fRs, '*', DisplayName = ['RDI fitting (mu = ' num2str(RDI_muR) ', sigma = ' num2str(RDI_sigmaR) ')'])
legend;
title(['fIC ' num2str(fIC)])
saveas(fig2, 'Figures/SimulatedFittedfRs.png')