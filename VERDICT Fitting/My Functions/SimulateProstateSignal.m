function [signals] = SimulateProstateSignal(scheme, tissue_params)

arguments
    scheme % scheme of scan paramaters

    tissue_params.fIC % Tissue parameters for signal simulation
    tissue_params.DIC = 2
    tissue_params.Rs = linspace(0.1,15.1,50)
    tissue_params.fRs = []
    tissue_params.fEES
    tissue_params.DEES = 2
    tissue_params.fVASC = 0
    tissue_params.DVASC = 8


end

% Default radii distribution
if isempty(tissue_params.fRs)
    tissue_params.fRs = normpdf(tissue_params.Rs, 8, 2);
end

% Normalise and scale fRs
tissue_params.fRs = tissue_params.fIC*tissue_params.fRs/sum(tissue_params.fRs);


% === Bits and bobs

nscheme = length(scheme);
nR = length(tissue_params.Rs);


% === Simulating signals

% Empty arrays for signals
sIC = zeros([nscheme nR]) ;
sEES = zeros([nscheme 1]) ;
sVASC = zeros([nscheme 1]) ;
stot = zeros([nscheme 1]) ;


for ischeme = 1:nscheme

    % IC signal for different radii
    for ir = 1:nR
        sIC(ischeme,ir) = tissue_params.fRs(ir)*sphereGPD(scheme(ischeme).delta, scheme(ischeme).DELTA, ...
            scheme(ischeme).G, tissue_params.Rs(ir), tissue_params.DIC);
    end

    % EES signal
    sEES(ischeme) = ball(scheme(ischeme).bval, tissue_params.DEES) ;
    
    % VASC signal
    sVASC(ischeme) = astrosticks(scheme(ischeme).bval, tissue_params.DVASC) ;

    % Total signal
    stot(ischeme) = sum(sIC(ischeme,:)) + tissue_params.fEES*sEES(ischeme) + ...
        tissue_params.fVASC*sVASC(ischeme);

    
end

signals = stot;


end