% Matlab function to evaluate fIC performance of specified model
function [bias, variance] = EvaluatefICfitting(scheme, ModelName, tissue_params, opts)

arguments

    scheme % Define scan parameters
    ModelName % Name of model used in fitting

    tissue_params.fIC % Tissue parameters for signal simulation
    tissue_params.DEES = 2
    tissue_params.Rs = linspace(0.1,15.1,50)
    tissue_params.fRs = []
    tissue_params.fEES

    opts.Nrep = 1000
    opts.NoiseSigma = 0.05 % Options for simulation/fitting
    opts.NoisyMatrix = false
    opts.solver
end


%% First simulate signals over scheme

fVASC = 1 - tissue_params.fIC - tissue_params.fEES;

signals = SimulateProstateSignal(scheme, ...
    fIC = tissue_params.fIC, ...
    Rs = tissue_params.Rs, ...
    fRs = tissue_params.fRs, ...
    fEES = tissue_params.fEES,...
    DEES = tissue_params.DEES, ...
    fVASC = fVASC);



%% Specify model fitting function

switch ModelName

    case 'Original VERDICT'

        % specify fitting function and options
        fitfunc = @verdict_fit;
        options.ncompart = 2;

    case 'no VASC VERDICT'

        % specify fitting function and options
        fitfunc = @verdict_fit;
        options.ncompart = 1;


    case 'RDI'

        % specify fitting function and options
        fitfunc = @RDI_fit;
        options.ncompart = 1;

end

% Further fitting specifications
options.NoisyMatrix = opts.NoisyMatrix;
options.NoiseSigma = opts.NoiseSigma;



%% Iterate over noise instances

% Empty array for fitted fIC values
fIC_fits = zeros(opts.Nrep,1);

for noiseIndx = 1:opts.Nrep

    noiseIndx

    % Add noise to signals
    SignalsNoisy = abs( AddNoise(signals, NoiseSigma = opts.NoiseSigma) );
    % Correct b=0 signals
    SignalsNoisy([scheme.bval] == 0) = 1;
    
    Y = zeros([1,1,length(SignalsNoisy)]);
    Y(1,1,:) = SignalsNoisy;

    % Apply fitting
    [fIC_fit, fEES_fit] = fitfunc( ...
        scheme, ...
        Y, ...
        ncompart = options.ncompart, ...
        NoisyMatrix = options.NoisyMatrix, ...
        NoiseSigma = options.NoiseSigma ...
        );

    % Append to array
    fIC_fits(noiseIndx) = fIC_fit;

end


%% Calculate bias and variance

bias = mean(fIC_fits - tissue_params.fIC);
variance = var(fIC_fits - tissue_params.fIC);

end


