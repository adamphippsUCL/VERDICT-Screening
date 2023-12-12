% MATLAB script to find scan parameters which maximise Fisher information 
% about a specified model parameter.
clear all;

% First, define model type
modeltype = 'RDI';

% Define model parameter index
paramIndx = 1; 

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


%% Optimisation

% == Define function to optimise
fun = @FI;

% == Specify starting parameters
delta0 = 20 + (30-20)*rand();
Delta0 = 30 + (40-30)*rand();
bval0 = 500 + (2000-500)*rand();

x0 = [delta0, Delta0, bval0];

% == Specify upper and lower paramater bounds

lb = [0, 0, 0];
ub = [40, 50, 3000];

% == Equality constraint (None)

Aeq = zeros(3,3);
beq = zeros(3,1);

% == Inequality constraint (delta < Delta)

A = [1, -1, 0; 0, 0, 0; 0, 0, 0];
b =  zeros(3,1);

% == Non-linear gradient constraint

nonlcon = @gradnonlcon;


% == Optimisation options
opts = optimoptions('fmincon');
opts.DiffMinChange = 1;
opts.Display = 'iter';
% opts.ScaleProblem = true;



[params, FishInf] = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, opts);

FishInfinit = FI(x0);


function f = FI(x)

    paramIndx = evalin('base', 'paramIndx');
    paramGrid = evalin('base', 'paramGrid');
    paramSpacings = evalin('base', 'paramSpacings');   
    modeltype = evalin('base', 'modeltype');      

    % Evaluate Fisher information at specified scan parameters
    f = - FisherInfoGaussian(paramIndx, paramGrid, paramSpacings, x, modeltype = modeltype);
    disp(x)
    disp(f)


end



