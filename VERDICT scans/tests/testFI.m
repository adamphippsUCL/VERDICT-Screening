

paramIndx  = 1;

% Define parameter grid
n = 11;

fmin = 0;
fmax = 1;
df = 0.1;

fRdist1s = linspace(fmin, fmax, n);
fRdist2s = linspace(fmin, fmax, n);
fRdist3s = linspace(fmin, fmax, n);
fEESs = linspace(fmin, fmax, n);

[F1s, F2s, F3s, F4s] = ndgrid(fRdist1s, fRdist2s, fRdist3s, fEESs);

paramGrid = cat(5, F1s, F2s, F3s, F4s);

paramSpacings = [df,df,df,df];

% Deinfe scan_params 
delta = 21%
Delta = 35.5;
b = 3000%;
scan_params = [delta, Delta, b];

FI = FisherInfoGaussian(paramIndx, paramGrid, paramSpacings, scan_params)

SNR = exp(-0.002*b)/0.05;


