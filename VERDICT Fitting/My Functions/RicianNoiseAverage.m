% MATLAB function to return average signal magnitude after application of
% Rician noise with a given SNR.

function signal = RicianNoiseAverage(signal, NoiseMagnitude)

% Construct Rician distribution
RiceDist = makedist('Rician','s',signal,'sigma',NoiseMagnitude);

signal = mean(RiceDist);


end