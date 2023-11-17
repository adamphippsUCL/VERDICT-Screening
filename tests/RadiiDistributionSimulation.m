% MATLAB script to simulate signal from given distribution of radii

% Simulate signal from distribution of spheres over VERDICT scheme, fit to
% signal (write my own fitting script just for radii) and assess fitting 
% performance.
% First define range of radii
Radii = linspace(0.1, 15.1, 20);

%% Define distribution

% % Normal distribution
% mu = 8;
% sigma = 3;
% fRs = normpdf(Radii, mu, sigma);

% Uniform distribution
fRs = (1/length(Radii))*ones(length(Radii), 1);


%% Simulate signal using sphereGPD


figure;
plot(Radii, fRs)
