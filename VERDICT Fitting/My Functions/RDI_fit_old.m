function [fIC, fEES, rmse] = RDI_fit(scheme,Y,t,opt)

% Perform Restricted Diffusion Imaging (RDI) fit for fIC and fEES
% 
% [fIC, fEES] = RDI_fit(scheme, Y, t, opt)
% 
% RDI fitting assumes normal cellular radii distribution as seen from histology 
% patches and only has two fitting parameters fIC and fEES
%     -> muR = 8µm
%     -> sigmaR = 2µm

% L2 Tikonhov regularisation used to penalise extreme measurements


arguments 
    scheme (1,:)  struct
    Y {mustBeNumeric,mustBeReal}  % [ny nx nmeas] or [ny nx nz nmeas] 
    t.dEES  (1,1) {mustBeNumeric,mustBeReal} = 2
    t.dIC   (1,1) {mustBeNumeric,mustBeReal} = 2
    t.muR (1,1) {mustBeNumeric,mustBeReal} = 10 % µm
    t.sigmaR (1,1) {mustBeNumeric,mustBeReal} = 4 % µm
    opt.mask {mustBeNumeric,mustBeReal} = []
    opt.solver = 'lsqnonnegTikhonov'
    opt.NoisyMatrix = false;
    opt.NoiseSigma = 0.05;

end

% Scheme, Data, and Image sizes
nscheme = length(scheme) ;
szY = size(Y) ;
szmap = szY(1:end-1) ;


if ~exist('opt','var') || isempty(opt.mask)
    opt.mask = ones(szmap) ;
end

if szY(end) ~= nscheme
    warning('MATLAB:verdict_fit:dataSizeInconsistencies', ...
        ['Number of measures (',num2str(szY(end)), ...
        ') must be same as size of schemes (', num2str(nscheme),')'])
end


% Number of unknowns
nUnknown = 2;


%% Construct Matrix

A = zeros([nscheme nUnknown]);
sIC = zeros([1 nscheme]);
sEES = zeros([1 nscheme]);


% == Define radii distribution
nR = 50;
minR = 0.1;
maxR = 15.1;
Rs = linspace(minR, maxR, nR);
muR = t.muR;
sigmaR = t.sigmaR;

fRs = normpdf(Rs, muR, sigmaR);
% Normalise
normfRs = fRs/sum(fRs);


% Simulate signal

for ischeme = 1:nscheme

    % === IC signal

    ICsignal = 0;

    for Rindx = 1:nR

        Rsignal = normfRs(Rindx)*sphereGPD(scheme(ischeme).delta, scheme(ischeme).DELTA, ...
           scheme(ischeme).G, Rs(Rindx), t.dIC);

        ICsignal = ICsignal + Rsignal;

    end

    % Noise 
    if opt.NoisyMatrix && (scheme(ischeme).bval ~= 0)
        ICsignal = RicianNoiseAverage(ICsignal, opt.NoiseSigma);
    end

    sIC(ischeme) = ICsignal;


    % === EES signal

    EESsignal = ball(scheme(ischeme).bval, t.dEES);

    % Noise
    if (opt.NoisyMatrix) && (scheme(ischeme).bval ~= 0)
        EESsignal = RicianNoiseAverage(EESsignal, opt.NoiseSigma);
    end

    sEES(ischeme)  = EESsignal;

    % Append to matrix
    A(ischeme,:) = [sIC(ischeme) sEES(ischeme)];

end

disp(A)
%% Fitting

fIC = zeros(szmap);
fEES = zeros(szmap);
rmse = zeros(szmap);

% Flatten
opt.mask = opt.mask(:) ;
fIC = fIC(:); fEES = fEES(:);
Y = reshape(Y,[prod(szmap) nscheme]) ;



locinmask = find(opt.mask) ;

for ip = 1:length(locinmask)  % can be parallelised

    ind = locinmask(ip) ;
    y = Y(ind,:) ;

    percent_done = (ip/length(locinmask))*100
    
    % Remove infinties and NaN
    y(y==inf) = 0;
    y(y==-inf) = 0;
    y(isnan(y)) = 0;

    switch opt.solver

        case 'lsqnonnegTikhonov'
    
            nUnknown = size(A,2) ;
            lambda = sqrt(0.00) ;
    
    
            L = eye(nUnknown) ;
    
            AT = [A; lambda*L] ;
            YT = [y(:) ; zeros([nUnknown 1])] ;

            x = lsqnonneg(AT, YT) ;


    fIC(ind) = x(1);
    fEES(ind) = x(2);
    rmse(ind) = norm(y(:) - A*x) / sqrt(nscheme) ;

    end
end

fIC = reshape(fIC,szmap) ;
fEES = reshape(fEES,szmap) ;
rmse = reshape(rmse, szmap) ;

end

