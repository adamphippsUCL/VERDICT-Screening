function [scan_params, CN] = minimiseCN(nscan, opts)

% Function to find set of scan parameters which minimise matrix condition
% number

arguments
    nscan % Number of scans

    opts.bvals = [] % b values of scans
    opts.bmin = 500;
    opts.bmax = 2000;

    % Bounds
    opts.deltaMin = 0
    opts.deltaMax = 40
    opts.DeltaMin = 0
    opts.DeltaMax = 50

    % Max gradient strength
    opts.Gmax = 60
    opts.t180 = 2 % 180 pulse duration

end


if isempty(opts.bvals)
    optimiseb = true;
else
    optimiseb = false;
end


%% Construct x0 vector

delta0 = 25;
delta0var = 5;
Delta0 = 40;
Delta0var = 5;

b0min = 500;
b0max = 2000;

switch optimiseb

    case false

        % Build x0 vector
        x0 = [];
        for scanIndx = 1:nscan
            x0 = [...
                x0,...
                delta0 - delta0var + 2*delta0var*rand(),...
                Delta0 - Delta0var + 2*Delta0var*rand(),...
                opts.bvals(scanIndx)...
                ];
        end


    case true

        % Build x0 vector
        x0 = [];
        for scanIndx = 1:nscan
            x0 = [...
                x0,...
                delta0 - delta0var + 2*delta0var*rand(),...
                Delta0 - Delta0var + 2*Delta0var*rand(),...
                b0min + (b0max-b0min)*rand()...
                ];
        end
end

x0 = transpose(x0);


%% Non-linear optimisation 

% Using 'fmincon'

% == Define function to optimise
fun = @conditionNumber;

% == Define bounds
deltaMin = 0;
deltaMax = 40;
DeltaMin = 0;
DeltaMax = 50;

% Define bound vectors
lb = [];
ub = [];
for scanIndx = 1:nscan

    switch optimiseb

        case false
            lb = [lb, deltaMin, DeltaMin, opts.bvals(scanIndx)];
            ub = [ub, deltaMax, DeltaMax, opts.bvals(scanIndx)];

        case true

            lb = [lb, deltaMin, DeltaMin, opts.bmin];
            ub = [ub, deltaMax, DeltaMax, opts.bmax];

    end
end

% == Equality constraint (None)

Aeq = zeros(nscan,3*nscan);
beq = zeros(nscan,1);


% == Linear constraint (delta < Delta)

% Define matrix
A = zeros(nscan,3*nscan);

for rowIndx = 1:nscan
    A(rowIndx, ((rowIndx-1)*3 + 1):rowIndx*3) = [1, -1, 0];
end

% Define bounding b vector
b = zeros(nscan,1);
for scanIndx = 1:nscan
    b(scanIndx) = -opts.t180;
end


% == Non-linear constraint (Gradient strength)

nonlcon = @gradnonlcon;


% == Apply optimisation

[scan_params, CN] = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon);









end