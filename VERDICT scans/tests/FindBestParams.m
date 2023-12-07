% MATLAB test script to find optimal parameters 

% Choose b values
bs = [500,1500,2000];

% number of scans
nscan = length(bs);

% Set default delta and Delta
delta0 = 25;
Delta0 = 40;


% Build x0 vector
x0 = [];
for scanIndx = 1:nscan
    x0 = [x0, delta0, Delta0, bs(scanIndx)];
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
    lb = [lb, deltaMin, DeltaMin, bs(scanIndx)];
    ub = [ub, deltaMax, DeltaMax, bs(scanIndx)];
end

% == Equality constraint
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


% == Non-linear constraint (Gradient strength)

nonlcon = @gradstrength;

% Try it...
fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon)





function [c, ceq] = gradstrength(x, opts)

% Function to determine if gradient strength is below maximum
% result < 0 for each scan indicates gradient strength is ok

    arguments
        x % scan parameter vector
        opts.Gmax = 60*10^(-6) % T/mm
        opts.gamma = 2.675221874*(10^8)% s/T* 1e+08;
    end
    
    nscan = (1/3)*length(x);
    
    c = zeros(nscan,1);
  
    for scanIndx =1:nscan

        b = x(3*(scanIndx-1)+3);
        delta = x(3*(scanIndx-1)+1)*(10^(-3));
        Delta = x(3*(scanIndx-1)+2)*(10^(-3));

        c(scanIndx) = (b/(opts.Gmax^2)) - ((opts.gamma*delta)^2)*(Delta - (1/3)*delta);
    end

    ceq = 0;
end


