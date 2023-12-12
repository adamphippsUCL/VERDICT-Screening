function [c, ceq] = gradnonlcon(x, opts)

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


