
function scheme = BuildScheme(scans)
% Build a scheme structure
% scans(indx, delta, Delta, bval)


for scanIndx = 1:size(scans,1)

    delta = scans(scanIndx, 1);
    Delta = scans(scanIndx,2);
    bval = scans(scanIndx,3);

    scheme(scanIndx).delta = delta;
    scheme(scanIndx).DELTA = Delta ;
    scheme(scanIndx).bval = bval;
    scheme(scanIndx).G = stejskal(delta,Delta,bval=bval);
end
    
end