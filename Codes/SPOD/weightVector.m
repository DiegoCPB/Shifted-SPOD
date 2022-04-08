function [wVec] = weightVector(X)
%   Function to calculate the weightvector os a 1D grid 
%   using the integration trapezoidal rule.
    
    dX = sqrt((X(2:end,:)-X(1:end-1,:)).^2);
    
    wX = [ 2*dX(1,:); dX(2:end,:)+dX(1:end-1,:); 2*dX(end,:)];
    
    wVec = abs(wX)/2;
end