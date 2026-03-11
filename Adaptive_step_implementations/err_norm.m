function err = err_norm(y_high, y_low, y_prev, atol, rtol)
% ERR_NORM Computes the weighted root-mean-square error norm.
%
% Syntax:
%   err = err_norm(y_high, y_low, y_prev, atol, rtol)
%
% Description:
%   Calculates the error norm based on formula (4.11) from Hairer & Nörsett.
%   This norm is used to measure the local error in numerical methods for ODEs.
%   The error is estimated as the difference between a higher-order solution
%   (y_high) and a lower-order solution (y_low).
%
% Input:
%   y_high - Higher-order solution vector at the end of the step.
%   y_low  - Lower-order solution vector at the end of the step.
%   y_prev - Solution vector at the beginning of the step.
%   atol   - Absolute tolerance (scalar).
%   rtol   - Relative tolerance (scalar).
%
% Output:
%   err - The computed scalar error norm.

    % Calculate scale for each component (formula 4.10)
    sc = atol + max(abs(y_prev), abs(y_high)) .* rtol;
    
    % Calculate the squared error for each component, weighted by the scale
    err_sq = ((y_high - y_low) ./ sc).^2;
    
    % Compute the root-mean-square of the weighted errors (formula 4.11)
    err = sqrt(mean(err_sq));
    
end
