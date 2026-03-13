function h_new = step_control(h, err, q, fac, facmin, facmax)
% STEP_CONTROL Implements the step size control formula.
%
% Syntax:
%   h_new = step_control(h, err, q, fac, facmin, facmax)
%
% Description:
%   Calculates the new step size 'h_new' based on the current step size 'h'
%   and the estimated error 'err'. The formula is taken from Hairer & Nörsett,
%   Solving Ordinary Differential Equations I, page 168, formula (4.13).
%
% Input:
%   h      - The current step size.
%   err    - The estimated local error (scalar, computed by a norm).
%   q      - Order of the error estimator (p in the book, where error is O(h^(p+1))).
%   fac    - Safety factor (e.g., 0.8 or 0.9).
%   facmin - Minimum factor for step size change.
%   facmax - Maximum factor for step size change.
%
% Output:
%   h_new  - The proposed new step size.

    % Formula (4.13)
    % h_new = h * min(facmax, max(facmin, fac * (1/err)^(1/(q+1))))
    
    if err == 0
        % If error is zero, increase step size by maximum factor
        h_new = h * facmax;
    else
        % Avoid division by zero and apply the formula
        step_ratio = fac * (1 / err)^(1 / (q + 1));
        h_new = h * min(facmax, max(facmin, step_ratio));
    end

end
