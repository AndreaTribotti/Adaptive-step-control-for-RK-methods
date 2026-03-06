function [y_high, y_low, k_next, err] = embedded_solver(fun, x, y, h, k1, atol, rtol, tableau)
% EMBEDDED_SOLVER Performs a single step of a generic embedded RK method.
%
% Syntax:
%   [y_high, y_low, k_next, err] = embedded_solver(fun, x, y, h, k1, atol, rtol, tableau)
%
% Description:
%   Implements a single step of an embedded Runge-Kutta method, which can
%   be specified via the 'tableau' parameter. It supports FSAL (First Same
%   As Last) functionality if the method has that property.
%
% Input:
%   fun     - Handle to the ODE function.
%   x       - Current time.
%   y       - Current value vector.
%   h       - Step size.
%   k1      - The result of the first function evaluation (from a previous
%             step for FSAL methods). Can be empty [].
%   atol    - Absolute tolerance.
%   rtol    - Relative tolerance.
%   tableau - Struct containing the Butcher tableau for the method. It must
%             contain fields A, c, b, and b_hat.
%
% Output:
%   y_high - Higher-order solution at time x + h.
%   y_low  - Lower-order solution at time x + h.
%   k_next - Function evaluation at (x+h, y_high), for FSAL.
%   err    - Estimated local error.

    A = tableau.A;
    c = tableau.c;
    b = tableau.b;
    b_hat = tableau.b_hat;
    
    s = size(A, 1);
    k = zeros(length(y), s);

    % --- Stage computations ---
    if isempty(k1)
        k(:, 1) = fun(x, y);
    else
        k(:, 1) = k1;
    end
    
    for i = 2:s
        y_stage = y;
        for j = 1:i-1
            y_stage = y_stage + h * A(i,j) * k(:,j);
        end
        k(:, i) = fun(x + c(i)*h, y_stage);
    end

    % --- Compute solutions ---
    y_high = y;
    for i = 1:length(b)
        y_high = y_high + h * b(i) * k(:,i);
    end
    
    % FSAL stage (if applicable)
    k_next = fun(x + h, y_high);
    
    y_low = y;
    for i = 1:length(b_hat)
        if i <= s
            y_low = y_low + h * b_hat(i) * k(:,i);
        else % This handles the FSAL part where b_hat has an extra element
            y_low = y_low + h * b_hat(i) * k_next;
        end
    end
    
    % --- Error estimation ---
    err = err_norm(y_high, y_low, y, atol, rtol);
    
end
