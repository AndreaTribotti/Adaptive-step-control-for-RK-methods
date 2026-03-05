function [y_high, y_low, k_next, err] = embedded_solver(fun, x, y, h, k1, atol, rtol)
% EMBEDDED_SOLVER Performs a single step of an embedded 4(3) RK method.
%
% Syntax:
%   [y_high, y_low, k_next, err] = embedded_solver(fun, x, y, h, k1, atol, rtol)
%
% Description:
%   Implements a single step of a 4(3) Runge-Kutta method with FSAL
%   (First Same As Last). The higher-order method is the classic "3/8 Rule"
%   of order 4. The embedded method of order 3 is derived from formula (4.9)
%   in Hairer & Nörsett.
%
% Input:
%   fun  - Handle to the ODE function, e.g., @brusselator.
%   x    - Current time.
%   y    - Current value vector.
%   h    - Step size.
%   k1   - The result of the first function evaluation (from the previous step
%          due to FSAL). Can be empty [] for the very first step.
%   atol - Absolute tolerance for error calculation.
%   rtol - Relative tolerance for error calculation.
%
% Output:
%   y_high - 4th order solution at time x + h.
%   y_low  - 3rd order solution (for error estimation) at time x + h.
%   k_next - The function evaluation at (x+h, y_high), to be used as k1 in
%            the next step (FSAL).
%   err    - The estimated local error for this step.

    % Butcher Tableau for the 3/8 rule (Order 4)
    c = [0, 1/3, 2/3, 1];
    A = [0, 0, 0, 0;
         1/3, 0, 0, 0;
         -1/3, 1, 0, 0;
         1, -1, 1, 0];
    b = [1/8, 3/8, 3/8, 1/8];

    % Coefficients for the embedded method (Order 3), from formula (4.9)
    % b_hat = [1/12, 1/2, 1/4, 0, 1/6]
    b_hat = [1/12, 1/2, 1/4, 0]; % Last one is for k5

    % --- Stage computations ---
    if isempty(k1)
        k1 = fun(x, y);
    end
    
    k2 = fun(x + c(2)*h, y + h*(A(2,1)*k1));
    k3 = fun(x + c(3)*h, y + h*(A(3,1)*k1 + A(3,2)*k2));
    k4 = fun(x + c(4)*h, y + h*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3));

    % --- Compute solutions ---
    % Higher-order solution (y1 in book)
    y_high = y + h * (b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4);
    
    % FSAL stage: k5 is the evaluation for the next step
    k_next = fun(x + h, y_high);
    
    % Lower-order solution (y_hat in book)
    y_low = y + h * (b_hat(1)*k1 + b_hat(2)*k2 + b_hat(3)*k3 + b_hat(4)*k4 + (1/6)*k_next);
    
    % --- Error estimation ---
    % The error is estimated using the norm of the difference
    err = err_norm(y_high, y_low, y, atol, rtol);
    
end
