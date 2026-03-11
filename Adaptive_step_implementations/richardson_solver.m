function [y_extrapolated, y_coarse] = richardson_solver(fun, x0, y0, H)
% RICHARDSON_SOLVER Performs a single step using Richardson extrapolation.
%
% Syntax:
%   [y_extrapolated, y_coarse] = richardson_solver(fun, x0, y0, H)
%
% Description:
%   Implements Richardson Extrapolation for a single step 'H', based on an
%   underlying Runge-Kutta method of order p. It follows Theorem 4.1 from
%   Hairer & Nörsett. The base method used is the classic 4th order
%   Runge-Kutta method implemented in 'rk_step.m'.
%
%   The solver computes:
%   1. y_fine: Solution at x0 + H by taking two small steps of size h = H/2.
%   2. y_coarse: Solution at x0 + H by taking one big step of size H.
%   3. y_extrapolated: The extrapolated solution of order p+1.
%
% Input:
%   fun  - Handle to the ODE function, e.g., @brusselator.
%   x0   - Initial time.
%   y0   - Initial value vector.
%   H    - The total step size to perform (the "2h" in the textbook).
%
% Output:
%   y_extrapolated - The extrapolated solution of order p+1 at time x0 + H.
%   y_coarse       - The coarse-grid solution of order p at time x0 + H.

    % The order of the underlying method (rk_step is classic RK4)
    p = 4;
    
    % h is the small step size
    h = H / 2;

    % 1. Compute fine-grid solution (y2 in the book)
    % Two steps with size h
    y_temp = rk_step(fun, x0, y0, h);
    y_fine = rk_step(fun, x0 + h, y_temp, h);

    % 2. Compute coarse-grid solution (w in the book)
    % One step with size H = 2h
    y_coarse = rk_step(fun, x0, y0, H);
    
    % 3. Apply Richardson Extrapolation (formula 4.5)
    % y_hat = y_fine + (y_fine - y_coarse) / (2^p - 1)
    y_extrapolated = y_fine + (y_fine - y_coarse) / (2^p - 1);
    
end
