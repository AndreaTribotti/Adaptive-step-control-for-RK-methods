function y_next = gauss_legendre_solver(fun, x0, y0, h)
% GAUSS_LEGENDRE_SOLVER Performs a single step of the 2-stage Gauss-Legendre method.
%
% Syntax:
%   y_next = gauss_legendre_solver(fun, x0, y0, h)
%
% Description:
%   Implements the 2-stage, 4th-order implicit Runge-Kutta method based on
%   Gauss-Legendre coefficients. The internal stages are computed by solving
%   a system of nonlinear equations using MATLAB's 'fsolve'.
%
% Input:
%   fun  - Handle to the ODE function, e.g., @brusselator.
%   x0   - Initial time.
%   y0   - Initial value vector.
%   h    - Step size.
%
% Output:
%   y_next - The solution vector at time x0 + h.

    % Ensure y0 is a column vector
    y0 = y0(:);
    n = length(y0);

    % Butcher Tableau for Gauss-Legendre (2-stage, 4th order)
    c = [1/2 - sqrt(3)/6, 1/2 + sqrt(3)/6];
    A = [1/4, 1/4 - sqrt(3)/6;
         1/4 + sqrt(3)/6, 1/4];
    b = [1/2, 1/2];
    s = 2; % Number of stages

    % --- Solve for the internal stages k_i ---
    
    % Initial guess for the stage vectors k_i
    k_guess1 = fun(x0, y0);
    k_guess2 = fun(x0, y0);
    K_guess = [k_guess1; k_guess2]; % K is a stacked vector [k1; k2]
    
    % fsolve options
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-9);
    
    % Define the system of equations for fsolve as a local function
    % K is a stacked vector of size (n*s x 1)
    F = @(K) residual_system(K, fun, x0, y0, h, n, s, c, A);
    
    % Solve for the stages
    K_sol = fsolve(F, K_guess, options);
    
    % Reshape the solution back into k1, k2
    k1 = K_sol(1:n);
    k2 = K_sol(n+1:end);
    
    % --- Compute final solution ---
    y_next = y0 + h * (b(1) * k1 + b(2) * k2);

end

function Res = residual_system(K, fun, x0, y0, h, n, s, c, A)
    % RESIDUAL_SYSTEM Defines the nonlinear system F(K) = 0 for the stages.
    % K is a stacked vector [k1; k2; ...; ks].
    
    Res = zeros(n * s, 1);
    
    % Unstack K into a cell array or a matrix for easier handling
    k_stages = reshape(K, n, s);
    
    for i = 1:s
        % Calculate the argument for the function evaluation
        y_arg = y0;
        for j = 1:s
            y_arg = y_arg + h * A(i, j) * k_stages(:, j);
        end
        
        % Evaluate the function
        f_eval = fun(x0 + c(i) * h, y_arg);
        
        % Compute the residual for the i-th block of equations
        Res((i-1)*n + 1 : i*n) = k_stages(:, i) - f_eval;
    end
end
