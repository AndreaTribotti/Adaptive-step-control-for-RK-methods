function h = initial_step_selector(fun, x0, y0, p, atol, rtol)
% INITIAL_STEP_SELECTOR Selects an initial step size for an ODE solver.
%
% Syntax:
%   h = initial_step_selector(fun, x0, y0, p, atol, rtol)
%
% Description:
%   Implements the algorithm by Gladwell, Shampine & Brankin (1987) for
%   choosing an initial step size, as described in Hairer & Nörsett,
%   Solving Ordinary Differential Equations I, page 169.
%
% Input:
%   fun  - Handle to the ODE function, e.g., @brusselator.
%   x0   - Initial time.
%   y0   - Initial value vector.
%   p    - Order of the numerical method to be used.
%   atol - Absolute tolerance (scalar).
%   rtol - Relative tolerance (scalar).
%
% Output:
%   h    - The suggested initial step size.

    % Ensure y0 is a column vector
    y0 = y0(:);

    % (a) Initial evaluations
    f0 = fun(x0, y0);
    
    % Scaling factor for norms
    sc = atol + abs(y0) .* rtol;
    
    d0 = vec_norm(y0, sc);
    d1 = vec_norm(f0, sc);

    % (b) First guess for the step size
    if d0 < 1e-5 || d1 < 1e-5
        h0 = 1e-6;
    else
        h0 = 0.01 * (d0 / d1);
    end

    % (c) Perform one explicit Euler step
    y1 = y0 + h0 * f0;
    f1 = fun(x0 + h0, y1);

    % (d) Estimate second derivative
    d2 = vec_norm(f1 - f0, sc) / h0;

    % (e) Compute a second step size guess h1
    max_d = max(d1, d2);
    if max_d <= 1e-15
        h1 = max(1e-6, h0 * 1e-3);
    else
        h1 = (0.01 / max_d)^(1 / (p + 1));
    end

    % (f) Final proposed step size
    h = min(100 * h0, h1);

end
