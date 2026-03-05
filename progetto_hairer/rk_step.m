function y_next = rk_step(fun, x, y, h)
% RK_STEP Performs a single step of the classic 4th order Runge-Kutta method.
%
% Syntax:
%   y_next = rk_step(fun, x, y, h)
%
% Description:
%   Integrates the ODE defined by 'fun' for a single step 'h' from the
%   current state (x, y) using the well-known RK4 method. This function
%   serves as a basic building block for more complex solvers.
%
% Input:
%   fun  - Handle to the ODE function, e.g., @brusselator.
%   x    - Current time.
%   y    - Current value vector.
%   h    - Step size.
%
% Output:
%   y_next - The solution vector at time x + h.

    k1 = fun(x, y);
    k2 = fun(x + 0.5 * h, y + 0.5 * h * k1);
    k3 = fun(x + 0.5 * h, y + 0.5 * h * k2);
    k4 = fun(x + h, y + h * k3);
    
    y_next = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

end
