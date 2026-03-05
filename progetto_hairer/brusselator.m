function dydt = brusselator(t, y)
% BRUSSELATOR Defines the Brusselator system of ODEs.
%
% Syntax:
%   dydt = brusselator(t, y)
%
% Description:
%   Implements the Brusselator system as defined in Hairer & Nörsett, 
%   Solving Ordinary Differential Equations I, page 170, formula (4.15).
%
%   The system is:
%   y1' = 1 + y1^2 * y2 - 4 * y1
%   y2' = 3 * y1 - y1^2 * y2
%
% Input:
%   t - Time (not used, as the system is autonomous).
%   y - A 2-element column vector [y1; y2] representing the state.
%
% Output:
%   dydt - A 2-element column vector representing the derivatives [y1'; y2'].

    dydt = zeros(2, 1);
    dydt(1) = 1 + y(1)^2 * y(2) - 4 * y(1);
    dydt(2) = 3 * y(1) - y(1)^2 * y(2);
end
