
% =========================================================================
% butcher_matrices.m
% =========================================================================
% This file defines the Butcher tableaux for various embedded Runge-Kutta
% methods. The desired method can be selected by changing the value of
% the 'default_method' variable.
%
% Available methods:
%   - 'ERK4(3)': A 4th order method with a 3rd order embedded error estimate.
%   - 'BS3(2)': Bogacki-Shampine 3(2) method.
%   - 'DP5(4)': Dormand-Prince 5(4) method.
% =========================================================================

% --- Method Selection ---
% Change this value to select the default RK method to be used.
default_method = 'BS3(2)'; 

% --- Butcher Tableaux Definitions ---

methods = struct();

% Method 1: ERK4(3) - 4th order method with 3rd order embedded error estimate
methods(1).name = 'ERK4(3)';
methods(1).A = [0, 0, 0, 0;
             1/3, 0, 0, 0;
             -1/3, 1, 0, 0;
             1, -1, 1, 0];
methods(1).c = [0, 1/3, 2/3, 1];
methods(1).b = [1/8, 3/8, 3/8, 1/8];
methods(1).b_hat = [1/12, 1/2, 1/4, 0, 1/6]; % Note: Last coefficient is for k_next (FSAL)
methods(1).p = 4;
methods(1).q = 3;

% Method 2: Bogacki-Shampine 3(2)
methods(2).name = 'BS3(2)';
methods(2).A = [0, 0, 0, 0;
             1/2, 0, 0, 0;
             0, 3/4, 0, 0;
             2/9, 1/3, 4/9, 0];
methods(2).c = [0, 1/2, 3/4, 1];
methods(2).b = [2/9, 1/3, 4/9, 0];
methods(2).b_hat = [7/24, 1/4, 1/3, 1/8];
methods(2).p = 3;
methods(2).q = 2;

% Method 3: Dormand-Prince 5(4) - This is the method used in MATLAB's ode45
methods(3).name = 'DP5(4)';
methods(3).c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
methods(3).A = [0, 0, 0, 0, 0, 0, 0;
             1/5, 0, 0, 0, 0, 0, 0;
             3/40, 9/40, 0, 0, 0, 0, 0;
             44/45, -56/15, 32/9, 0, 0, 0, 0;
             19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
             9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0;
             35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
methods(3).b = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
methods(3).b_hat = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
methods(3).p = 5;
methods(3).q = 4;
