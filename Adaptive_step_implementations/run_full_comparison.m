% =========================================================================
% run_full_comparison.m
%
% This script performs a comparison between different Runge-Kutta methods,
% including both adaptive-step and fixed-step methods.
%
% The objective is to generate an efficiency graph (work-precision diagram)
% that plots the global error against the number of function evaluations.
%
% Adaptive Methods (based on embedded_solver.m):
% 1. BS3(2) - Bogacki-Shampine 3(2)
% 2. ERK4(3) - Custom 4(3) pair
% 3. DP5(4) - Dormand-Prince 5(4)
%
% Fixed-Step Methods:
% 4. Richardson Extrapolation (based on RK4)
% 5. ERK3 (fixed-step implementation)
% 6. ERK4 (fixed-step implementation)
% =========================================================================

clear; close all; clc;

% --- Problem and Solver Settings ---

fprintf('Problem setup...\n');
butcher_matrices; % Load Butcher matrices

% Problem data
fun = @brusselator;
t_span = [0, 20];
y0 = [1.5; 3];

% --- Calculation of the Reference Solution ---
fprintf('Calculating reference solution with ode45...\n');
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
sol_ref = ode45(fun, t_span, y0, opts);
y_ref_final = deval(sol_ref, t_span(2));
fprintf('Reference solution calculated.\n');

% Structure to store results
results = struct();

% --- Part 1: Testing Adaptive Methods ---
fprintf('--- Starting Adaptive Methods Test ---\n');

methods_adaptive = {'BS3(2)', 'ERK4(3)', 'DP5(4)'};
tols = 10.^(-4:-1:-10);

for i = 1:length(methods_adaptive)
    method_name = methods_adaptive{i};
    fprintf('Testing adaptive method: %s\n', method_name);
    
    global_errors = zeros(1, length(tols));
    function_evals = zeros(1, length(tols));
    
    % Find the correct tableau
    idx = find(strcmp({methods.name}, method_name));
    tableau = methods(idx);

    for j = 1:length(tols)
        tol = tols(j);
        [y_final, fevals] = solve_adaptive_embedded(fun, t_span, y0, tol, tol, tableau);
        global_errors(j) = norm(y_final - y_ref_final);
        function_evals(j) = fevals;
    end
    
    results.(matlab.lang.makeValidName(method_name)).errors = global_errors;
    results.(matlab.lang.makeValidName(method_name)).fevals = function_evals;
    results.(matlab.lang.makeValidName(method_name)).type = 'adaptive';
end


% --- Part 2: Testing Fixed-Step Methods ---
fprintf('--- Starting Fixed-Step Methods Test ---\n');

% For ERK3 and ERK4, we need their matrices (without the embedded part)
b3=[1/4,0,3/4]; c3=[0,1/3,2/3]; 
A3=[0 0
   1/3 0
   0 2/3];

erk4_idx = find(strcmp({methods.name}, 'ERK4(3)'));
erk4_b = methods(erk4_idx).b;
erk4_c = methods(erk4_idx).c;
erk4_A = methods(erk4_idx).A;

methods_fixed = {'Richardson', 'ERK3', 'ERK4'};
step_counts = [100, 200, 500, 750, 1000, 1250, 1500];

for i = 1:length(methods_fixed)
    method_name = methods_fixed{i};
    fprintf('Testing fixed-step method: %s\n', method_name);

    global_errors = zeros(1, length(step_counts));
    function_evals = zeros(1, length(step_counts));
    
    for j = 1:length(step_counts)
        n_steps = step_counts(j);
        h = (t_span(2) - t_span(1)) / n_steps;
        
        y_final = zeros(size(y0));
        fevals = 0;
        
        switch method_name
            case 'Richardson'
                y_current = y0;
                for k = 1:n_steps
                    [y_current, ~] = richardson_solver(fun, t_span(1) + (k-1)*h, y_current, h);
                end
                y_final = y_current;
                fevals = 12 * n_steps; % 8 (fine) + 4 (coarse)
                
            case 'ERK3'
                [~, u, fevals_erk] = ERK3(fun, t_span, y0, h, b3, c3, A3);
                y_final = u(:, end);
                fevals = fevals_erk;

            case 'ERK4'
                [~, u, fevals_erk] = ERK4(fun, t_span, y0, h, erk4_b, erk4_c, erk4_A);
                y_final = u(:, end);
                fevals = fevals_erk;
        end
        
        global_errors(j) = norm(y_final - y_ref_final);
        function_evals(j) = fevals;
    end
    
    results.(matlab.lang.makeValidName(method_name)).errors = global_errors;
    results.(matlab.lang.makeValidName(method_name)).fevals = function_evals;
    results.(matlab.lang.makeValidName(method_name)).type = 'fixed';
end


% --- Part 3: Plot Generation ---
figure;
hold on;
grid on;

method_names_plot = [methods_adaptive, methods_fixed];
colors = lines(length(method_names_plot)); % Generate distinct colors

for i = 1:length(method_names_plot)
    method_name = method_names_plot{i};
    valid_name = matlab.lang.makeValidName(method_name);
    
    res = results.(valid_name);
    valid_indices = res.errors > 0;
    
    if any(valid_indices)
        if strcmp(res.type, 'adaptive')
            % Plot adaptive with line
            loglog(res.errors(valid_indices), res.fevals(valid_indices), ...
                   '-', 'DisplayName', method_name, 'LineWidth', 1.5, 'Color', colors(i,:));
        else
            % Plot fixed with dashed line
            loglog(res.errors(valid_indices), res.fevals(valid_indices), ...
                   '--', 'DisplayName', [method_name, ' (fixed)'], 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
end

% Plot settings
title('Efficiency Comparison: Adaptive vs. Fixed-Step Methods');
xlabel('Global Error');
ylabel('Number of Function Evaluations');
legend('show', 'Location', 'northwest');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XDir', 'reverse');

% =========================================================================
% LOCAL SOLVER FUNCTIONS (reused from before)
% =========================================================================

function [y_final, total_fevals] = solve_adaptive_embedded(fun, t_span, y0, rtol, atol, tableau)
    % Step size control parameters
    fac = 0.8; facmin = 0.1; facmax = 5.0; max_steps = 50000;

    % Initialization
    t_start = t_span(1); t_end = t_span(2);
    t_current = t_start; y_current = y0;
    
    % Initial step selector
    p = tableau.p;
    h = initial_step_selector(fun, t_current, y_current, p, atol, rtol);
    total_fevals = 2; % initial_step_selector makes 2 evaluations

    k1_fsal = []; % k1 for the first step (empty)
    step_accepted = false; % Flag to track if the previous step was accepted
    
    num_stages = length(tableau.c);
    is_fsal = isfield(tableau, 'fsal') && tableau.fsal;

    while t_current < t_end
        if t_current + h > t_end, h = t_end - t_current; end

        % For an FSAL method, k1 is reused only if the previous step was accepted
        if is_fsal && step_accepted
            k1_to_pass = k1_fsal;
            evals_this_step = num_stages - 1; % One evaluation is saved
        else
            k1_to_pass = [];
            evals_this_step = num_stages; % All evaluations are performed
        end
        
        % Perform one step and calculate k_next (which costs 1 extra evaluation)
        [y_high, ~, k_next, err] = embedded_solver(fun, t_current, y_current, h, k1_to_pass, atol, rtol, tableau);
        evals_this_step = evals_this_step + 1; % Add the evaluation for k_next

        if err < 1
            % Step accepted
            t_current = t_current + h;
            y_current = y_high;
            if is_fsal
                k1_fsal = k_next; % Save k_next for the next step
            end
            step_accepted = true;
            facmax_current = facmax;
        else
            % Step rejected
            step_accepted = false;
            facmax_current = 1.0;
        end
        
        total_fevals = total_fevals + evals_this_step;
        
        % Adjust the step for the next attempt
        if t_current < t_end
            h = step_control(h, err, tableau.q, fac, facmin, facmax_current);
        end
    end
    y_final = y_current;
end
