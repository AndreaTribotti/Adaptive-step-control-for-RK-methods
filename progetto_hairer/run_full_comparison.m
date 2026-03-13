% =========================================================================
% run_full_comparison.m
%
% This script performs a comparison between different Runge-Kutta methods,
% including both adaptive-step and fixed-step methods.
%
% The objective is to generate an efficiency graph (work-precision diagram)
% that plots the global error against the number of function evaluations
% and a comparison table.
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
erk3_idx = find(strcmp({methods.name}, 'BS3(2)'));
erk3_b = methods(erk3_idx).b;
erk3_c = methods(erk3_idx).c;
erk3_A = methods(erk3_idx).A;

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
                [~, u, fevals_erk] = ERK3(fun, t_span, y0, h, erk3_b, erk3_c, erk3_A);
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


%% Comparison table

fprintf("\n")
n_steps_fixed = 600;
tol_adaptive = 1e-7;

results_table = {'Method', 'h_min', 'h_max', 'FEvals', 'Error', 'Total Steps', 'Rejected Steps'};
method_names = {'ERK3', 'ERK4', 'Richardson', 'BS3(2)', 'ERK4(3)', 'DP5(4)'};

% --- Run Solvers and Collect Data ---
for i = 1:length(method_names)
    method_name = method_names{i};
    
    h_min_str = 'N/A';
    h_max_str = 'N/A';
    
    switch method_name
        case {'ERK3', 'ERK4'}
            % --- Fixed-Step ERK Methods ---
            h = (t_span(2) - t_span(1)) / n_steps_fixed;
            
            if strcmp(method_name, 'ERK3')
                erk_idx = find(strcmp({methods.name}, 'BS3(2)')); % Use BS32 table for ERK3
                tableau = methods(erk_idx);
                [~, u, fevals] = ERK3(fun, t_span, y0, h, tableau.b, tableau.c, tableau.A);
            else % ERK4
                erk_idx = find(strcmp({methods.name}, 'ERK4(3)')); % Use ERK43 table for ERK4
                tableau = methods(erk_idx);
                [~, u, fevals] = ERK4(fun, t_span, y0, h, tableau.b, tableau.c, tableau.A);
            end
            
            y_final = u(:, end);
            err = norm(y_final - y_ref_final);
            total_steps = n_steps_fixed;
            rejected_steps = 0;

        case 'Richardson'
            % --- Richardson Extrapolation ---
            h = (t_span(2) - t_span(1)) / n_steps_fixed;
            y_current = y0;
            for k = 1:n_steps_fixed
                t_current = t_span(1) + (k-1)*h;
                [y_current, ~] = richardson_solver(fun, t_current, y_current, h);
            end
            y_final = y_current;
            fevals = 11 * n_steps_fixed; % From richardson_solver implementation
            err = norm(y_final - y_ref_final);
            total_steps = n_steps_fixed;
            rejected_steps = 0;
            
        case {'BS3(2)', 'ERK4(3)', 'DP5(4)'}
            % --- Adaptive-Step Methods ---
            idx = find(strcmp({methods.name}, method_name));
            tableau = methods(idx);
            
            [y_final, fevals, total_steps, rejected_steps, h_min, h_max] = ...
                solve_adaptive_with_stats(fun, t_span, y0, tol_adaptive, tol_adaptive, tableau);
            
            err = norm(y_final - y_ref_final);
            h_min_str = sprintf('%.2e', h_min);
            h_max_str = sprintf('%.2e', h_max);
    end
    
    % Store results
    results_table(i+1, :) = {method_name, h_min_str, h_max_str, fevals, sprintf('%.2e', err), total_steps, rejected_steps};
end

% --- Display Results Table ---
fprintf('Generating final comparison table:\n\n');
fprintf('%-15s | %-10s | %-10s | %-10s | %-12s | %-12s | %-15s\n', results_table{1,:});
fprintf('%s\n', repmat('-', 1, 105));
for i = 2:size(results_table, 1)
    fprintf('%-15s | %-10s | %-10s | %-10d | %-12s | %-12d | %-15d\n', results_table{i,:});
end
fprintf('\n');


% =========================================================================
% LOCAL FUNCTION for adaptive solving with statistics gathering
% =========================================================================
function [y_final, total_fevals, total_steps, rejected_steps, h_min, h_max] = solve_adaptive_with_stats(fun, t_span, y0, rtol, atol, tableau)
    % Step size control parameters
    fac = 0.8; facmin = 0.1; facmax = 5.0;

    % Initialization
    t_start = t_span(1); t_end = t_span(2);
    t_current = t_start; y_current = y0;
    
    % Statistics initialization
    total_steps = 0;
    rejected_steps = 0;
    h_min = inf;
    h_max = 0;
    
    % Initial step selector
    p = tableau.p;
    h = initial_step_selector(fun, t_current, y_current, p, atol, rtol);
    total_fevals = 2; % From initial_step_selector

    % FSAL setup
    num_stages = length(tableau.c);
    is_fsal = isfield(tableau, 'fsal') && tableau.fsal;
    k1_fsal = [];
    step_accepted_prev = false;

    while t_current < t_end
        if t_current + h > t_end, h = t_end - t_current; end

        % Track min/max step size
        h_min = min(h_min, h);
        h_max = max(h_max, h);
        
        % --- Perform one step ---
        if is_fsal && step_accepted_prev
            k1_to_pass = k1_fsal;
            fevals_step = num_stages - 1;
        else
            k1_to_pass = [];
            fevals_step = num_stages;
        end

        [y_next, ~, k_next, err] = embedded_solver(fun, t_current, y_current, h, k1_to_pass, atol, rtol, tableau);
        total_fevals = total_fevals + fevals_step;
        
        % The error 'err' is already scaled by err_norm.
        err_scaled = err;

        % --- Step acceptance/rejection ---
        if err_scaled < 1
            % ACCEPTED
            t_current = t_current + h;
            y_current = y_next;
            total_steps = total_steps + 1;
            step_accepted_prev = true;
            if is_fsal, k1_fsal = k_next; end
            
            % Use a more aggressive facmax for accepted steps
            facmax_current = facmax;
        else
            % REJECTED
            rejected_steps = rejected_steps + 1;
            step_accepted_prev = false;
            
            % Use a facmax of 1 for rejected steps to avoid overly large step increases
            facmax_current = 1.0;
        end
        
        % --- Step size control for next step ---
        h = step_control(h, err_scaled, tableau.q, fac, facmin, facmax_current);
    end
    
    y_final = y_current;
end

