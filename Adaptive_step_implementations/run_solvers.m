clear; close all; clc;
addpath(pwd); % Ensure that files in the current folder are visible

% Load ODEs and analytical solutions
odes_examples

% Load Butcher tables
butcher_matrices

% Select the default method defined in the matrix file
selected_method_index = find(strcmp({methods.name}, default_method));
if isempty(selected_method_index)
    error('Default method "%s" not found.', default_method);
end
tableau = methods(selected_method_index);
fprintf('--- Selected Runge-Kutta method: %s ---\n\n', tableau.name);

% TEST AND VALIDATION SCRIPT
% ----------------------------------------------------
functions = {@brusselator, f1, f2};
solutions = {[], sol1, sol2}; % Added solution array
T = [20, 20, 2];
y0s = {[1.5; 3], 1, 5};
labels = {'Brusselator', 'y''=3+cos(t)-y', 'y''=y(2+y)(2-y)'};

% Common parameters for step control
atol = 1e-4;
rtol = 1e-4;
fac = 0.8;
facmin = 0.1;
facmax = 5.0;
max_steps = 5000;


% For simplicity, we use a fixed number of steps for non-adaptive solvers
n_steps_fixed = 200;

for j = 1:length(functions)
    fprintf('--- Calculation in progress for ODE: %s ---\n', labels{j})
    fun = functions{j};
    sol = solutions{j};
    t_span = [0, T(j)];
    x0 = t_span(1);
    y0 = y0s{j};
    H_fixed = (t_span(2) - t_span(1)) / n_steps_fixed;

    % Initial step calculation
    p_embedded = tableau.p; 
    h0 = initial_step_selector(fun, x0, y0, p_embedded, atol, rtol);
    fprintf('Suggested initial step: h = %.6f\n', h0);

    % If an analytical solution is not available, calculate a reference one
    if isempty(sol)
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        sol_struct = ode15s(fun, t_span, y0, options);
        sol = @(t) deval(sol_struct, t);
    end

    % --- Embedded Solver (Adaptive Step) ---
    fprintf('Running Embedded Solver (adaptive)...\n');
    t_current = x0;
    y_current = y0;
    h = h0;
    k1 = [];

    % Arrays for results
    t_embedded = t_current;
    y_embedded = y_current;
    
    % Vectors for logging
    h_history = [];
    t_history = [];
    is_accepted = [];
    err_history = [];
    global_err_history = [];
    exact_local_err_history = [];

    steps_taken = 0;
    rejected_steps = 0;
    
    while t_current < t_span(2) && steps_taken < max_steps
        steps_taken = steps_taken + 1;
        
        if t_current + h > t_span(2)
            h = t_span(2) - t_current;
        end
        
        [y_high, y_low, k_next, err] = embedded_solver(fun, t_current, y_current, h, k1, atol, rtol, tableau);
    
        t_history(end+1) = t_current;
        h_history(end+1) = h;
    
        if err < 1
            % STEP ACCEPTED
            is_accepted(end+1) = 1;
            err_history(end+1) = norm(y_low-y_high);
            
            t_prev = t_current;
            t_current = t_current + h;
            y_current = y_high;
            k1 = k_next;
            
            t_embedded(:, end+1) = t_current;
            y_embedded(:, end+1) = y_current;
            facmax_current = facmax;

            if ~isempty(sol)
                % Global error calculation
                global_err = norm(y_current - sol(t_current));
                global_err_history(end+1) = global_err;
                
                % Exact local error calculation
                y_true_start = sol(t_prev);
                [y_num_onestep, ~, ~, ~] = embedded_solver(fun, t_prev, y_true_start, h, [], atol, rtol, tableau);
                exact_local_err = norm(y_num_onestep - sol(t_current));
                exact_local_err_history(end+1) = exact_local_err;
            end
            
        else
            % STEP REJECTED
            is_accepted(end+1) = 0;
            rejected_steps = rejected_steps + 1;
            facmax_current = 1;
        end
        
        h = step_control(h, err, tableau.q, fac, facmin, facmax_current);
    end

    % Plot step size h evolution
    figure;
    t_acc = t_history(is_accepted == 1);
    h_acc = h_history(is_accepted == 1);
    t_rej = t_history(is_accepted == 0);
    h_rej = h_history(is_accepted == 0);
    
    semilogy(t_acc, h_acc, 'b.-', 'DisplayName', 'Accepted Steps'); hold on;
    semilogy(t_rej, h_rej, 'rx', 'DisplayName', 'Rejected Steps');
    grid on; xlabel('Time t'); ylabel('Step size h');
    title(['Step Size h Evolution for ', labels{j}]);
    legend('show');

    fprintf('Embedded Solver: %d total steps, %d accepted, %d rejected.\n\n', steps_taken, length(h_acc), rejected_steps);
    
    % Plot errors if analytical solution is available
    if ~isempty(sol)
        figure;
        semilogy(t_embedded(2:end), err_history, 'b-', 'DisplayName', 'Local Error Estimate');
        hold on;
        semilogy(t_embedded(2:end), global_err_history, 'r-', 'DisplayName', 'Global Error');
        semilogy(t_embedded(2:end), exact_local_err_history, 'g--', 'DisplayName', 'Exact Local Error');
        grid on;
        xlabel('Time t');
        ylabel('Error');
        title(['Error Analysis for ', labels{j}]);
        legend('show', 'Location', 'best');
    end

    % Plot of the solution
    figure;
    hold on; grid on;
    plot(t_embedded, y_embedded(1,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'y1 (Embedded)');
    if size(y0, 1) == 2
        plot(t_embedded, y_embedded(2,:), 'c-', 'LineWidth', 1.5, 'DisplayName', 'y2 (Embedded)');
    end
    if ~isempty(sol)
        % 1. Create a high-resolution time vector
        t_fine = linspace(t_span(1), t_span(end), 1000);
        
        % 2. Evaluate the reference solution (returns a matrix)
        y_fine = sol(t_fine); 
        
        % 3. Use plot by extracting individual rows (variables) from the matrix
        plot(t_fine, y_fine(1,:), 'r--', 'DisplayName', 'Reference Solution y_1');
        hold on;
        if length(y0)==2
            plot(t_fine, y_fine(2,:), 'b--', 'DisplayName', 'Reference Solution y_2');
        end
    end
    title(['ODE Solution: ', labels{j}]);
    xlabel('Time t');
    ylabel('y-value');
    legend('show', 'Location', 'best');
    hold off;
end
