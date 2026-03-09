clear; close all; clc;
addpath(pwd); % Assicura che i file nella cartella corrente siano visibili

% Carico ODEs e soluzioni analitiche
odes_examples

% Carico tabelle di Butcher
butcher_matrices

% Seleziona il metodo di default definito nel file delle matrici
selected_method_index = find(strcmp({methods.name}, default_method));
if isempty(selected_method_index)
    error('Il metodo di default "%s" non è stato trovato.', default_method);
end
tableau = methods(selected_method_index);
fprintf('--- Metodo Runge-Kutta selezionato: %s ---\n\n', tableau.name);

% SCRIPT DI TEST E VALIDAZIONE
% ----------------------------------------------------
functions = {@brusselator, f1, f2};
solutions = {[], sol1, sol2}; % Aggiunto array di soluzioni
T = [20, 20, 2];
y0s = {[1.5; 3], 1, 5};
labels = {'Brusselator', 'y''=3+cos(t)-y', 'y''=y(2+y)(2-y)'};

% Parametri comuni per il controllo del passo
atol = 1e-4;
rtol = 1e-4;
fac = 0.8;
facmin = 0.1;
facmax = 5.0;
max_steps = 5000;


% Per semplicità, usiamo un numero fisso di passi per i solutori non adattivi
n_steps_fixed = 200;

for j = 1:length(functions)
    fprintf('--- Calcolo in corso per l''ODE: %s ---\n', labels{j})
    fun = functions{j};
    sol = solutions{j};
    t_span = [0, T(j)];
    x0 = t_span(1);
    y0 = y0s{j};
    H_fixed = (t_span(2) - t_span(1)) / n_steps_fixed;

    fprintf('Esecuzione Richardson Solver...\n');
    y_richardson = zeros(2, n_steps_fixed + 1);
    y_richardson(:, 1) = y0;
    t_richardson = linspace(t_span(1), t_span(2), n_steps_fixed + 1);
    for i = 1:n_steps_fixed
        y_richardson(:, i+1) = richardson_solver(fun, t_richardson(i), y_richardson(:, i), H_fixed);
    end

    % --- Gauss-Legendre Solver (Fixed Step) ---
    fprintf('Esecuzione Gauss-Legendre Solver...\n');
    y_gauss = zeros(2, n_steps_fixed + 1);
    y_gauss(:, 1) = y0;
    t_gauss = linspace(t_span(1), t_span(2), n_steps_fixed + 1);
    for i = 1:n_steps_fixed
        y_gauss(:, i+1) = gauss_legendre_solver(fun, t_gauss(i), y_gauss(:, i), H_fixed);
    end
    % Calcolo del passo iniziale
    p_embedded = tableau.p; 
    h0 = initial_step_selector(fun, x0, y0, p_embedded, atol, rtol);
    fprintf('Passo iniziale suggerito: h = %.6f\n', h0);

    % Se non è disponibile una soluzione analitica, calcolane una di riferimento
    if isempty(sol)
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        sol_struct = ode15s(fun, t_span, y0, options);
        sol = @(t) deval(sol_struct, t);
    end

    % --- Embedded Solver (Adaptive Step) ---
    fprintf('Esecuzione Embedded Solver (adattivo)...\n');
    t_current = x0;
    y_current = y0;
    h = h0;
    k1 = [];

    % Array per i risultati
    t_embedded = t_current;
    y_embedded = y_current;
    
    % Vettori per il log
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
            % PASSO ACCETTATO
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
                % Calcolo errore globale
                global_err = norm(y_current - sol(t_current));
                global_err_history(end+1) = global_err;
                
                % Calcolo errore locale esatto
                y_true_start = sol(t_prev);
                [y_num_onestep, ~, ~, ~] = embedded_solver(fun, t_prev, y_true_start, h, [], atol, rtol, tableau);
                exact_local_err = norm(y_num_onestep - sol(t_current));
                exact_local_err_history(end+1) = exact_local_err;
            end
            
        else
            % PASSO RIFIUTATO
            is_accepted(end+1) = 0;
            rejected_steps = rejected_steps + 1;
            facmax_current = 1;
        end
        
        h = step_control(h, err, tableau.q, fac, facmin, facmax_current);
    end

    % Plot evoluzione passo h
    figure;
    t_acc = t_history(is_accepted == 1);
    h_acc = h_history(is_accepted == 1);
    t_rej = t_history(is_accepted == 0);
    h_rej = h_history(is_accepted == 0);
    
    semilogy(t_acc, h_acc, 'b.-', 'DisplayName', 'Passi Accettati'); hold on;
    semilogy(t_rej, h_rej, 'rx', 'DisplayName', 'Passi Rifiutati');
    grid on; xlabel('Tempo t'); ylabel('Ampiezza passo h');
    title(['Evoluzione del passo h per ', labels{j}]);
    legend('show');

    fprintf('Embedded Solver: %d passi totali, %d accettati, %d rifiutati.\n\n', steps_taken, length(h_acc), rejected_steps);
    
    % Plot degli errori se la soluzione analitica è disponibile
    if ~isempty(sol)
        figure;
        semilogy(t_embedded(2:end), err_history, 'b-', 'DisplayName', 'Stima Errore Locale');
        hold on;
        semilogy(t_embedded(2:end), global_err_history, 'r-', 'DisplayName', 'Errore Globale');
        semilogy(t_embedded(2:end), exact_local_err_history, 'g--', 'DisplayName', 'Errore Locale Esatto');
        grid on;
        xlabel('Tempo t');
        ylabel('Errore');
        title(['Analisi degli Errori per ', labels{j}]);
        legend('show', 'Location', 'best');
    end

    % Plot della soluzione
    figure;
    hold on; grid on;
    plot(t_embedded, y_embedded(1,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'y1 (Embedded)');
    if size(y0, 1) == 2
        plot(t_embedded, y_embedded(2,:), 'c-', 'LineWidth', 1.5, 'DisplayName', 'y2 (Embedded)');
    end
    if ~isempty(sol)
        % 1. Crea un vettore di tempi ad alta risoluzione
        t_fine = linspace(t_span(1), t_span(end), 1000);
        
        % 2. Valuta la soluzione di riferimento (restituisce una matrice)
        y_fine = sol(t_fine); 
        
        % 3. Usa plot estraendo le singole righe (variabili) della matrice
        plot(t_fine, y_fine(1,:), 'r--', 'DisplayName', 'Soluzione Riferimento y_1');
        hold on;
        if length(y0)==2
            plot(t_fine, y_fine(2,:), 'b--', 'DisplayName', 'Soluzione Riferimento y_2');
        end
    end
    title(['Soluzione ODE: ', labels{j}]);
    xlabel('Tempo t');
    ylabel('Valore y');
    legend('show', 'Location', 'best');
    hold off;
end

%% --- CODICE PER LA RIGENERAZIONE DELLA FIGURA 4.2 DI HAIRER ---
fprintf('\n--- Inizio della rigenerazione della Figura 4.2 di Hairer ---\n');

% --- Setup del problema e della soluzione di riferimento ---
fun_brusselator = @brusselator;
t_span_brusselator = [0, 20];
y0_brusselator = [1.5; 3];

fprintf('Calcolo della soluzione di riferimento ad alta precisione...\n');
options_ref = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);
sol_struct_ref = ode15s(fun_brusselator, t_span_brusselator, y0_brusselator, options_ref);
y_ref_end = deval(sol_struct_ref, t_span_brusselator(2));

% --- Parametri per i cicli ---
tols = logspace(-2, -10, 9);
n_solvers = 3;
work = cell(n_solvers, 1);
precision = cell(n_solvers, 1);
solver_names = {'Embedded Solver', 'Richardson Extrapolation', 'Gauss-Legendre (fixed step)'};
max_steps_fig = 1e6;

% --- Wrapper per il conteggio delle chiamate alla funzione ---
global nfe_counter;
f_wrapped = @(t,y) f_eval_wrapper(fun_brusselator, t, y);

% --- 1. Embedded Solver (Adattivo) ---
fprintf('Avvio calcoli per: %s\n', solver_names{1});
work{1} = zeros(size(tols));
precision{1} = zeros(size(tols));
for i = 1:length(tols)
    atol_current = tols(i);
    rtol_current = tols(i);
    nfe_counter = 0;
    
    t_current = t_span_brusselator(1);
    y_current = y0_brusselator;
    h = initial_step_selector(f_wrapped, t_current, y_current, tableau.p, atol_current, rtol_current);
    k1 = [];
    
    steps = 0;
    while t_current < t_span_brusselator(2) && steps < max_steps_fig
        if t_current + h > t_span_brusselator(2), h = t_span_brusselator(2) - t_current; end
        
        [y_high, ~, k_next, err] = embedded_solver(f_wrapped, t_current, y_current, h, k1, atol_current, rtol_current, tableau);
        
        if err < 1
            t_current = t_current + h;
            y_current = y_high;
            k1 = k_next; % FSAL
        end
        h = step_control(h, err, tableau.q, fac, facmin, facmax);
        steps = steps + 1;
    end
    work{1}(i) = nfe_counter;
    precision{1}(i) = norm(y_current - y_ref_end);
    fprintf('Tol: %.E, Work: %d, Precision: %.E\n', atol_current, work{1}(i), precision{1}(i));
end

% --- 2. Richardson Extrapolation (Adattivo) ---
fprintf('Avvio calcoli per: %s\n', solver_names{2});
work{2} = zeros(size(tols));
precision{2} = zeros(size(tols));
p_richardson = 4; % Ordine del metodo base (RK4)
q_richardson = 4; % Ordine dell'errore (p)
for i = 1:length(tols)
    atol_current = tols(i);
    rtol_current = tols(i);
    nfe_counter = 0;
    
    t_current = t_span_brusselator(1);
    y_current = y0_brusselator;
    h = initial_step_selector(fun_brusselator, t_current, y_current, p_richardson, atol_current, rtol_current);

    steps = 0;
    while t_current < t_span_brusselator(2) && steps < max_steps_fig
        if t_current + h > t_span_brusselator(2), h = t_span_brusselator(2) - t_current; end

        % Il solver richardson non è stato scritto per ritornare l'errore,
        % quindi lo calcoliamo qui basandoci sui suoi risultati.
        [y_extrap, y_coarse] = richardson_solver(fun_brusselator, t_current, y_current, h);
        nfe_counter = nfe_counter + 12; % 3 chiamate a RK4 che ha 4 stadi = 12 eval
        
        err_est = norm(y_extrap - y_coarse);
        sc = atol_current + rtol_current * max(norm(y_current), norm(y_coarse));
        err = err_est / sc;

        if err < 1
            t_current = t_current + h;
            y_current = y_extrap; % Avanziamo con la soluzione migliore
        end
        h = step_control(h, err, q_richardson, fac, facmin, facmax);
        steps = steps + 1;
    end
    work{2}(i) = nfe_counter;
    precision{2}(i) = norm(y_current - y_ref_end);
    fprintf('Tol: %.E, Work: %d, Precision: %.E\n', atol_current, work{2}(i), precision{2}(i));
end

% --- 3. Gauss-Legendre (Passo Fisso) ---
fprintf('Avvio calcoli per: %s\n', solver_names{3});
n_steps_array = round(logspace(2.5, 4.5, 9));
work{3} = zeros(size(n_steps_array));
precision{3} = zeros(size(n_steps_array));
for i = 1:length(n_steps_array)
    n_steps = n_steps_array(i);
    h = (t_span_brusselator(2) - t_span_brusselator(1)) / n_steps;
    nfe_counter = 0;
    
    y_current = y0_brusselator;
    for step = 1:n_steps
        y_current = gauss_legendre_solver(f_wrapped, 0, y_current, h);
    end
    
    work{3}(i) = nfe_counter;
    precision{3}(i) = norm(y_current - y_ref_end);
    fprintf('Steps: %d, Work: %d, Precision: %.E\n', n_steps, work{3}(i), precision{3}(i));
end


% --- Plotting del diagramma Work-Precision ---
figure;
hold on;
colors = {'r', 'b', 'k'};
markers = {'-o', '--s', ':d'};

for i = 1:n_solvers
    % Rimuoviamo eventuali punti con precisione zero che causano problemi nel log plot
    valid_points = precision{i} > 0;
    if any(valid_points)
        loglog(precision{i}(valid_points), work{i}(valid_points), markers{i}, 'Color', colors{i}, 'DisplayName', solver_names{i}, 'LineWidth', 1.5);
    end
end

title('Work-Precision Diagram (Brusselator Problem)');
xlabel('Global Error at t=20');
ylabel('Number of Function Evaluations');
legend('show', 'Location', 'SouthWest');
grid on;
set(gca, 'XDir', 'reverse'); % L'errore diminuisce verso destra
xlim([1e-11, 1e-1]);
hold off;


function dy = f_eval_wrapper(fun, t, y)
    global nfe_counter;
    dy = fun(t, y);
    nfe_counter = nfe_counter + 1;
end
