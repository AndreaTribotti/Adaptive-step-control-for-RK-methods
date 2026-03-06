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
atol = 1e-5;
rtol = 1e-5;
fac = 0.8;
facmin = 0.1;
facmax = 5.0;
max_steps = 10000;

for j = 1:length(functions)
    fprintf('--- Calcolo in corso per l''ODE: %s ---\n', labels{j})
    fun = functions{j};
    sol = solutions{j};
    t_span = [0, T(j)];
    x0 = t_span(1);
    y0 = y0s{j};

    % Calcolo del passo iniziale
    p_richardson = 2; % Ordine per Richardson
    h0 = initial_step_selector(fun, x0, y0, p_richardson, atol, rtol);
    fprintf('Passo iniziale suggerito: h = %.6f\n', h0);

    % Se non è disponibile una soluzione analitica, calcolane una di riferimento
    if isempty(sol)
        fprintf('Calcolo della soluzione di riferimento ad alta precisione con ode15s...\n');
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        sol_struct = ode15s(fun, t_span, y0, options);
        sol = @(t) deval(sol_struct, t);
        fprintf('Soluzione di riferimento calcolata.\n');
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
        %fplot(sol, t_span, 'r--', 'DisplayName', 'Soluzione Esatta');
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
