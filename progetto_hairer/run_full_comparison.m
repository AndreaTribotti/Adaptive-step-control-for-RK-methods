% =========================================================================
% run_full_comparison.m
%
% Questo script esegue un confronto tra diversi metodi Runge-Kutta,
% includendo sia metodi con passo adattivo che metodi con passo fisso.
%
% L'obiettivo è generare un grafico di efficienza (work-precision diagram)
% che plotti l'errore globale rispetto al numero di valutazioni di funzione.
%
% Metodi Adattivi (basati su embedded_solver.m):
% 1. BS3(2) - Bogacki-Shampine 3(2)
% 2. ERK4(3) - Coppia custom 4(3)
% 3. DP5(4) - Dormand-Prince 5(4)
%
% Metodi a Passo Fisso:
% 4. Richardson Extrapolation (basato su RK4)
% 5. ERK3 (implementazione a passo fisso)
% 6. ERK4 (implementazione a passo fisso)
% =========================================================================

clear; close all; clc;

% --- Impostazioni del Problema e dei Solutori ---

fprintf('Setup del problema...\n');
butcher_matrices; % Carica le matrici di Butcher

% Dati del problema
fun = @brusselator;
t_span = [0, 20];
y0 = [1.5; 3];

% --- Calcolo della Soluzione di Riferimento ---
fprintf('Calcolo della soluzione di riferimento con ode45...\n');
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
sol_ref = ode45(fun, t_span, y0, opts);
y_ref_final = deval(sol_ref, t_span(2));
fprintf('Soluzione di riferimento calcolata.\n');

% Struttura per memorizzare i risultati
results = struct();

% --- Parte 1: Test dei Metodi Adattivi ---
fprintf('--- Inizio Test Metodi Adattivi ---\n');

methods_adaptive = {'BS3(2)', 'ERK4(3)', 'DP5(4)'};
tols = 10.^(-4:-1:-10);

for i = 1:length(methods_adaptive)
    method_name = methods_adaptive{i};
    fprintf('Testing metodo adattivo: %s\n', method_name);
    
    global_errors = zeros(1, length(tols));
    function_evals = zeros(1, length(tols));
    
    % Trova il tableau corretto
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


% --- Parte 2: Test dei Metodi a Passo Fisso ---
fprintf('--- Inizio Test Metodi a Passo Fisso ---\n');

% Per ERK3 e ERK4, abbiamo bisogno delle loro matrici (senza parte embedded)
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
    fprintf('Testing metodo a passo fisso: %s\n', method_name);

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


% --- Parte 3: Generazione del Grafico ---
figure;
hold on;
grid on;

method_names_plot = [methods_adaptive, methods_fixed];
colors = lines(length(method_names_plot)); % Genera colori distinti

for i = 1:length(method_names_plot)
    method_name = method_names_plot{i};
    valid_name = matlab.lang.makeValidName(method_name);
    
    res = results.(valid_name);
    valid_indices = res.errors > 0;
    
    if any(valid_indices)
        if strcmp(res.type, 'adaptive')
            % Plotta adattivi con linea e marker
            loglog(res.errors(valid_indices), res.fevals(valid_indices), ...
                   '-', 'DisplayName', method_name, 'LineWidth', 1.5, 'Color', colors(i,:));
        else
            % Plotta fissi solo con marker
            loglog(res.errors(valid_indices), res.fevals(valid_indices), ...
                   '--', 'DisplayName', [method_name, ' (fisso)'], 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', colors(i,:));
        end
    end
end

% Impostazioni del grafico
title('Confronto di Efficienza Metodi Adattivi vs. Passo Fisso');
xlabel('Errore Globale');
ylabel('Numero di Valutazioni di Funzione');
legend('show', 'Location', 'northwest');
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XDir', 'reverse');

% =========================================================================
% FUNZIONI LOCALI DEI SOLUTORI (riutilizzate da prima)
% =========================================================================

function [y_final, total_fevals] = solve_adaptive_embedded(fun, t_span, y0, rtol, atol, tableau)
    % Parametri di controllo del passo
    fac = 0.8; facmin = 0.1; facmax = 5.0; max_steps = 50000;

    % Inizializzazione
    t_start = t_span(1); t_end = t_span(2);
    t_current = t_start; y_current = y0;
    
    % Selettore del passo iniziale
    p = tableau.p;
    h = initial_step_selector(fun, t_current, y_current, p, atol, rtol);
    total_fevals = 2; % initial_step_selector fa 2 valutazioni

    k1_fsal = []; % k1 per il primo passo (vuoto)
    step_accepted = false; % Flag per tracciare se il passo precedente è stato accettato
    
    num_stages = length(tableau.c);
    is_fsal = isfield(tableau, 'fsal') && tableau.fsal;

    while t_current < t_end
        if t_current + h > t_end, h = t_end - t_current; end

        % Per un metodo FSAL, k1 è riutilizzato solo se il passo precedente è stato accettato
        if is_fsal && step_accepted
            k1_to_pass = k1_fsal;
            evals_this_step = num_stages - 1; % Si risparmia una valutazione
        else
            k1_to_pass = [];
            evals_this_step = num_stages; % Si eseguono tutte le valutazioni
        end
        
        % Esegui un passo e calcola k_next (che costa 1 valutazione extra)
        [y_high, ~, k_next, err] = embedded_solver(fun, t_current, y_current, h, k1_to_pass, atol, rtol, tableau);
        evals_this_step = evals_this_step + 1; % Aggiungi la valutazione per k_next

        if err < 1
            % Passo accettato
            t_current = t_current + h;
            y_current = y_high;
            if is_fsal
                k1_fsal = k_next; % Salva k_next per il prossimo passo
            end
            step_accepted = true;
            facmax_current = facmax;
        else
            % Passo rifiutato
            step_accepted = false;
            facmax_current = 1.0;
        end
        
        total_fevals = total_fevals + evals_this_step;
        
        % Adatta il passo per il prossimo tentativo
        if t_current < t_end
            h = step_control(h, err, tableau.q, fac, facmin, facmax_current);
        end
    end
    y_final = y_current;
end
