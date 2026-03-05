clear; close all; clc;
addpath(pwd); % Assicura che i file nella cartella corrente siano visibili

% Carico due ODEs
odes_examples

% Carico tablue
butcher_matrices

% SCRIPT DI TEST E VALIDAZIONE PER I SOLUTORI ODE
% ----------------------------------------------------
% Questo script esegue e confronta i diversi solutori ODE implementati
% su tre differenti ODEs.

%% 1. Impostazione del Problema

functions = {@brusselator, f1, f2};
% y0_2=5; tspan2=[0 2];  per f2;
% y0=1;   tspan=[0 20]; per f1
T=[20,20,2];
y0s = {[1.5; 3],1,5};

% Parametri per il controllo del passo
atol = 1e-4;
rtol = 1e-4;
p_embedded = 4; % Ordine del metodo embedded (4(3))


labels={'brusselator',func2str(f1),func2str(f2)};
for j=1:3
    fprintf('--- Calcolo in corso per la funzione: %s ---\n', labels{j});
    fun=functions{j};
    t_span = [0, T(j)];
    x0 = t_span(1);
    y0=y0s{j};
    % Calcolo del passo iniziale suggerito
    h0 = initial_step_selector(fun, x0, y0, p_embedded, atol, rtol);
    fprintf('Passo iniziale suggerito: h = %.6f\n', h0);
    
    
    %% 2. Esecuzione dei Solutori
    % Per semplicità, usiamo un numero fisso di passi per i solutori non adattivi
    n_steps_fixed = 200;
    H_fixed = (t_span(2) - t_span(1)) / n_steps_fixed;
    
    % --- Richardson Solver (Fixed Step) ---
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
    
    % --- Embedded Solver (Adaptive Step) ---
    fprintf('Esecuzione Embedded Solver (adattivo)...\n');
    t_current = x0;
    y_current = y0;
    h = h0; % Usa il passo iniziale calcolato
    k1 = []; % Per FSAL, inizialmente vuoto
    
    % Array per salvare i risultati
    t_embedded = t_current;
    y_embedded = y_current;
    
    % Parametri per il controllo del passo
    q = 3; % Ordine dell'errore per il metodo 4(3)
    fac = 0.8;
    facmin = 0.1;
    facmax = 5.0;
    max_steps = 5000;
    steps_taken = 0;
    rejected_steps = 0;
    
    % Inizializza vettori di log prima del loop
    h_history = [];
    t_history = [];
    is_accepted = []; % 1 per accettato, 0 per rifiutato
    
    while t_current < t_span(2) && steps_taken < max_steps
        steps_taken = steps_taken + 1;
        
        if t_current + h > t_span(2)
            h_tentativo = t_span(2) - t_current; 
        else
            h_tentativo = h;
        end
        
        [y_high, ~, k_next, err] = embedded_solver(fun, t_current, y_current, h_tentativo, k1, atol, rtol);
    
        % --- LOG DI OGNI TENTATIVO ---
        t_history(end+1) = t_current;
        h_history(end+1) = h_tentativo;
    
        if err < 1
            % PASSO ACCETTATO
            is_accepted(end+1) = 1;
            
            t_current = t_current + h_tentativo;
            y_current = y_high;
            k1 = k_next;
            
            t_embedded(:, end+1) = t_current;
            y_embedded(:, end+1) = y_current;
            facmax_current = facmax;
        else
            % PASSO RIFIUTATO
            is_accepted(end+1) = 0;
            
            rejected_steps = rejected_steps + 1;
            facmax_current = 1;
        end
        
        h = step_control(h_tentativo, err, q, fac, facmin, facmax_current);
    end


    figure;

    % Separiamo i dati per il plot
    t_acc = t_history(is_accepted == 1);
    h_acc = h_history(is_accepted == 1);
    
    t_rej = t_history(is_accepted == 0);
    h_rej = h_history(is_accepted == 0);
    
    % Plot in scala semilogaritmica
    semilogy(t_acc, h_acc, 'b.-', 'LineWidth', 1, 'MarkerSize', 10); 
    hold on;
    semilogy(t_rej, h_rej, 'rx', 'LineWidth', 1.5, 'MarkerSize', 8);
    
    grid on;
    xlabel('Tempo t');
    ylabel('Ampiezza passo h');
    title(['Evoluzione del passo h per ',labels{j}]);
    legend('Step Accettati', 'Step Rifiutati', 'Location', 'best');
    
    % Aggiunta estetica: limite superiore/inferiore
    y_lims = ylim;
    ylim([y_lims(1)*0.5, y_lims(2)*1.5]);

    fprintf('Embedded Solver: %d passi accettati, %d passi rifiutati.\n\n', length(h_acc), length(h_rej));
    
    
    %% 3. Visualizzazione dei Risultati
    figure;
    hold on;
    grid on;
    
    plot(t_richardson, y_richardson(1,:), 'r--', 'DisplayName', 'y1 (Richardson)');
    plot(t_gauss, y_gauss(1,:), 'b-.', 'DisplayName', 'y1 (Gauss-Legendre)');
    plot(t_embedded, y_embedded(1,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'y1 (Embedded Adattivo)');
    
    if length(y0)==2
        plot(t_richardson, y_richardson(2,:), 'm--', 'DisplayName', 'y2 (Richardson)');
        plot(t_gauss, y_gauss(2,:), 'c-.', 'DisplayName', 'y2 (Gauss-Legendre)');
        plot(t_embedded, y_embedded(2,:), 'k-', 'LineWidth', 1.5, 'DisplayName', 'y2 (Embedded Adattivo)');
    end
    
    title(['ODE: ', labels{j} ]);
    xlabel('Tempo t');
    ylabel('Valore y');
    legend('show', 'Location', 'best');
    hold off;
    
    %% Confronto y1 vs y2
    if length(y0)==2
        figure('Name', 'Spazio delle Fasi del Brusselator');
        hold on;
        grid on;
        
        plot(y_richardson(1,:), y_richardson(2,:), 'r--', 'DisplayName', 'Richardson');
        plot(y_gauss(1,:), y_gauss(2,:), 'b-.', 'DisplayName', 'Gauss-Legendre');
        plot(y_embedded(1,:), y_embedded(2,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Embedded Adattivo');
        
        title('Spazio delle Fasi (y1 vs y2)');
        xlabel('y1');
        ylabel('y2');
        legend('show', 'Location', 'best');
        hold off;
    end
    
    % ERK3
    [t,u,num_fval] = ERK3(fun,t_span,y0,H_fixed,b3H,c3H,A3H);
end

