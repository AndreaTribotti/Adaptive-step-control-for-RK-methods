% SCRIPT DI TEST E VALIDAZIONE PER I SOLUTORI ODE
% ----------------------------------------------------
% Questo script esegue e confronta i diversi solutori ODE implementati
% sul problema del Brusselator.

clear; close all; clc;
addpath(pwd); % Assicura che i file nella cartella corrente siano visibili

%% 1. Impostazione del Problema
fun = @brusselator;
y0 = [1.5; 3];
t_span = [0, 20];
x0 = t_span(1);

% Parametri per il controllo del passo
atol = 1e-4;
rtol = 1e-4;
p_embedded = 4; % Ordine del metodo embedded (4(3))

% Calcolo del passo iniziale suggerito
h0 = initial_step_selector(fun, x0, y0, p_embedded, atol, rtol);
fprintf('Passo iniziale suggerito: h = %.6f', h0);


%% 2. Esecuzione dei Solutori
% Per semplicità, usiamo un numero fisso di passi per i solutori non adattivi
n_steps_fixed = 200;
H_fixed = (t_span(2) - t_span(1)) / n_steps_fixed;

% --- Richardson Solver (Fixed Step) ---
fprintf('Esecuzione Richardson Solver...');
y_richardson = zeros(2, n_steps_fixed + 1);
y_richardson(:, 1) = y0;
t_richardson = linspace(t_span(1), t_span(2), n_steps_fixed + 1);
for i = 1:n_steps_fixed
    y_richardson(:, i+1) = richardson_solver(fun, t_richardson(i), y_richardson(:, i), H_fixed);
end

% --- Gauss-Legendre Solver (Fixed Step) ---
fprintf('Esecuzione Gauss-Legendre Solver...');
y_gauss = zeros(2, n_steps_fixed + 1);
y_gauss(:, 1) = y0;
t_gauss = linspace(t_span(1), t_span(2), n_steps_fixed + 1);
for i = 1:n_steps_fixed
    y_gauss(:, i+1) = gauss_legendre_solver(fun, t_gauss(i), y_gauss(:, i), H_fixed);
end

% --- Embedded Solver (Adaptive Step) ---
fprintf('Esecuzione Embedded Solver (adattivo)...');
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

while t_current < t_span(2) && steps_taken < max_steps
    steps_taken = steps_taken + 1;
    if t_current + h > t_span(2)
        h = t_span(2) - t_current; % Adatta l'ultimo passo
    end
    
    % Esegui un passo e calcola l'errore
    [y_high, ~, k_next, err] = embedded_solver(fun, t_current, y_current, h, k1, atol, rtol);

    % Logica di accettazione/rifiuto del passo
    if err < 1
        % Passo accettato
        t_current = t_current + h;
        y_current = y_high;
        k1 = k_next; % FSAL
        
        % Salva il risultato
        t_embedded(:, end+1) = t_current;
        y_embedded(:, end+1) = y_current;
        
        % Calcola il nuovo passo (aumenta)
        facmax_current = facmax;
    else
        % Passo rifiutato
        rejected_steps = rejected_steps + 1;
        % Calcola il nuovo passo (diminuisce)
        % Come da istruzioni, dopo un rifiuto, facmax per il passo successivo è 1
        facmax_current = 1;
    end
    
    h = step_control(h, err, q, fac, facmin, facmax_current);
end
fprintf('Embedded Solver: %d passi accettati, %d passi rifiutati.', steps_taken - rejected_steps, rejected_steps);


%% 3. Visualizzazione dei Risultati
figure('Name', 'Confronto Solutori ODE sul Problema del Brusselator');
hold on;
grid on;

plot(t_richardson, y_richardson(1,:), 'r--', 'DisplayName', 'y1 (Richardson)');
plot(t_gauss, y_gauss(1,:), 'b-.', 'DisplayName', 'y1 (Gauss-Legendre)');
plot(t_embedded, y_embedded(1,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'y1 (Embedded Adattivo)');

plot(t_richardson, y_richardson(2,:), 'm--', 'DisplayName', 'y2 (Richardson)');
plot(t_gauss, y_gauss(2,:), 'c-.', 'DisplayName', 'y2 (Gauss-Legendre)');
plot(t_embedded, y_embedded(2,:), 'k-', 'LineWidth', 1.5, 'DisplayName', 'y2 (Embedded Adattivo)');

title('Sistema del Brusselator');
xlabel('Tempo t');
ylabel('Valore y');
legend('show', 'Location', 'best');
hold off;

%% Confronto y1 vs y2
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
