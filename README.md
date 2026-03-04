# Analisi e Implementazione di Metodi Runge-Kutta Adattivi

Questo progetto Matlab implementa le tecniche di controllo del passo e stima dell'errore descritte nel testo di Hairer & Nörsett, "Solving Ordinary Differential Equations I". Il focus è sulla precisione computazionale e sull'efficienza del solutore attraverso l'adattività del passo.

## Funzionalità Implementate

### 1. Tecniche di Stima dell'Errore
* **Richardson Extrapolation**: Confronto tra due passi di ampiezza $h$ e un passo grande $2h$ per stimare l'errore e aumentare l'ordine di precisione.
* **Embedded Runge-Kutta**: Implementazione di coppie di metodi (es. 3(2) e 4(3) FSAL) che condividono i valori della funzione per ridurre il costo computazionale.
* **Metodi Impliciti di Gauss-Legendre**: Implementazione di solutori impliciti (punto medio e ordine 4) per una stabilità superiore su problemi rigidi (stiff).

### 2. Controllo Adattivo del Passo
* **Algoritmo di Gladwell**: Calcolo automatico del passo iniziale $h$ basato sulla norma della soluzione e delle sue derivate.
* **Logica di Accettazione/Rifiuto**: Regolazione del passo basata sulle tolleranze `Atol` e `Rtol` fornite dall'utente.
* **Local Extrapolation**: Il solutore avanza sempre utilizzando l'approssimazione di ordine superiore per garantire la massima precisione.

### 3. Interfaccia Grafica (App Designer)
Un'applicazione interattiva che permette di:
* Selezionare il sistema ODE (es. Brusselator).
* Configurare i parametri di tolleranza.
* Visualizzare graficamente i passi accettati, quelli rigettati e l'andamento dell'errore.

## Requisiti
* MATLAB (R2021a o superiore)
* Control System Toolbox (opzionale per analisi di stabilità)

## Esercizi Risolti
* **Esercizio 1**: Interpretazione del metodo di Runge come estrapolazione di Richardson.
* **Esercizio 3**: Verifica dell'invarianza della strategia di controllo rispetto al riscalamento della variabile indipendente.
