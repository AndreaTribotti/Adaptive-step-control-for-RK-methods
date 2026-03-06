# Istruzioni per lo Sviluppo Matlab - Progetto Hairer RK

Agisci come un esperto di analisi numerica. Implementa i seguenti moduli Matlab seguendo fedelmente le formule del testo di Hairer & Nörsett.
Nella cartella Project c'è un progetto simile a quello che andrai a fare. Puoi guardare le funzioni per farti un idea di quello che andrai a fare, ma NON devi assolutamente modificarle o usarle nel progetto che andrai a creare. Se tra quelle c'è una funzione che ti serve la dovrai riscrivere. Inoltre NON considerare la cartella ERK_trib.

## 1. Modulo Controllo Passo (step_control.m)
Implementa la formula (4.13):
- Input: `err`, `h`, `q`, `fac`, `facmax`, `facmin`.
- Logica: $h_{new} = h \cdot \min(facmax, \max(facmin, fac \cdot (1/err)^{1/(q+1)}))$.
- Applica la sicurezza: se un passo viene rigettato, imposta `facmax = 1` per il tentativo successivo.

## 2. Passo Iniziale (initial_step_selector.m)
Implementa l'algoritmo di Gladwell, Shampine & Brankin:
- Segui i passi (a-f) a pagina 168-169 del documento.
- Usa la norma definita in (4.11).

## 3. Solutori Specifici
- **richardson_solver.m**: Implementa l'estrapolazione seguendo il Teorema 4.1. Calcola $y_2$ (due passi $h$) e $w$ (un passo $2h$), poi calcola il valore estrapolato $\hat{y}_2$.
- **embedded_solver.m**: Implementa il metodo 4(3) con tecnica FSAL (First Same As Last) descritta a pagina 167.
- **gauss_legendre_solver.m**: Implementa i metodi impliciti di Gauss. Usa `fsolve` per il calcolo degli stadi $k_i$.

## 4. Test Case: Brusselator
Implementa il sistema (4.15):
$y_1' = 1 + y_1^2 y_2 - 4 y_1$
$y_2' = 3 y_1 - y_1^2 y_2$
Con valori iniziali $y(0) = [1.5, 3]$.

## 5. Vincoli Tecnici
- Usa sempre la **Local Extrapolation**: avanza nel tempo usando il valore di ordine $p+1$.
- Mantieni le funzioni separate e ben documentate.
- Non aggiungere messaggi di cortesia nei commenti, sii conciso e tecnico.

## 6. primo aggiornamento progetto
Ho messo mano alla cartella progetto_hairer dopo che avevi creato tutte le funzioni. ho creato i file app1.mlapp, ERK3.m che però devi ignorare. I file nuovi che devi considerare è butcher_matrices.m e odes_examples.m. butcher_matrices.m ti permette di caricare le butcher matrices. odes_examples.m ti carica altri esempi di equazioni differenziali.  HO MODIFICATO file run_solvers.m in modo tale da avere più esempi e in modo tale che ti carichi le butcher matrices tramite il file butcher_matrices.m.
Il tuo compito adesso è aggiornare il file run_solvers.m aggiungendo quanto segue:
1 voglio che per ogni esempio venga aggiunto un grafico come l'ultimo di pagina 7 del documento Hairer_Norsett_v1_Step_Control.pdf, quindi deve contenere local error estimate, global error, exat local error.
2 vorrei aggiungere nel file butcher_matrices.m altre due butcher matrices: quella relativa al metodo 4(“5”) e quella relativa al metodo 3(2). Come scelta di default abbiamo sempre il metodo che stiamo utilizzando adesso 4(3), ma questa deve poter essere semplicemente cambiata, cambiando il nome di una variabile.
