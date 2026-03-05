# Istruzioni per lo Sviluppo Matlab - Progetto Hairer RK

Agisci come un esperto di analisi numerica. Implementa i seguenti moduli Matlab seguendo fedelmente le formule del testo di Hairer & Nörsett.
Nella cartella Project c'è un progetto simile a quello che andrai a fare. Puoi guardare le funzioni per farti un idea di quello che andrai a fare, ma NON devi assolutamente modificarle o usarle nel progetto che andrai a creare. Se tra quelle c'è una funzione che ti serve la dovrai riscrivere.

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
