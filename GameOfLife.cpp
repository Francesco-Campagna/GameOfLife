#include <allegro5/allegro.h>
#include "jallegro/Allegro.h"
#include "jallegro/Frame.h"
#include "GameOfLifePanel.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <iostream>
using namespace std;
using namespace jallegro;

#define NROWS 100
#define NCOLS 100

#define v(r, c) (r * NCOLS + c) //accede agli indici di matrix (dove matrix è la matrice NROWS x NCOLS)
#define m(r, c) (r * (colsWithoutGhost+2) + c) //accede a readM e writeM (sempre locali alla sottomatrice del processo x) incluse le colonne fantasma
#define w(r,c) (r * colsWithoutGhost + c) //accede a buff (dove buff è la sottomatrice del processo x)
#define h(r,c)(r * (colsWithoutGhost + (NCOLS % partizioni_x)) + c) //accede a buff (dove buff è la sottomatrice del processo x nel caso in cui abbiamo righe/colonne extra)

int* readM; // Matrice di read
int* writeM; // Matrice di write
int* matrix; // Matrice di partenza
int* buf; // Buffer per l'invio della matrice al rank finale
int nproc; // nproc è pari a partizione_x * partizione_y + 1
int myrank; // Rank del processo
int rankLeft, rankRight, rankUp, rankDown, angleTR, angleTL, angleBR, angleBL;
int partizioni_x, partizioni_y;

int step, nthread; // Numero di step e numero di thread
int colsWithoutGhost, rowsWithoutGhost; // nrows/partizioni_y   ncols/partizioni_x

int col_start, col_end; // Variabili per il calcolo delle colonne
int row_start, row_end; // Variabili per il calcolo delle righe

double start_time, end_time; // Variabili per il calcolo del tempo

/**
 * @brief Datatype column e row, entrambi di tipo vector (vedi giù).
 */
MPI_Datatype column;
MPI_Datatype row;
MPI_Status status;

/**
 * @brief Datatype necessario per inviare a writeM.
 * 3 tipi di datatype:
 * - Contiguous: usato per costruire dataype derivato da elementi adiacenti di un array
 * - Vector: usato per costruire dataype derivato da elementi equidistanti in un array
 * - Indexed: usato per costruire dataype derivato da elementi arbitrari di un array
 * 
 * mat_shared è di tipo vector
 * 
 */
MPI_Datatype mat_shared; //vector utilizzato dall'ultimo processo per riceve le sottomatrici dai vari processi

/**
 * @brief Panel e frame per Allegro
 * 
 */
Panel* p;
Frame* f;

/**
 * @brief Struct che definisce la barriera:
 * - Mutex (lock)
 * - Cond (condizione)
 * - Count 
 * - Max_count 
 */
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count; 
    int max_count;
} Barrier;

/**
 * @brief Barriera necessaria per la sincronizzazione dei threads
 */
Barrier myBarrier;

/**
 * @brief Funzione che inizializza la barriera con il numero massimo di thread che devono aspettare su di essa.
 * 
 * @param barrier 
 * @param max_count
 *  
 * @return valore invalido/errore nella creazione del mutex/creazione variabile di condizione -1, se riuscito 0
 */
int barrier_init(Barrier *barrier, int max_count) {
    if (max_count <= 0) {
        return -1;
    }

    barrier->count = 0;
    barrier->max_count = max_count;
    if (pthread_mutex_init(&barrier->mutex, NULL) != 0) { 
        return -1; // Errore nella creazione del mutex
    }
    if (pthread_cond_init(&barrier->cond, NULL) != 0) { 
        pthread_mutex_destroy(&barrier->mutex);
        return -1; // Errore nella creazione della variabile di condizione
    }

    return 0;
}

/**
 * @brief Attende fino a quando tutti i thread hanno raggiunto la barrier.
 * Nell'else, c'è il caso in cui c'è l'ultimo thread ad essere arrivato, 
 * e di conseguenza sblocca tutti gli altri.
 * 
 * @param barrier 
 * @return int 
 */
int barrier_wait(Barrier *barrier) {
    pthread_mutex_lock(&barrier->mutex);
    barrier->count++;

    if (barrier->count < barrier->max_count) {
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
    } else {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond); //sveglio tutti i thread
    }

    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}

/**
 * @brief Distrugge la barriera e rilascia tutte le risorse precedentemente allocate.
 * 
 * @param barrier 
 * @return int 
 */
int barrier_destroy(Barrier *barrier) {
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

/**
 * @brief Funzione init: inizializza a 0 tutte le celle di writeM/readM,
 * copiando nella readM il blocco ricavato dalla matrice di partenza, escludendo righe/colonne fantasma.
 */
void init(){

    for (int i = 0; i < rowsWithoutGhost + 2; i++){  //faccio +2 così mi garantisco di copiare le celle fantasma
        for (int j = 0; j < colsWithoutGhost + 2; j++){
            readM[m(i, j)] = 0;
            writeM[m(i,j)] = 0;
        }
    }

    //parto da 1 così da non includere le celle fantasma
    for (int i = 1, ii = row_start; i <= rowsWithoutGhost; i++, ii++){
        for (int j = 1, jj = col_start; j <= colsWithoutGhost; j++, jj++){
            readM[m(i,j)] = matrix[v(ii, jj)];
            writeM[m(i,j)] = matrix[v(ii, jj)];
        }
    }
}

/**
 * @brief Funzione di swap: swappa la matrice di lettura con quella di scrittura, 
 * usando p come puntatore placeholder.
 * 
 */
void swap(){
    int *p = readM;
    readM = writeM;
    writeM = p;
}

/**
 * @brief Funzione di scambio dei bordi: 8 send, 8 recieve, BLOCCANTI.
 * 
 */
void exchBord(){
    MPI_Send(&readM[m(1,1)], 1, column, rankLeft, 0, MPI_COMM_WORLD);   //il processo x manda al rankLeft la prima colonna della sua sottomatrice
    MPI_Send(&readM[m(1, colsWithoutGhost)], 1, column, rankRight, 0, MPI_COMM_WORLD);   //il processo x manda al rankRight l'ultima colonna della sua sottomatrice
    MPI_Send(&readM[m(1,1)], 1, row, rankUp, 0, MPI_COMM_WORLD);    //il processo x manda al rankUp la prima riga della sua sottomatrice 
    MPI_Send(&readM[m(rowsWithoutGhost, 1)], 1, row, rankDown, 0, MPI_COMM_WORLD);    //il processo x manda al rankDown la sua ultima riga della sua sottomatrice 
    
    MPI_Send(&readM[m(1,colsWithoutGhost)],1, MPI_INT, angleTR, 0, MPI_COMM_WORLD);   //il processo x manda l'angolo alto destro della sua sottomatrice 
    MPI_Send(&readM[m(1,1)], 1, MPI_INT, angleTL, 0, MPI_COMM_WORLD);                 //il processo x manda l'angolo alto sinitro della sua sottomatrice
    MPI_Send(&readM[m(rowsWithoutGhost,colsWithoutGhost)], 1, MPI_INT, angleBR, 0, MPI_COMM_WORLD);   //il processo x manda l'angolo basso destro della sua sottomatrice
    MPI_Send(&readM[m(rowsWithoutGhost,1)], 1, MPI_INT, angleBL, 0, MPI_COMM_WORLD);     //il processo x manda l'angolo basso sinistro della sua sottomatrice
     
    MPI_Recv(&readM[m(1,0)], 1, column, rankLeft, 0, MPI_COMM_WORLD, &status);  // si salva la prima colonna passata da rankLeft 
    MPI_Recv(&readM[m(1, (colsWithoutGhost + 1))], 1, column, rankRight, 0, MPI_COMM_WORLD, &status);  // si salva l'ultima colonna passata da rankRight 
    MPI_Recv(&readM[m(0,1)], 1, row, rankUp, 0, MPI_COMM_WORLD, &status);    // si salva la prima riga passata da rankUp 
    MPI_Recv(&readM[m((rowsWithoutGhost + 1), 1)], 1, row, rankDown, 0, MPI_COMM_WORLD, &status);   // si salva l'ultima riga passata da rankDown 

    MPI_Recv(&readM[m(0,(colsWithoutGhost+1))], 1, MPI_INT, angleTR, 0, MPI_COMM_WORLD, &status);  // si salva l'angolo alto destro ricevuto 
    MPI_Recv(&readM[m(0,0)], 1, MPI_INT, angleTL, 0, MPI_COMM_WORLD, &status);   // si salva l'angolo alto sinistro ricevuto 
    MPI_Recv(&readM[m((rowsWithoutGhost+1), (colsWithoutGhost+1))], 1, MPI_INT, angleBR, 0, MPI_COMM_WORLD, &status);   // si salva l'angolo basso destro ricevuto 
    MPI_Recv(&readM[m((rowsWithoutGhost+1),0)], 1, MPI_INT, angleBL, 0, MPI_COMM_WORLD, &status);   // si salva l'angolo basso sinistro ricevuto 
}

/**
 * @brief Funzione che gestisce la transizione delle celle, usata dal thread del processo corrente.
 * 
 * @param i Posizione i della cella
 * @param j Posizione j della cella
 */
void transitionFunction(int i, int j){ 
    /**
     * @brief Ciclo for che controlla il numero di celle vive adiacenti, 
     * da usare per lo switch case in basso per il vicinato di Von-Neumann.
     * Escludendo la cella stessa.
     */
    int cont = 0;
    for (int di = -1; di < 2; di++){
        for (int dj = -1; dj < 2; dj++){
            //controlla che non siamo nella posizione centrale "(di != 0 || dj != 0)" E controlla che il vicinato è vivo, se si aumenta cont
            if ((di != 0 || dj != 0) && (readM[m(((i+di+rowsWithoutGhost+2) % (rowsWithoutGhost+2)), 
            ((j + dj+ colsWithoutGhost+2) % (colsWithoutGhost+2)))]) == 1){
                cont++;
            }
        }
    }

    /**
     * @brief Switch case: per ogni cella passata, 
     * devo controllare tutti i vicini utilizzando il vicinato di Von-Neumann, 2 principali casi:
     * 1. Se la cella è viva, e ci sono 2/3 celle adiacenti vive, allora resta viva, altrimenti muore.
     * 2. Se la cella è morta, e ci sono esattamente 3 celle adiacenti vive, allora viene impostata come viva, altrimenti resta morta.
     */
    switch (readM[m(i,j)]) {
        case 1:
            switch (cont) {
                case 2:
                case 3:
                    writeM[m(i,j)] = 1;
                    break;
                default:
                    writeM[m(i,j)] = 0;
                    break;
            }
            break;
        default:
            switch (cont) {
                case 3:
                    writeM[m(i,j)] = 1;
                    break;
                default:
                    writeM[m(i,j)] = 0;
                    break;
            }
            break;
    }
}

/**
 * @brief Funzione run per ogni thread:
 * Ogni thread si occupa di un numero di colonne; 
 * nel caso in cui il numero di colonne non è divisibile per il numero di thread, 
 * l'ultimo thread si occuperà anche delle colonne rimanenti aggiuntive (resto).
 * 
 * @param arg 
 * @return void* 
 * 
 */

void *run(void *arg){
    int thread = *(int*) arg;

    int startIndex = thread * ((colsWithoutGhost + 2) / nthread); //calcola la colonna di partenza
    int endIndex = startIndex + ((colsWithoutGhost + 2) / nthread); //calcola l'ultima colonna assegnata
    
    if (thread == nthread-1){
        endIndex = startIndex + ((colsWithoutGhost + 2) / nthread + (colsWithoutGhost + 2) % nthread); //se nell'ulitmo thread c'è resto dalla divisone intera, gli assegna tutte le colonne restanti
    } 

    for (int i = 1; i <= rowsWithoutGhost; i++){
        for (int j = startIndex; j <= endIndex; j++){
            if (j > 0 && j < colsWithoutGhost+1){
                transitionFunction(i,j);
            }
        }
    }
   
    barrier_wait(&myBarrier);
  
    return NULL;
}

/**
 * @brief Funzione di transizione delle celle, steps:
 * 1. Si crea il thread di tipo pthread_t (posix thread)
 * 2. Si inizializza la barriera con quel thread
 * 3. Prova a crearsi il thread
 * 4. Prova a joinare il thread
 * 5. Memorizzo la writeM escludendo righe e colonne fantasma in un buffer da inviare 
 *    all'ultimo rank che si occuperà di copiare il blocco nella propria matrice
 * 6. Distruggo la barrier appena creata
 * 7. Invio il buffer all'ultimo rank
 */
void transFunc(){
    pthread_t th[nthread];
    barrier_init(&myBarrier, nthread);
    for(int i = 0; i < nthread; i++){ 
        int *thread = new int;
        *thread = i;

        if (pthread_create(&th[i], NULL, &run, thread) != 0){
            perror("Thread create failed");
            exit(EXIT_FAILURE);
        }
    }
    
    for(int i = 0; i < nthread; i++){
        if (pthread_join(th[i], NULL) != 0){
            perror("Thread join failed");
            exit(EXIT_FAILURE);
        }
    }
    
    
    for (int i = 1, p = 0; i <= rowsWithoutGhost; i++, p++){
        for(int j = 1, t = 0; j <= colsWithoutGhost; j++, t++){
            buf[w(p,t)] = writeM[m(i,j)];
        }
    }

    barrier_destroy(&myBarrier);
    MPI_Send(&buf[w(0,0)], 1, mat_shared, (nproc-1), 0, MPI_COMM_WORLD); 
}

/**
 * @brief Funzione di calcolo delle righe e colonne di partenza e di fine blocco assegnato al processo corrente.
 * 
 * @param col_s = column start
 * @param row_s = row start
 * @param i = rank del processo
 * 
 */
void calc(int& col_s, int& row_s, int& i) {
    int diff;
    col_s = i % partizioni_x * colsWithoutGhost;
    if ((i + 1) % partizioni_x == 0){
        colsWithoutGhost +=  NCOLS % partizioni_x;
    }

    row_s = 0;
    if (i < partizioni_x){
        row_s = 0;
    }
    else{
        int count = i;
        while (count >= partizioni_x){
            row_s += rowsWithoutGhost;
            count -= partizioni_x;
        }
    }

    diff = nproc-1-i; 
    if(diff <= partizioni_x){                               //entra nell'ultimo set di sottomatrici
        rowsWithoutGhost += NROWS % partizioni_y;           //calcolo utile se la matrice non è divisibile in sottomatrici uguali
    }
}

int main(int argc, char *argv[]){
    //lettura file di configurazione contenente il numero di partizioni, numero di thread per processo e numero di step 
    FILE *fileC = fopen("configuration.txt", "r");
    if (fileC == NULL) {
        printf("Impossibile aprire il file %s\n", "configuration.txt");
        return 1;
    }  
    fscanf(fileC, "%d %d %d %d", &partizioni_x, &partizioni_y, &nthread, &step);
    fclose(fileC);
    
    //numero di righe e colonne escludendo quelle fantasma
    rowsWithoutGhost = NROWS / partizioni_y;
    colsWithoutGhost = NCOLS / partizioni_x;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Barrier(MPI_COMM_WORLD);

    // Registra il tempo di inizio
    start_time = MPI_Wtime();

    //lettura file input contenente la matrice di partenza la quale verrà copiatta nella matrix
    matrix = (int *)malloc(NROWS * NCOLS * sizeof(int)); 
    FILE *fileI = fopen("input.txt", "r");
    if (fileI == NULL) {
        printf("Impossibile aprire il file %s\n", "prova.txt");
        return 1;
    }
    
    for (int i = 0; i < NROWS; i++) {
        for (int j = 0; j < NCOLS; j++) {
            fscanf(fileI, "%d", &matrix[v(i,j)]);  // Leggi gli elementi della matrice
        }
    }

    fclose(fileI);


    //l'ultimo rank si occupa soltanto della proiezione del display, gli altri rank mandano la propria writeM che viene copiata nella matrix dell'ultimo rank
    Allegro::init();
    if (myrank == nproc - 1){
        f = new Frame(NROWS * 10, NCOLS * 10, "Game of Life");
        p = new GameOfLifePanel(matrix, NROWS, NCOLS, partizioni_x,  partizioni_y, nthread);
        f->add(p);
        f->show();

        int sizeBuffer = (colsWithoutGhost + (NCOLS % partizioni_x)) * (rowsWithoutGhost + (NROWS % partizioni_y));
        buf = (int *)malloc((sizeBuffer) * sizeof(int));
            for (int i = 0; i < sizeBuffer; i++){
               buf[i] = 0;      
            }

        
        for(int s = 0; s < step; s++){
            for(int i = 0; i < nproc - 1; i++){
                rowsWithoutGhost = NROWS / partizioni_y;
                colsWithoutGhost = NCOLS / partizioni_x;

                int col_s, row_s;
                calc(col_s, row_s, i);
               
                MPI_Type_vector(colsWithoutGhost*rowsWithoutGhost, 1, 1, MPI_INT, &mat_shared);
                MPI_Type_commit(&mat_shared); 

                MPI_Recv(&buf[0], 1, mat_shared, i, 0, MPI_COMM_WORLD, &status); 

                for (int l = 0, ii = row_s; l < rowsWithoutGhost; l++, ii++){
                    for (int g = 0, jj = col_s; g < colsWithoutGhost; g++, jj++){
                        matrix[v(ii, jj)] = buf[h(l,g)];
                    }
                }

                MPI_Type_free(&mat_shared);
            }
        
            p->repaint();
        }
        delete p;
        f->close();
        delete f;
        free(buf);
    }
    
    else{
        /*mi calcolo gli indici di riga e di colonna di partenza e di fine del blocco assegnato al processo corrente, 
        ricalcolando le righe e le colonne fantasma se necessario*/
        
        calc(col_start, row_start, myrank);
        col_end = col_start + colsWithoutGhost;
        row_end = row_start + rowsWithoutGhost;
              
        MPI_Type_vector(colsWithoutGhost*rowsWithoutGhost, 1, 1, MPI_INT, &mat_shared); //usato nell'invio-recezione del buf , tratta la sottomatrice come una sequenza lineare di interi
        MPI_Type_vector(rowsWithoutGhost, 1, colsWithoutGhost+2, MPI_INT, &column); //necessario a memorizzare le colonne delle sottomatrici escludendo le colonne fantasma
        MPI_Type_vector(colsWithoutGhost, 1, 1, MPI_INT, &row); //necessario a memorizzare le righe delle sottomatrici escludendo le righe fantasma
        MPI_Type_commit(&column);
        MPI_Type_commit(&row);
        MPI_Type_commit(&mat_shared); 

        if (nthread > colsWithoutGhost){
            nthread = colsWithoutGhost;
        }
       
        int sizeM = (rowsWithoutGhost + 2) * (colsWithoutGhost + 2);
        readM = (int*)malloc(sizeM * sizeof(int));
        writeM = (int*)malloc(sizeM * sizeof(int));

        //definiamo i vicini
        rankDown = (myrank + partizioni_x) % (nproc-1);
        rankUp = (myrank - partizioni_x + (nproc-1)) % (nproc-1);
        
        //calcolo rank left/right per le sottomatrici sul bordo destro
        if ((myrank + 1) % partizioni_x == 0){
            rankLeft = (myrank - 1 + (nproc-1)) % (nproc-1);
            rankRight = myrank + 1 - partizioni_x;
        }
        //calcolo rank left/right per le sottomatrici sul bordo sinistro
        else if (myrank % partizioni_x == 0)
        {
            rankRight = (myrank + 1) % (nproc-1);
            rankLeft = myrank + (partizioni_x - 1);
        }
        //calcolo rank left/right per le altre sottomatrici 
        else{
            rankLeft = (myrank - 1 + (nproc-1)) % (nproc-1);
            rankRight = (myrank + 1) % (nproc-1);
        }
        //calcolo per l'angolo sinistro/destro alto sul bordo destro
        if ((rankUp + 1) % partizioni_x == 0){
             angleTL = (rankUp - 1 + (nproc-1)) % (nproc-1);
             angleTR = rankUp + 1 - partizioni_x;
        }
        //calcolo per l'angolo sinistro/destro alto sul bordo sinistro
        else if (rankUp % partizioni_x == 0)
        {
            angleTR = (rankUp + 1) % (nproc-1);
            angleTL = rankUp + (partizioni_x - 1);
        }
        //calcolo per l'angolo sinistro/destro alto delle altre sottomatrici
        else{
            angleTR = (rankUp + 1) % (nproc-1);
            angleTL = (rankUp - 1 + (nproc-1)) % (nproc-1);
        }
        //calcolo per l'angolo sinistro/destro basso sul bordo destro
        if ((rankDown + 1) % partizioni_x == 0){
             angleBL = (rankDown - 1 + (nproc-1)) % (nproc-1);
             angleBR = rankDown + 1 - partizioni_x;
        }
        //calcolo per l'angolo sinistro/destro basso sul bordo sinistro
        else if (rankDown % partizioni_x == 0)
        {
            angleBR = (rankDown + 1) % (nproc-1);
            angleBL = rankDown + (partizioni_x - 1);
        }
        //calcolo per l'angolo sinistro/destro basso delle altre sottomatrici
        else{
            angleBR = (rankDown + 1) % (nproc-1);
            angleBL = (rankDown - 1 + (nproc-1)) % (nproc-1);
        }

    
        buf = (int *)malloc((colsWithoutGhost * rowsWithoutGhost) * sizeof(int));
        for (int i = 0; i < (rowsWithoutGhost); i++){
            for(int j = 0; j < (colsWithoutGhost); j++){
                buf[w(i,j)] = 0;
            }       
        }

        init();
    
        for (int s = 0; s < step; s++){           
            exchBord();
            transFunc();
            swap(); 
        }       
        
        free(readM);
        free(writeM);
        free(buf);

        MPI_Type_free(&row);
        MPI_Type_free(&column);
        MPI_Type_free(&mat_shared);
              
    }

    free(matrix);

    // Registra il tempo di fine
    end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("MPI_Finalize() encountered an error.\n");
    } else {
        if(myrank == nproc - 1) {
            printf("MPI_Finalize() successful.\n");
            printf("Tempo di esecuzione con %d processi: %.6f secondi.\n", nproc, elapsed_time);
        }
    }

    exit(0);
}