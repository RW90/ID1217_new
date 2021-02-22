// gcc -o st jacobiST.c -lpthread
// GridSize: 14, NumIter: 30000000, maxError: 0.000000, Time: 31.9954
// GridSize: 26, NumIter: 7500000, maxError: 0.000000, Time: 31.2199
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>
#include <sched.h>

#define DEFAULT_NUM_WORKERS 4
#define DEFAULT_GRID_SIZE 14
#define DEFAULT_NUM_ITER 300
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE_NW 1
#define BOUNDARY_VALUE_SE 1

struct workerArgs {
    int id;
    int *gridSizes;
    int numWorkers;
    int numIter;
    double **oldGrids;
    double **newGrids;
    int *barrierFlags;
};

double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void initGrid(int gridSize, double grid[][gridSize]) {
    for(int i = 0; i < gridSize; i++) {
        for(int j = 0; j < gridSize; j++) {
            if(i == 0 || j == 0) {
                grid[i][j] = BOUNDARY_VALUE_NW;
            } else if(i == gridSize - 1 || j == gridSize - 1) {
                grid[i][j] = BOUNDARY_VALUE_SE;
            } else {
                grid[i][j] = GRID_START_VALUE;
            }
        }
    }
}

void printGrid(int gridSize, double grid[][gridSize], bool toFile) {
    
    FILE *file;
    if(toFile){
        file = fopen("./seqOutput.txt", "w");
    }
    for(int i = 0; i < gridSize; i++) {
        for(int j = 0; j < gridSize; j++) {
            (toFile) ? fprintf(file, "%f ", grid[i][j]) : printf("%.2f ", grid[i][j]);
        }
        (toFile) ? fprintf(file, "\n") : printf("\n");
    }
    if(!toFile) {
        for(int i = 0; i < gridSize; i++) {
            printf("----");
        }
        printf("\n");
    } else {
        fclose(file);
    }
}

void disBarrier(int numWorkers, int id, int *flags) {

    int numRounds = (int) ceil(log(numWorkers) / log(2)); // log(n) / log(2) is log2
    int partner;

    for(int i = 0; i < numRounds; i++) {
        flags[id]++;
        partner = ((int)(id + pow(2, i))) % numWorkers;
        while(flags[partner] < flags[id]) {
            //sched_yield();
        }
    }
}

void calcGrid(int numIter, int gridSize, double *oldGrid, double *newGrid, int id, int numWorkers, int *barrierFlags) {

    //printf("Tråd %d hoppar in i forloop\n", id);
    int row;

    for(int iter = 0; iter < numIter; iter++){
        for(int i = 1 + id; i < gridSize - 1; i += numWorkers) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[row + j] = (oldGrid[row - gridSize + j] + oldGrid[row + gridSize + j] + oldGrid[row + j - 1] + oldGrid[row + j + 1])*0.25;
            }
        }

        //printf("Tråd %d hoppar in i första barriären\n", id);
        disBarrier(numWorkers, id, barrierFlags);

        for(int i = 1 + id; i < gridSize - 1; i += numWorkers) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                oldGrid[row + j] = (newGrid[row - gridSize + j] + newGrid[row + gridSize + j] + newGrid[row + j - 1] + newGrid[row + j + 1])*0.25;
            }
        }

        //printf("Hoppar in i andra barriären\n");
        disBarrier(numWorkers, id, barrierFlags);

        //printf("Loop %d klar\n", iter);
    }
}

void collapseGrid(int coarseGridSize, int fineGridSize, double *coarseGrid, double *fineGrid, int id, int numWorkers, int *barrierFlags) {

    int coarseRow, fineRow, fineCol;

    for(int i = 1 + id; i < coarseGridSize - 1; i += numWorkers) {
        coarseRow = i * coarseGridSize;
        fineRow = (2 * i) * fineGridSize;
        fineCol = 2;
        for(int j = 1; j < coarseGridSize - 1; j++) {
            coarseGrid[coarseRow + j] = fineGrid[fineRow + fineCol] * 0.5 + 
            (fineGrid[fineRow - fineGridSize + fineCol] + fineGrid[fineRow + fineGridSize + fineCol] + fineGrid[fineRow + fineCol - 1] + fineGrid[fineRow + fineCol + 1]) * 0.125;
            fineCol += 2;
        }
    }
    
    disBarrier(numWorkers, id, barrierFlags);

   /*
    printGrid(fineGridSize, fineGrid, false);
    printf("--------------------------------\n");
    printGrid(coarseGridSize, coarseGrid, false);
    */
}

void expandGrid(int coarseGridSize, int fineGridSize, double *coarseGrid, double *fineGrid, int id, int numWorkers, int *barrierFlags) {

    int coarseRow, fineRow, fineCol;

    //Map values from coarse grid to fine grid
    for(int i = 1 + id; i < coarseGridSize - 1; i += numWorkers) {
        coarseRow = i * coarseGridSize;
        fineRow = (2 * i) * fineGridSize;
        fineCol = 2;
        for(int j = 1; j < coarseGridSize - 1; j++) {
            fineGrid[fineRow + fineCol] = coarseGrid[coarseRow + j];
            fineCol += 2;
        }
    }

    disBarrier(numWorkers, id, barrierFlags);

    //Calculate values for the columns for points inbetween mapped points
    for(int i = 1 + id; i < fineGridSize - 1; i += 2 * numWorkers) {
        fineRow = i * fineGridSize;
        for(int j = 2; j < fineGridSize - 1; j += 2) {
            fineGrid[fineRow + j] = (fineGrid[fineRow - fineGridSize + j] + fineGrid[fineRow + fineGridSize + j]) * 0.5;
        }
    }

    disBarrier(numWorkers, id, barrierFlags);

    // Calculate values for the rest of the columns
    for(int i = 1 + id; i < fineGridSize - 1; i += numWorkers) {
        fineRow = i * fineGridSize;
        for(int j = 1; j < fineGridSize - 1; j += 2) {
            fineGrid[fineRow + j] = (fineGrid[fineRow + j - 1] + fineGrid[fineRow + j + 1]) * 0.5;
        }
    }

    disBarrier(numWorkers, id, barrierFlags);

    /*
    printf("--------------------------------\n");
    printGrid(coarseGridSize, coarseGrid, false);
    printGrid(fineGridSize, fineGrid, false);
    printf("--------------------------------\n");
    */
}

//WORKER FUNCTION

void *workerFunc(void *args) {

    int numIter = ((struct workerArgs *) args)->numIter;
    int *gridSizes = ((struct workerArgs *) args)->gridSizes;
    int id = ((struct workerArgs *) args)->id;
    int numWorkers = ((struct workerArgs *) args)->numWorkers;
    double **oldGrids = ((struct workerArgs *) args)->oldGrids;
    double **newGrids = ((struct workerArgs *) args)->newGrids;
    int *barrierFlags = ((struct workerArgs *) args)->barrierFlags;

    for(int i = 0; i < 3; i++) {
        calcGrid(4, gridSizes[i], oldGrids[i], newGrids[i], id, numWorkers, barrierFlags);
        collapseGrid(gridSizes[i + 1], gridSizes[i], oldGrids[i + 1], oldGrids[i], id, numWorkers, barrierFlags);
    }

    calcGrid(numIter, gridSizes[3], oldGrids[3], newGrids[3], id, numWorkers, barrierFlags);

    for(int i = 2; i > -1; i--) {
        expandGrid(gridSizes[i + 1], gridSizes[i], oldGrids[i + 1], oldGrids[i], id, numWorkers, barrierFlags);
        calcGrid(4, gridSizes[i], oldGrids[i], newGrids[i], id, numWorkers, barrierFlags);
    }
}

//MAIN FUNCTION

int main(int argc, char *argv[]) {

    //Initiation
    int numIter, numWorkers;
    double startTime, endTime, maxError;
    int *gridSizes = malloc(4 * sizeof(int *));
    int *barrierFlags = calloc(numWorkers, sizeof(int));

    gridSizes[3] = (argc > 1) ? atoi(argv[1]) : DEFAULT_GRID_SIZE;
    numIter = (argc > 2) ? atoi(argv[2]) : DEFAULT_NUM_ITER;
    numWorkers = (argc > 3) ? atoi(argv[3]) : DEFAULT_NUM_WORKERS;

    for(int i = 2; i > -1; i--) {
        gridSizes[i] = 2 * gridSizes[i + 1] - 1;
    }

    double **oldGrids = malloc(4 * sizeof(double **));
    double **newGrids = malloc(4 * sizeof(double **));

    for(int i = 0; i < 4; i++) {
        oldGrids[i] = (double *) malloc(gridSizes[i] * gridSizes[i] * sizeof(double));
        newGrids[i] = (double *) malloc(gridSizes[i] * gridSizes[i] * sizeof(double));
        initGrid(gridSizes[i], oldGrids[i]);
        initGrid(gridSizes[i], newGrids[i]);
    }

    //Thread init
    pthread_attr_t attr;
    pthread_t workers[numWorkers];
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    struct workerArgs *args = malloc(numWorkers * sizeof(struct workerArgs));
    for(int i = 0; i < numWorkers; i++) {
        args[i].id = i;
        args[i].gridSizes = gridSizes;
        args[i].numIter = numIter;
        args[i].numWorkers = numWorkers;
        args[i].oldGrids = oldGrids;
        args[i].newGrids = newGrids;
        args[i].barrierFlags = barrierFlags;
    }
    

    // Iterate solution
    startTime = read_timer();

    for (int i = 0; i < numWorkers; i++){

        // create returns 0 on success.
        if(pthread_create(&workers[i], &attr, workerFunc, (void *) &args[i])){
            printf("Thread unable to be created. Exiting program.");
            exit(0);
        }
    }
    for(int i = 0; i < numWorkers; i++){
        pthread_join(workers[i], NULL);
    }

    endTime = read_timer();

    // Calculate max error
    maxError = 0.0;
    double cmpVal = 0;
    for(int i = 1; i < gridSizes[0] - 1; i++) {
        for(int j = 1; j < gridSizes[0] - 1; j++) {
            cmpVal = 1.0 - oldGrids[0][i * gridSizes[0] + j];
            if(cmpVal > maxError) {
                maxError = cmpVal;
            }
        }
    }

    /* For calculating epsilon
    double maxDiff = 0;
    for(int i = 1; i < gridSize - 1; i++) {
        for(int j = 1; j < gridSize - 1; j++) {
            cmpVal = oldGrid[i][j] - newGrid[i][j];
            if(cmpVal < 0) {
                cmpVal = -cmpVal;
            }
            if(cmpVal > maxDiff) {
                maxDiff = cmpVal;
            }
        }
    }
    */

    printGrid(gridSizes[0], oldGrids[0], false);
    printf("GridSize: %d, NumIter: %d, maxError: %f, numWorkers: %d, Time: %g \n", gridSizes[3], numIter, maxError, numWorkers, endTime - startTime);

    return 0;
}