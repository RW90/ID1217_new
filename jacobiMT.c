// gcc -o mt jacobiMT.c -lpthread -lm
// GridSize: 200, NumIter: 100000, maxError: 0.000000, numWorkers: 4, Time: 13.184
// GridSize: 100, NumIter: 500000, maxError: 0.000000, numWorkers: 4, Time: 16.7525

#ifndef _REENTRANT 
#define _REENTRANT 
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>
#include <sched.h>

#define DEFAULT_NUM_WORKERS 4
#define DEFAULT_GRID_SIZE 22
#define DEFAULT_NUM_ITER 300
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE_NW 1
#define BOUNDARY_VALUE_SE 1

struct workerArgs {
    int id;
    int gridSize;
    int numWorkers;
    int numIter;
    double *oldGrid;
    double *newGrid;
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

void *calcGrid(void *args) {

    int row;
    int numIter = ((struct workerArgs *) args)->numIter;
    int gridSize = ((struct workerArgs *) args)->gridSize;
    int id = ((struct workerArgs *) args)->id;
    int numWorkers = ((struct workerArgs *) args)->numWorkers;
    double *oldGrid = ((struct workerArgs *) args)->oldGrid;
    double *newGrid = ((struct workerArgs *) args)->newGrid;
    int *barrierFlags = ((struct workerArgs *) args)->barrierFlags;

    for(int iter = 0; iter < numIter; iter++){
        for(int i = 1 + id; i < gridSize - 1; i += numWorkers) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[row + j] = (oldGrid[row - gridSize + j] + oldGrid[row + gridSize + j] + oldGrid[row + j - 1] + oldGrid[row + j + 1])*0.25;
            }
        }

        disBarrier(numWorkers, id, barrierFlags);

        for(int i = 1 + id; i < gridSize - 1; i += numWorkers) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                oldGrid[row + j] = (newGrid[row - gridSize + j] + newGrid[row + gridSize + j] + newGrid[row + j - 1] + newGrid[row + j + 1])*0.25;
            }
        }

        disBarrier(numWorkers, id, barrierFlags);
    }
}




int main(int argc, char *argv[]) {

    //Initiation
    int gridSize, numIter, numWorkers;
    double startTime, endTime, maxError;

    gridSize = (argc > 1) ? atoi(argv[1]) : DEFAULT_GRID_SIZE;
    numIter = (argc > 2) ? atoi(argv[2]) : DEFAULT_NUM_ITER;
    numWorkers = (argc > 3) ? atoi(argv[3]) : DEFAULT_NUM_WORKERS;

    double *oldGrid = malloc(gridSize * gridSize * sizeof(double));
    double *newGrid = malloc(gridSize * gridSize * sizeof(double));
    int *barrierFlags = calloc(numWorkers, sizeof(int));

    initGrid(gridSize, oldGrid);
    initGrid(gridSize, newGrid);

    //Init threads

    pthread_attr_t attr;
    pthread_t workers[numWorkers];
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    struct workerArgs *args = malloc(numWorkers * sizeof(struct workerArgs));
    for(int i = 0; i < numWorkers; i++) {
        args[i].id = i;
        args[i].gridSize = gridSize;
        args[i].numIter = numIter;
        args[i].numWorkers = numWorkers;
        args[i].oldGrid = oldGrid;
        args[i].newGrid = newGrid;
        args[i].barrierFlags = barrierFlags;
    }

    // Iterate solution
    startTime = read_timer();

    for (int i = 0; i < numWorkers; i++){

        // create returns 0 on success.
        if(pthread_create(&workers[i], &attr, calcGrid, (void *) &args[i])){
            printf("Thread unable to be created. Exiting program.");
            exit(0);
        }
    }
    for(int i = 0; i < numWorkers; i++){
        pthread_join(workers[i], NULL);
    }

    endTime = read_timer();

    maxError = 0.0;
    double cmpVal = 0;
    for(int i = 1; i < gridSize - 1; i++) {
        for(int j = 1; j < gridSize - 1; j++) {
            cmpVal = 1.0 - oldGrid[i * gridSize + j];
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

    printGrid(gridSize, oldGrid, true);
    printf("GridSize: %d, NumIter: %d, maxError: %f, numWorkers: %d, Time: %g \n", gridSize, numIter, maxError, numWorkers, endTime - startTime);

    return 0;
}