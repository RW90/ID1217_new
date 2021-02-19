// 
// 

#ifndef _REENTRANT 
#define _REENTRANT 
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

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

void *calcGrid(void *args) {

    //Till en början
    

    int iter = ((struct workerArgs *) args)->numIter;
    int gridSize = ((struct workerArgs *) args)->gridSize;
    int id = ((struct workerArgs *) args)->id;
    int numWorkers = ((struct workerArgs *) args)->numWorkers;
    double *oldGrid = ((struct workerArgs *) args)->oldGrid;
    double *newGrid = ((struct workerArgs *) args)->newGrid;

    printf("Hi! I'm thread %d, with numIter %d, gridSize %d, numWorkers %d\n", id, iter, gridSize, numWorkers);


    return;
    /*

    for(int iter = 0; iter < numIter; iter++){
        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[i][j] = (oldGrid[i - 1][j] + oldGrid[i + 1][j] + oldGrid[i][j - 1] + oldGrid[i][j + 1])*0.25;
            }
        }

        //Barrier neccessary

        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                oldGrid[i][j] = (newGrid[i - 1][j] + newGrid[i + 1][j] + newGrid[i][j - 1] + newGrid[i][j + 1])*0.25;
            }
        }

        //Barrier neccessary
    }

    */

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

    printGrid(gridSize, oldGrid, false);
    printf("GridSize: %d, NumIter: %d, maxError: %f, numWorkers: %d, Time: %g \n", gridSize, numIter, maxError, numWorkers, endTime - startTime);

    return 0;
}