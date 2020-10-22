#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>

#define DEFAULT_GRID_SIZE 202
#define DEFAULT_NUM_ITER 30000
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE_NW 1000
#define BOUNDARY_VALUE_SE -10

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

void printGrid(int gridSize, double grid[][gridSize]) {
    for(int i = 0; i < gridSize; i++) {
        for(int j = 0; j < gridSize; j++) {
            printf("%.2f ", grid[i][j]);
        }
        printf("\n");
    }
    for(int i = 0; i < gridSize; i++) {
        printf("----");
    }
    printf("\n");
}


int main(int argc, char *argv[]) {

    //Initiation
    int gridSize, numIter;
    double startTime, endTime, maxDiff;

    gridSize = (argc > 1) ? atoi(argv[1]) : DEFAULT_GRID_SIZE;
    numIter = (argc > 2) ? atoi(argv[2]) : DEFAULT_NUM_ITER;

    double oldGrid[gridSize][gridSize];
    double newGrid[gridSize][gridSize];

    initGrid(gridSize, oldGrid);
    initGrid(gridSize, newGrid);

    // Iterate solution
    startTime = read_timer();

    for(int iter = 0; iter < numIter; iter++){
        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[i][j] = (oldGrid[i - 1][j] + oldGrid[i + 1][j] + oldGrid[i][j - 1] + oldGrid[i][j + 1])*0.25;
            }
        }

        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                oldGrid[i][j] = (newGrid[i - 1][j] + newGrid[i + 1][j] + newGrid[i][j - 1] + newGrid[i][j + 1])*0.25;
            }
        }
    }

    double cmpVal = 0;
    maxDiff = 0;
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

    endTime = read_timer();
    //printGrid(gridSize, oldGrid);
    printf("GridSize: %d, NumIter: %d, maxDiff: %f, Time: %g \n", gridSize, numIter, maxDiff, endTime - startTime);

    return 0;
}