#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#define DEFAULT_GRID_SIZE 22
#define DEFAULT_NUM_ITER 100
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE_NW 7
#define BOUNDARY_VALUE_SE 5


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
}

int main(int argc, char *argv[]) {

    //Initiation
    int gridSize, numIter;

    gridSize = (argc > 1) ? atoi(argv[1]) : DEFAULT_GRID_SIZE;
    numIter = (argc > 2) ? atoi(argv[2]) : DEFAULT_NUM_ITER;

    double oldGrid[gridSize][gridSize];
    double newGrid[gridSize][gridSize];

    initGrid(gridSize, oldGrid);
    printGrid(gridSize, oldGrid);


    return 0;
}