// gcc -o st jacobiST.c -lpthread
// GridSize: 14, NumIter: 30000000, maxError: 0.000000, Time: 31.9954
// GridSize: 26, NumIter: 7500000, maxError: 0.000000, Time: 31.2199

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>

#define DEFAULT_GRID_SIZE 14
#define DEFAULT_NUM_ITER 300
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE_NW 1
#define BOUNDARY_VALUE_SE 1

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

void calcGrid(int numIter, int gridSize, double *oldGrid, double *newGrid) {

    int row;

    for(int iter = 0; iter < numIter; iter++){
        for(int i = 1; i < gridSize - 1; i++) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[row + j] = (oldGrid[row - gridSize + j] + oldGrid[row + gridSize + j] + oldGrid[row + j - 1] + oldGrid[row + j + 1])*0.25;
            }
        }

        for(int i = 1; i < gridSize - 1; i++) {
            row = i * gridSize;
            for(int j = 1; j < gridSize - 1; j++) {
                oldGrid[row + j] = (newGrid[row - gridSize + j] + newGrid[row + gridSize + j] + newGrid[row + j - 1] + newGrid[row + j + 1])*0.25;
            }
        }
    }
}

void collapseGrid(int coarseGridSize, int fineGridSize, double *coarseGrid, double *fineGrid) {

    int coarseRow, fineRow, fineCol;

    for(int i = 1; i < coarseGridSize - 1; i++) {
        coarseRow = i * coarseGridSize;
        fineRow = (2 * i) * fineGridSize;
        fineCol = 2;
        for(int j = 1; j < coarseGridSize - 1; j++) {
            coarseGrid[coarseRow + j] = fineGrid[fineRow + fineCol] * 0.5 + 
            (fineGrid[fineRow - fineGridSize + fineCol] + fineGrid[fineRow + fineGridSize + fineCol] + fineGrid[fineRow + fineCol - 1] + fineGrid[fineRow + fineCol + 1]) * 0.125;
            fineCol += 2;
        }
    }
}

void expandGrid(int coarseGridSize, int fineGridSize, double *coarseGrid, double *fineGrid) {

    int coarseRow, fineRow, fineCol;

    //Map values from coarse grid to fine grid
    for(int i = 1; i < coarseGridSize - 1; i++) {
        coarseRow = i * coarseGridSize;
        fineRow = (2 * i) * fineGridSize;
        fineCol = 2;
        for(int j = 1; j < coarseGridSize - 1; j++) {
            fineGrid[fineRow + fineCol] = coarseGrid[coarseRow + j];
            fineCol += 2;
        }
    }

    //Calculate values for the columns for points inbetween mapped points
    for(int i = 1; i < fineGridSize - 1; i += 2) {
        fineRow = i * fineGridSize;
        for(int j = 2; j < fineGridSize - 1; j += 2) {
            fineGrid[fineRow + j] = (fineGrid[fineRow - fineGridSize + j] + fineGrid[fineRow + fineGridSize + j]) * 0.5;
        }
    }

    // Calculate values for the rest of the columns
    for(int i = 1; i < fineGridSize - 1; i++) {
        fineRow = i * fineGridSize;
        for(int j = 1; j < fineGridSize - 1; j += 2) {
            fineGrid[fineRow + j] = (fineGrid[fineRow + j - 1] + fineGrid[fineRow + j + 1]) * 0.5;
        }
    }
}


int main(int argc, char *argv[]) {

    //Initiation
    int numIter;
    double startTime, endTime, maxError;
    int *gridSizes = malloc(4 * sizeof(int *));

    gridSizes[3] = (argc > 1) ? atoi(argv[1]) : DEFAULT_GRID_SIZE;
    numIter = (argc > 2) ? atoi(argv[2]) : DEFAULT_NUM_ITER;

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
    

    // Iterate solution
    startTime = read_timer();

    for(int i = 0; i < 3; i++) {
        calcGrid(4, gridSizes[i], oldGrids[i], newGrids[i]);
        collapseGrid(gridSizes[i + 1], gridSizes[i], oldGrids[i + 1], oldGrids[i]);
    }

    calcGrid(numIter, gridSizes[3], oldGrids[3], newGrids[3]);

    for(int i = 2; i > -1; i--) {
        expandGrid(gridSizes[i + 1], gridSizes[i], oldGrids[i + 1], oldGrids[i]);
        calcGrid(4, gridSizes[i], oldGrids[i], newGrids[i]);
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

    printGrid(gridSizes[0], oldGrids[0], true);
    printf("GridSize: %d, NumIter: %d, maxError: %f, Time: %g \n", gridSizes[3], numIter, maxError, endTime - startTime);

    return 0;
}