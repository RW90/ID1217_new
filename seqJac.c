#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#define DEFAULT_GRID_SIZE 22
#define DEFAULT_NUM_ITER 300
#define GRID_START_VALUE 0
#define BOUNDARY_VALUE 1
#define MAX(x, y) ((x) > (y)) ? (x) : (y)

double read_timer(void);
void printGrid(int gridSize, double grid[][gridSize], bool toFile);
void initGridValues(int gridSize, double grid[][gridSize]);
double jacobiSolver(int gridSize, double grid[][gridSize], double newGrid[][gridSize], int loops);

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

}

double jacobiSolver(int gridSize, double grid[][gridSize], double newGrid[][gridSize], int numLoops) {
    
    for(int loop = 0; loop < numLoops; loop += 2) {
        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                newGrid[i][j] = (grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1])*0.25;
            }  
        }

        for(int i = 1; i < gridSize - 1; i++) {
            for(int j = 1; j < gridSize - 1; j++) {
                grid[i][j] = (newGrid[i-1][j] + newGrid[i+1][j] + newGrid[i][j-1] + newGrid[i][j+1])*0.25;
            }  
        }
    }

    double maxDiff = 0;
    for(int i = 1; i < gridSize - 1; i++) {
        for(int j = 1; j < gridSize - 1; j++) {
            maxDiff = MAX(maxDiff, abs(grid[i][j] - newGrid[i][j]));
        }  
    }

    return maxDiff;
}

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

void printGrid(int gridSize, double grid[][gridSize], bool toFile) {
    
    FILE *file;
    if(toFile){
        file = fopen("./seqjacout.txt", "w");
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

void initGridValues(int gridSize, double grid[][gridSize]) {
    for(int i = 0; i < gridSize; i++) {
        for(int j = o; j < gridSize; j++) {
            if(i != 0 && j != 0 && i != gridSize - 1 && j != gridSize - 1) {
                grid[i][j] = 0f;
            } else {
                grid[i][j] = 1f;
            }
        }
    }
}