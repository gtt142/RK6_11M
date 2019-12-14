#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <unistd.h>

#define TEMP_BOT 50
#define TEMP_LEFT 200
#define TEMP_UP 100
#define TEMP_RIGHT 100
#define TEMP_ARC 100
#define LENGTH_BOT 8
#define LENGTH_RIGHT 3
#define LENGTH_LEFT 6
#define LENGTH_UP 5
#define RAD_ARC 3
#define APROX_X_NODES_NUM 100
#define ANIMATION_FRAME_DELAY 100
#define TIME 25
#define ALFA1 1
#define ALFA2 1



double deltaX;
double* tempCur;
double* tempNext; 
int x_nodes_num = APROX_X_NODES_NUM;
int y_nodes_num;
double time_cur;
int arc_center_x;
int arc_center_y;
FILE *fp;
FILE * gnuplotPipe;


// void testPrint() {
//     for (int j = y_nodes_num - 1; j >= 0 ; j--) {
//         for (int i = 0; i < x_nodes_num; i++)
//             printf("%10lf ", tempCur[i][j]);
//         printf("\n");
//     }
//     printf("\n");
//     for (int i = 0; i < x_nodes_num; i++) {
//         for (int j = 0; j < y_nodes_num; j++)
//             printf("%lf ", tempNext[i][j]);
//         printf("\n");
//     }
// }

void printNext() {
    for (int j = y_nodes_num - 1; j >= 0 ; j--) {
        for (int i = 0; i < x_nodes_num; i++)
            fprintf(fp, "%10lf ", tempNext[i*y_nodes_num+j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

void printCur() {
    for (int j = y_nodes_num - 1; j >= 0 ; j--) {
        for (int i = 0; i < x_nodes_num; i++)
            fprintf(fp, "%10lf ", tempCur[i*y_nodes_num+j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

// void printToGnuplot(double** temp) {
//     for(int i = 0; i < x_nodes_num; i++) {
//         for(int j = 0; j < y_nodes_num; j++) {
//             double xCoor = i*deltaX;
//             double y_coor = j*deltaX;
//             fprintf(gnuplotPipe, "%lf %lf %lf\n", xCoor, y_coor, temp[i][j]);
//         }
//         fprintf(gnuplotPipe, "\n");
//     }
//     fprintf(gnuplotPipe, "e\n");    
// }

void swapRes() {
    double* temp;
    temp = tempCur;
    tempCur = tempNext;
    tempNext = temp;
}

void initGranUsl(double* temp) {
    for (int i = 1; i < x_nodes_num - 1; i++) {
        temp[i*y_nodes_num+0] = TEMP_BOT;
    }
    for (int i = 1; i < y_nodes_num - 1; i++) {
        temp[0+i] = TEMP_LEFT;
    }
    for (int i = 1; i <= LENGTH_UP*(x_nodes_num - 1)/LENGTH_BOT; i++) {
        temp[i*y_nodes_num+y_nodes_num-1] = TEMP_UP;
    }
    for (int i = 1; i <= LENGTH_RIGHT*(x_nodes_num - 1)/LENGTH_BOT; i++) {
        temp[(x_nodes_num - 1)*y_nodes_num+i] = TEMP_RIGHT;
    }
    // temp[0][y_nodes_num/2] = temp[1][y_nodes_num/2]/(1 + deltaX);
}

void memAlloc(double** temp) {
    if ((*temp = (double*)calloc(x_nodes_num*y_nodes_num, sizeof(double*))) == NULL) {
        printf("Не хватает памяти\n");
        exit(1);
    }
    // for (int i = 0; i < x_nodes_num; i++) {
    //     if (((*temp)[i] = (double*)calloc(y_nodes_num, sizeof(double))) == NULL) {
    //         printf("Не хватает памяти\n");
    //         exit(1);
    //     }
    // }
}

double fArc(double x, double y) {
    double res = (x - LENGTH_UP)*(x - LENGTH_UP) + (y - LENGTH_RIGHT)*(y - LENGTH_RIGHT);
    return res;
}

double arcCoorX(double y) {
    double x = sqrt(RAD_ARC*RAD_ARC - (y - LENGTH_RIGHT)*(y - LENGTH_RIGHT)) + 5;
    return x;
}

double arcCoorY(double x) {
    double y = sqrt(RAD_ARC*RAD_ARC - (x - LENGTH_UP)*(x - LENGTH_UP)) + 3;
}

int main(int argc, char** argv) {
    float CPUtime = 0.0f;
    x_nodes_num -= (x_nodes_num - 1) % LENGTH_BOT;
    y_nodes_num = (x_nodes_num - 1) * LENGTH_LEFT/LENGTH_BOT + 1;
    if (x_nodes_num < 9) {
        printf("nodes number should be more\n");
        exit(0);
    }
    deltaX = (double)LENGTH_BOT / (x_nodes_num - 1);

    fp = fopen("result_125.txt", "w");
    if (fp==NULL)
        printf("Open failed\n");

    clock_t CPUstart = clock();
    memAlloc(&tempCur);
    memAlloc(&tempNext);

    initGranUsl(tempCur);
    initGranUsl(tempNext);

    arc_center_x = (x_nodes_num - 1) * LENGTH_UP/LENGTH_BOT;
    arc_center_y = (y_nodes_num - 1) * LENGTH_RIGHT/LENGTH_LEFT;


    printCur();
    time_cur = 0;
    int krt = (x_nodes_num-1)/LENGTH_BOT;
    double dt = deltaX/(krt*krt*10.0);
    // gnuplotPipe = popen("gnuplot -persistent", "w");
    // fprintf(gnuplotPipe, "set terminal gif animate delay %d\n", ANIMATION_FRAME_DELAY);
    // fprintf(gnuplotPipe, "set output 'animate.gif'\n");
    // fprintf(gnuplotPipe, "set pm3d map interpolate 0,0\n");
    // fprintf(gnuplotPipe, "set palette defined (0 'white', 0.01 '#002137', %d '#8b00ff', %d '#8b0000', %d 'yellow')\n", TEMP_LEFT/4, TEMP_LEFT/2, TEMP_LEFT);
    int iterCount = (int)(TIME/dt);
    // fprintf(gnuplotPipe, "do for [i=1:%d] {\n", iterCount+1);       
    // fprintf(gnuplotPipe, "splot '-' with pm3d\n");
    // fprintf(gnuplotPipe, "}\n");
    double nextFrameTime = ANIMATION_FRAME_DELAY/100;
    // printToGnuplot(tempCur);
    while (time_cur <= TIME) {
        double d2Tdx2 = 0;
        double d2Tdy2 = 0;
        for (int y_node = 1; y_node < y_nodes_num - 1; y_node++) {
            int inArc = 1;
            for (int x_node = 1; (x_node < x_nodes_num - 1) && inArc; x_node++) {
                if(x_node >= arc_center_x && y_node > arc_center_y) {
                    double x_coor = x_node * deltaX;
                    double y_coor = y_node * deltaX;
                    
                    if (fArc(x_coor + deltaX, y_coor) >= RAD_ARC*RAD_ARC) {
                        double arcX = arcCoorX(y_coor);
                        double mu = (arcX - x_coor)/deltaX;
                        d2Tdx2 = 2 * (mu*tempCur[(x_node-1)*y_nodes_num + y_node] - (mu + 1)*tempCur[x_node*y_nodes_num+y_node] + TEMP_ARC) / (mu*(mu+1)*deltaX*deltaX);
                        inArc = 0;
                    }
                    else {
                        d2Tdx2 = (tempCur[(x_node+1)*y_nodes_num+y_node] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[(x_node-1)*y_nodes_num+y_node])/(deltaX*deltaX);
                    }

                    if (fArc(x_coor, y_coor+deltaX) >= RAD_ARC*RAD_ARC) {
                        double arcY = arcCoorY(x_coor);
                        double lambda = (arcY - y_coor)/deltaX;
                        d2Tdy2 = 2 * (lambda*tempCur[x_node*y_nodes_num + y_node - 1] - (lambda+1)*tempCur[x_node*y_nodes_num + y_node] + TEMP_ARC) / (lambda*(lambda+1)*deltaX*deltaX);
                    }
                    else {
                        d2Tdy2 = (tempCur[x_node*y_nodes_num+y_node+1] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[x_node*y_nodes_num+y_node-1])/(deltaX*deltaX);
                    }
                }
                else {
                    d2Tdx2 = (tempCur[(x_node+1)*y_nodes_num + y_node] - 2*tempCur[x_node*y_nodes_num + y_node] + tempCur[(x_node-1)*y_nodes_num+y_node])/(deltaX*deltaX);                    
                    d2Tdy2 = (tempCur[x_node*y_nodes_num+y_node+1] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[x_node*y_nodes_num+y_node-1])/(deltaX*deltaX);
                }

                tempNext[x_node*y_nodes_num+y_node] = dt*(ALFA1*d2Tdx2 + ALFA2*d2Tdy2) + tempCur[x_node*y_nodes_num+y_node];
            }
        }
        // tempNext[0][y_nodes_num/2] = tempNext[1][y_nodes_num/2]/(1 + deltaX);        
        if(time_cur >= nextFrameTime) {
            printNext();            
            // printToGnuplot(tempNext);
            // printf("Time: %lf\t0 - %lf   1 - %lf\n",time_cur, tempNext[0][y_nodes_num/2], tempNext[1][y_nodes_num/2]);
            nextFrameTime += (double)ANIMATION_FRAME_DELAY/100;
        }
        // printf("%lf\n", time_cur);
        time_cur += dt;
        swapRes();
    }
    clock_t CPUstop = clock();
    CPUtime = 1000.*(CPUstop - CPUstart) / CLOCKS_PER_SEC;
    printf("CPU time : %.3f ms\n", CPUtime);

    // swapRes();
    // testPrint();

    // pclose(gnuplotPipe);
    fclose(fp);
    return 0;
}