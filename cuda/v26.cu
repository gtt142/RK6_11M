#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#define ANIMATION_FRAME_DELAY 100
#define TIME 25
#define ALFA1 1
#define ALFA2 1

#define APROX_X_NODES_NUM 125
#define THREADS_X 32
#define THREADS_Y 32



double deltaX;
double* tempCur;
double* tempNext;
double* tempCurDev;
double* tempNextDev;  
int x_nodes_num = APROX_X_NODES_NUM;
int y_nodes_num;
// __global__ int y_nodes_numDev;
double time_cur;
int arc_center_x;
int arc_center_y;
FILE *fp;
FILE * gnuplotPipe;

__device__ double fArc(double, double);
__device__ double arcCoorX(double);
__device__ double arcCoorY(double);

__global__ void Kernel(double *current, double *next, size_t x_size, size_t y_size, int arc_center_x, int arc_center_y, double dt, double deltaX) {
    unsigned int x_node = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int y_node = threadIdx.y + blockIdx.y * blockDim.y;

    double d2Tdx2 = 0;
    double d2Tdy2 = 0;
    
    double x_coor = x_node * deltaX;
    double y_coor = y_node * deltaX;
    
    if ( x_node > 0 && y_node > 0 && (x_node < x_size-1) && (y_node < y_size-1) ) {

        if (!(x_node >= arc_center_x && y_node > arc_center_y && (fArc(x_coor, y_coor) >= RAD_ARC*RAD_ARC))) {
           

            if(x_node >= arc_center_x && y_node > arc_center_y) {
                
                if (fArc(x_coor + deltaX, y_coor) >= RAD_ARC*RAD_ARC) {
                    double arcX = arcCoorX(y_coor);
                    double mu = (arcX - x_coor)/deltaX;
                    d2Tdx2 = 2 * (mu*current[(x_node-1)*y_size + y_node] - (mu + 1)*current[x_node*y_size+y_node] + TEMP_ARC) / (mu*(mu+1)*deltaX*deltaX);
                    // inArc = 0;
                }
                else {
                    d2Tdx2 = (current[(x_node+1)*y_size+y_node] - 2*current[x_node*y_size+y_node] + current[(x_node-1)*y_size+y_node])/(deltaX*deltaX);
                }

                if (fArc(x_coor, y_coor+deltaX) >= RAD_ARC*RAD_ARC) {
                    double arcY = arcCoorY(x_coor);
                    double lambda = (arcY - y_coor)/deltaX;
                    d2Tdy2 = 2 * (lambda*current[x_node*y_size + y_node - 1] - (lambda+1)*current[x_node*y_size + y_node] + TEMP_ARC) / (lambda*(lambda+1)*deltaX*deltaX);
                }
                else {
                    d2Tdy2 = (current[x_node*y_size+y_node+1] - 2*current[x_node*y_size+y_node] + current[x_node*y_size+y_node-1])/(deltaX*deltaX);
                }
            }
            else {
                d2Tdx2 = (current[(x_node+1)*y_size + y_node] - 2*current[x_node*y_size + y_node] + current[(x_node-1)*y_size+y_node])/(deltaX*deltaX);                    
                d2Tdy2 = (current[x_node*y_size+y_node+1] - 2*current[x_node*y_size+y_node] + current[x_node*y_size+y_node-1])/(deltaX*deltaX);
            }

            next[x_node*y_size+y_node] = dt*(ALFA1*d2Tdx2 + ALFA2*d2Tdy2) + current[x_node*y_size+y_node];
            // next[idx_X*y_size+idx_Y] = 1.0;//20.+X[idx_X]*X[idx_X]+Y[idx_Y]*Y[idx_Y]-10.*(__cosf(2.*M_PI*X[idx_X]) + __cosf(2.*M_PI*Y[idx_Y]));
        }
    }
}


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

void printNext_Field() {
    for (int j = y_nodes_num - 1; j >= 0 ; j--) {
        for (int i = 0; i < x_nodes_num; i++)
            if (tempNext[i*y_nodes_num+j] == 1.0)
                fprintf(fp, "%d;%d ", i, j);
            else
                fprintf(fp, "*;*  ");
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

void swapResDev() {
    double* temp;
    temp = tempCurDev;
    tempCurDev = tempNextDev;
    tempNextDev = temp;
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
    if ((*temp = (double*)malloc(x_nodes_num*y_nodes_num * sizeof(double))) == NULL) {
        printf("Не хватает памяти\n");
        exit(1);
    }
    for (unsigned int i = 0; i < x_nodes_num*y_nodes_num; i++) {
        (*temp)[i] = 0.0;
    }
    // for (int i = 0; i < x_nodes_num; i++) {
    //     if (((*temp)[i] = (double*)calloc(y_nodes_num, sizeof(double))) == NULL) {
    //         printf("Не хватает памяти\n");
    //         exit(1);
    //     }
    // }
}

__device__ double fArc(double x, double y) {
    double res = (x - LENGTH_UP)*(x - LENGTH_UP) + (y - LENGTH_RIGHT)*(y - LENGTH_RIGHT);
    return res;
}

__device__ double arcCoorX(double y) {
    double x = sqrt(RAD_ARC*RAD_ARC - (y - LENGTH_RIGHT)*(y - LENGTH_RIGHT)) + 5;
    return x;
}

__device__ double arcCoorY(double x) {
    double y = sqrt(RAD_ARC*RAD_ARC - (x - LENGTH_UP)*(x - LENGTH_UP)) + 3;
    return y;
}

int main(int argc, char** argv) {
    cudaError cudaStatus;
	cudaEvent_t GPUstart, GPUstop;
	float GPUtime = 0.0f;
    
    x_nodes_num -= (x_nodes_num - 1) % LENGTH_BOT;
    y_nodes_num = (x_nodes_num - 1) * LENGTH_LEFT/LENGTH_BOT + 1;
    
    // cudaMemcpy(&y_nodes_numDev, &y_nodes_num, 1, cudaMemcpyHostToDevice);
    // y_nodes_numDev = y_nodes_num;
    cudaStatus = cudaGetLastError();
    if(cudaStatus != cudaSuccess) {
        printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
        return 0;
    }
    printf("y_nodes\n");

    if (x_nodes_num < 9) {
        printf("nodes number should be more\n");
        exit(0);
    }
    deltaX = (double)LENGTH_BOT / (x_nodes_num - 1);

    fp = fopen("result_cuda.txt", "w");
    if (fp==NULL)
        printf("Open failed\n");

    memAlloc(&tempCur);
    memAlloc(&tempNext);
    initGranUsl(tempCur);
    initGranUsl(tempNext);
    
	cudaEventCreate(&GPUstart);
	cudaEventCreate(&GPUstop);

	cudaEventRecord(GPUstart, 0);

    printf("=====");
    cudaMalloc((void**)&tempCurDev, x_nodes_num*y_nodes_num * sizeof(double));
    printf("Mem = %ld\n", x_nodes_num*y_nodes_num * sizeof(double));
    cudaStatus = cudaGetLastError();
    if(cudaStatus != cudaSuccess) {
        printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
        return 0;
    }
    printf("malloc 1");
    cudaMalloc((void**)&tempNextDev, x_nodes_num*y_nodes_num * sizeof(double));
    cudaStatus = cudaGetLastError();
    if(cudaStatus != cudaSuccess) {
        printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
        return 0;
    }
    printf("malloc 2");
    cudaMemcpy(tempCurDev, tempCur, x_nodes_num*y_nodes_num*sizeof(double), cudaMemcpyHostToDevice);
    cudaStatus = cudaGetLastError();
    if(cudaStatus != cudaSuccess) {
        printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
        return 0;
    }
    printf("memcpy 1");
    cudaMemcpy(tempNextDev, tempCurDev, x_nodes_num*y_nodes_num*sizeof(double), cudaMemcpyDeviceToDevice);
    cudaStatus = cudaGetLastError();
    if(cudaStatus != cudaSuccess) {
        printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
        return 0;
    }
    printf("memcpy 2\n");
    
    printf("x = %d; y = %d\n", x_nodes_num, y_nodes_num);

    int blocks_X, blocks_Y;

	if ((x_nodes_num % THREADS_X) == 0) {
		blocks_X = (x_nodes_num / THREADS_X);
	}
	else {
		blocks_X = (x_nodes_num / THREADS_X) + 1;
    }
    if ((y_nodes_num % THREADS_Y) == 0) {
		blocks_Y = (y_nodes_num / THREADS_Y);
	}
	else {
		blocks_Y = (y_nodes_num / THREADS_Y) + 1;
	}

    dim3 N_Block (THREADS_X, THREADS_Y, 1);
    dim3 N_Grid (blocks_X, blocks_Y, 1);

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
    // int iterCount = (int)(TIME/dt);
    // fprintf(gnuplotPipe, "do for [i=1:%d] {\n", iterCount+1);       
    // fprintf(gnuplotPipe, "splot '-' with pm3d\n");
    // fprintf(gnuplotPipe, "}\n");
    double nextFrameTime = ANIMATION_FRAME_DELAY/100;
    // printToGnuplot(tempCur);
    while (time_cur <= TIME) {

        Kernel <<< N_Grid, N_Block >>> (tempCurDev, tempNextDev, x_nodes_num, y_nodes_num, arc_center_x, arc_center_y, dt, deltaX); //////////////////
        cudaStatus = cudaGetLastError();
        if(cudaStatus != cudaSuccess) {
            printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
            return 0;
        }

        // tempNext[0][y_nodes_num/2] = tempNext[1][y_nodes_num/2]/(1 + deltaX);        
        if(time_cur >= nextFrameTime) {
            cudaMemcpy(tempNext, tempNextDev, x_nodes_num*y_nodes_num*sizeof(double), cudaMemcpyDeviceToHost);
            // printf("%lf\n", tempNext[0]);
            printNext();
            // printNext_Field();
            // printToGnuplot(tempNext);
            // printf("Time: %lf\t0 - %lf   1 - %lf\n",time, tempNext[0][y_nodes_num/2], tempNext[1][y_nodes_num/2]);
            nextFrameTime += (double)ANIMATION_FRAME_DELAY/100;
        }
        // printf("%lf\n", time);
        time_cur += dt;
        swapResDev();
    }

    cudaEventRecord(GPUstop, 0);
	cudaEventSynchronize(GPUstop);

	cudaEventElapsedTime(&GPUtime, GPUstart, GPUstop);
	printf("GPU time : %.3f ms\n", GPUtime);

    // swapRes();
    // testPrint();

    // pclose(gnuplotPipe);
    fclose(fp);
    return 0;
}

// for (int y_node = 1; y_node < y_nodes_num - 1; y_node++) {
//     int inArc = 1;
//     for (int x_node = 1; (x_node < x_nodes_num - 1) && inArc; x_node++) {
//         if(x_node >= arc_center_x && y_node > arc_center_y) {
//             double x_coor = x_node * deltaX;
//             double y_coor = y_node * deltaX;
            
//             if (fArc(x_coor + deltaX, y_coor) >= RAD_ARC*RAD_ARC) {
//                 double arcX = arcCoorX(y_coor);
//                 double mu = (arcX - x_coor)/deltaX;
//                 d2Tdx2 = 2 * (mu*tempCur[(x_node-1)*y_nodes_num + y_node] - (mu + 1)*tempCur[x_node*y_nodes_num+y_node] + TEMP_ARC) / (mu*(mu+1)*deltaX*deltaX);
//                 inArc = 0;
//             }
//             else {
//                 d2Tdx2 = (tempCur[(x_node+1)*y_nodes_num+y_node] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[(x_node-1)*y_nodes_num+y_node])/(deltaX*deltaX);
//             }

//             if (fArc(x_coor, y_coor+deltaX) >= RAD_ARC*RAD_ARC) {
//                 double arcY = arcCoorY(x_coor);
//                 double lambda = (arcY - y_coor)/deltaX;
//                 d2Tdy2 = 2 * (lambda*tempCur[x_node*y_nodes_num + y_node - 1] - (lambda+1)*tempCur[x_node*y_nodes_num + y_node] + TEMP_ARC) / (lambda*(lambda+1)*deltaX*deltaX);
//             }
//             else {
//                 d2Tdy2 = (tempCur[x_node*y_nodes_num+y_node+1] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[x_node*y_nodes_num+y_node-1])/(deltaX*deltaX);
//             }
//         }
//         else {
//             d2Tdx2 = (tempCur[(x_node+1)*y_nodes_num + y_node] - 2*tempCur[x_node*y_nodes_num + y_node] + tempCur[(x_node-1)*y_nodes_num+y_node])/(deltaX*deltaX);                    
//             d2Tdy2 = (tempCur[x_node*y_nodes_num+y_node+1] - 2*tempCur[x_node*y_nodes_num+y_node] + tempCur[x_node*y_nodes_num+y_node-1])/(deltaX*deltaX);
//         }

//         tempNext[x_node*y_nodes_num+y_node] = dt*(ALFA1*d2Tdx2 + ALFA2*d2Tdy2) + tempCur[x_node*y_nodes_num+y_node];
//     }
// }