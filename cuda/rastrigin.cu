#include <stdio.h>

const int DIM = 2500;

__global__ void Kernel(float *X, float *Y, float *Z) {
    unsigned int idx_X = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int idx_Y = threadIdx.y + blockIdx.y * blockDim.y;
    // printf("%d\n", idx_X*DIM+idx_Y);
    Z[idx_X*DIM+idx_Y] = 20.+X[idx_X]*X[idx_X]+Y[idx_Y]*Y[idx_Y]-10.*(__cosf(2.*M_PI*X[idx_X]) + __cosf(2.*M_PI*Y[idx_Y]));
}

void initialization(float min, float max, float* mem, int dim) {
    float delta = (max - min) / (dim - 1);
    for (int i = 0; i < dim; i++) {
        mem[i] = min + delta * i;
    }
}

int main() {
    cudaError cudaStatus;
    size_t mem_size = sizeof(float)*DIM;
    float *hostX, *hostY, *hostZ;
    float *devX, *devY, *devZ;
    hostX = (float*)malloc(mem_size);
    hostY = (float*)malloc(mem_size);
    hostZ = (float*)malloc(mem_size*DIM);
    cudaMalloc((void**)&devX, mem_size);
    cudaMalloc((void**)&devY, mem_size);
    cudaMalloc((void**)&devZ, mem_size*DIM);
    initialization(-5, 5, hostX, DIM);
    cudaMemcpy(devX, hostX, mem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(devY, devX, mem_size, cudaMemcpyDeviceToDevice);
    dim3 N_Block (32, 32, 1);
    dim3 N_Grid (DIM/32,DIM/32,1);
    Kernel <<< N_Grid, N_Block >>> (devX,devY,devZ);
    cudaMemcpy(hostZ, devZ, mem_size*DIM, cudaMemcpyDeviceToHost);
	
    cudaStatus = cudaGetLastError();
	if(cudaStatus != cudaSuccess) {
		printf("Last error: %s\n", cudaGetErrorString(cudaStatus));
		return 0;
	}
    
    cudaFree(devX);
    cudaFree(devY);
    cudaFree(devZ);

    for (int i = 0; i < DIM; i++) {
        for (int j = 0; j < DIM; j++) {
            printf("%lf;", hostZ[i*DIM + j]);
        }
        printf("\n");
    }

    free(hostX);
    free(hostY);
    free(hostZ);
    return 0;
}

