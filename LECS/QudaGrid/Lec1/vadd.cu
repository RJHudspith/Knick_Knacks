#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include <sys/time.h>

struct timeval t1, t2;

#define N (32*1024*1024)
#define MAX_ERR 1e-6

// single thread and block
__global__ void vector_add(float *out, float *a, float *b, int n)
{
    for(size_t i = 0; i < n; i ++){
        out[i] = a[i] + b[i];
    }
}

// multi block
__global__ void vadd_block(float *out, float *a, float *b, int n)
{
      out[blockIdx.x] = a[blockIdx.x] + b[blockIdx.x] ;
}

// multi thread
__global__ void vadd_thread(float *out, float *a, float *b, int n)
{
      out[threadIdx.x] = a[threadIdx.x] + b[threadIdx.x] ;
}

// combined
__global__ void vadd_blth(float *out, float *a, float *b, int n)
{
     // unique id
     int i = blockIdx.x * blockDim.x + threadIdx.x;
     if( i < n ) {
        out[i] = a[i] + b[i] ;
     }
}

int main(){
    float *a, *b, *out;
    float *d_a, *d_b, *d_out; 

    // Allocate host memory
    cudaMallocHost( &a , N*sizeof(float) ) ;
    cudaMallocHost( &b , N*sizeof(float) ) ;
    cudaMallocHost( &out , N*sizeof(float) ) ;

    // Initialize host arrays
    for(size_t i = 0; i < N; i++){
        a[i] = 1.0f;
        b[i] = 2.0f;
    }

    // Allocate device memory
    cudaMalloc((void**)&d_a, sizeof(float) * N);
    cudaMalloc((void**)&d_b, sizeof(float) * N);
    cudaMalloc((void**)&d_out, sizeof(float) * N);

    gettimeofday(&t1, 0);	       

    // Transfer data from host to device memory
    cudaMemcpy(d_a, a, sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(float) * N, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize() ;
    gettimeofday(&t2, 0) ;
    double time = (t2.tv_usec-t1.tv_usec)/1000.0;
    printf("Time to copy:  %3.1f ms \n", time);

    for( int th = 1 ; th<2048 ; th*=2 ) {
        if( N%th == 1 ) continue ;
        gettimeofday(&t1, 0);
    	// Executing kernel 
    	//vector_add<<<1,1>>>(d_out, d_a, d_b, N);
    	//vadd_block<<<N,1>>>(d_out, d_a, d_b, N);
    	//vadd_thread<<<1,N>>>(d_out, d_a, d_b, N);
    	vadd_blth<<<N/th,th>>>(d_out, d_a, d_b, N);
    	gettimeofday(&t2, 0) ;
    	time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    	printf("Time to execute th = %d :  %e ms \n", th , time);
    }
    cudaDeviceSynchronize() ;
    
    gettimeofday(&t1, 0) ;

    // Transfer data back to host memory
    cudaMemcpy(out, d_out, sizeof(float) * N, cudaMemcpyDeviceToHost);

    gettimeofday(&t2, 0) ;
    time = (t2.tv_usec-t1.tv_usec)/1000.0;
    printf("Time to copy:  %3.1f ms \n", time);

    // Verification
    for(size_t i = 0; i < N; i++){
        assert(fabs(out[i] - a[i] - b[i]) < MAX_ERR);
    }
    printf("out[0] = %f\n", out[0]);
    printf("PASSED\n");

    // Deallocate device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_out);

    // Deallocate host memory
    cudaFreeHost(a) ;
    cudaFreeHost(b) ;	
    cudaFreeHost(out) ;
}
