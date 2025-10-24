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

/*
// single block sum over threads
__global__ void vsum(float *out, float *a, int n)
{
    float sum = 0 ;
    for( int i = 0 ; i < n ; i++ ) {
        sum += a[i] ;
    }
    __syncthreads() ;
    out[0] = sum ;
}
*/

// multi-block and multithread
__global__ void vsum(float *out, float *a, int n)
{
    extern __shared__ float sdata[] ; // shared memory
    int tid = threadIdx.x ;
    int idx = blockIdx.x*blockDim.x+tid ;
    sdata[tid] = a[idx] ; // global to local
    __syncthreads() ;
    for( int s = blockDim.x/2 ; s > 0 ; s>>=1 ) {
    	 if( tid < s ) {
	     sdata[tid] += sdata[tid+s] ;
	 }
	     __syncthreads() ;
    }
    if( tid == 0 ) out[blockIdx.x] = sdata[0] ;
}

__device__ inline void wred( volatile float *sdata , const int tid )
{
   sdata[tid] += sdata[tid+32] ;
   sdata[tid] += sdata[tid+16] ;
   sdata[tid] += sdata[tid+8] ;
   sdata[tid] += sdata[tid+4] ;
   sdata[tid] += sdata[tid+2] ;
   sdata[tid] += sdata[tid+1] ;
}

// multi-block and multithread
__global__ void vsum2(float *out, float *a, int n)
{
    extern __shared__ float sdata[] ; // shared memory
    int tid = threadIdx.x ;
    int idx = blockIdx.x*blockDim.x+tid ;
    sdata[tid] = a[idx] ; // global to local
    __syncthreads() ;
    for( int s = blockDim.x/2 ; s > 32 ; s>>=1 ) {
    	 if( tid < s ) {
	     sdata[tid] += sdata[tid+s] ;
	 }
	     __syncthreads() ;
    }
    if( tid < 32 ) wred( sdata , tid ) ;
    if( tid == 0 ) out[blockIdx.x] = sdata[0] ;
}

int main(){
    float *a ;
    float *d_a, *d_out; 

    // Allocate host memory
    cudaMallocHost( &a , sizeof(float) * N);

    // Initialize host arrays
    for(size_t i = 0; i < N; i++){
        a[i] = 1.0f;
    }
    // Allocate device memory
    cudaMalloc((void**)&d_a, sizeof(float) * N);

    gettimeofday(&t1, 0);	       

    // Transfer data from host to device memory
    cudaMemcpy(d_a, a, sizeof(float) * N, cudaMemcpyHostToDevice);

	 cudaDeviceSynchronize() ;   
    gettimeofday(&t2, 0) ;
    double time = (t2.tv_usec-t1.tv_usec)/1000.;
    printf("Time to copy:  %3.1f ms \n", time);

    for( int th = 1 ; th < 512 ; th*=2 ) {
    	 if( N%th != 0 ) continue ;
	 const int Block = N/th ;
    	 float *out ; 
	 cudaMallocHost( &out , Block*sizeof(float) ); 
         cudaMalloc( (void**)&d_out, sizeof(float) * Block);

	 gettimeofday(&t1, 0);	       
	 vsum2<<<Block,th>>>(d_out, d_a, N);
	 //vsum<<<Block,th>>>(d_out, d_a, N);
	 gettimeofday(&t2, 0) ;
	 double time = (t2.tv_usec-t1.tv_usec)/1000. ;
    	 printf("Time to execute: Nth %d  %e ms \n", th , time);

	 cudaDeviceSynchronize() ;
	 
	 gettimeofday(&t1, 0);	       
	 // Transfer data back to host memory
    	 cudaMemcpy(out, d_out, sizeof(float) * Block, cudaMemcpyDeviceToHost);
	 gettimeofday(&t2, 0) ;
	 time = (t2.tv_usec-t1.tv_usec)/1000. ;
    	 printf("Copy back  %e ms \n", time);

    	 // traditional sum outside
    	 float sum = 0 ;
    	 for( int i = 0 ; i < Block ; i++ ) {
    	      sum += out[i] ;
    	 }
    	 printf( "Sum %f N %d\n" , sum , N ) ;

         cudaFree(d_out);
    	 cudaFreeHost(out) ;
    }

    // Deallocate device memory
    cudaFree(d_a);

    // Deallocate host memory
    cudaFreeHost(a);
}
