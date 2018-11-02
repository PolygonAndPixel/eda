/** Some unit tests for the file linear_algebra.cu
 *
 */
#include "helper/cuda_helper.cuh"
#include "helper/abbreviations.h"
#include "../src/Minimizer/MN_gpu/linear_algebra.cu"

#include "Minimizer/MN_gpu/linear_algebra.cuh"

#include <iostream>

template<
    typename value_t,
    uint block_size>
__global__ 
void test_reduction(
    value_t * A,
    uint n_dims,
    uint n,
    value_t * B)
{
    extern __shared__ value_t cache[];
    sum_reduce_2D<value_t, block_size>(A, n_dims, n, B, 1.0, cache);
}

int main()
{
    // Test reduction
    uint d = 5;
    uint n = 100;
    float * a_f = nullptr, * A_f = nullptr, * b_f = nullptr, * B_f = nullptr;
    double * a_d = nullptr, * A_d = nullptr, * b_d = nullptr, * B_d = nullptr;
    cudaMalloc(&A_f, sizeof(float)*d*n);                                CUERR 
    cudaMalloc(&B_f, sizeof(float)*d*n);                                CUERR 
    cudaMalloc(&A_d, sizeof(double)*d*n);                               CUERR 
    cudaMalloc(&B_d, sizeof(double)*d*n);                               CUERR

    cudaMallocHost(&a_f, sizeof(float)*d*n);                            CUERR
    cudaMallocHost(&b_f, sizeof(float)*d*n);                            CUERR
    cudaMallocHost(&a_d, sizeof(double)*d*n);                           CUERR
    cudaMallocHost(&b_d, sizeof(double)*d*n);                           CUERR

    // initialize arrays 
    for(uint i=0; i<d*n; i++)
    {
        a_f[i] = 1.0;
        a_d[i] = 1.0;
    }

    cudaMemcpy(A_f, a_f, sizeof(float)*n*d, H2D);                       CUERR 
    cudaMemcpy(A_d, a_d, sizeof(double)*n*d, H2D);                      CUERR 

    // Start tests with different block sizes and compare results
    auto smem_size = d*SDIV(n, 1024)*sizeof(double);
    double resultd_1024 = 0;
    test_reduction<double, 1024><<<1,1024,smem_size>>>(A_d, d, n, B_d); CUERR
    cudaMemcpy(&resultd_1024, B_d, sizeof(double), D2H);                CUERR
    
    smem_size = d*SDIV(n, 512)*sizeof(double);
    double resultd_512 = 0;
    test_reduction<double, 512><<<1,512,smem_size>>>(A_d, d, n, B_d);   CUERR
    cudaMemcpy(&resultd_512, B_d, sizeof(double), D2H);                 CUERR

    smem_size = d*SDIV(n, 256)*sizeof(double);
    double resultd_256 = 0;
    test_reduction<double, 256><<<1,256,smem_size>>>(A_d, d, n, B_d);   CUERR
    cudaMemcpy(&resultd_256, B_d, sizeof(double), D2H);                 CUERR

    smem_size = d*SDIV(n, 128)*sizeof(double);
    double resultd_128 = 0;
    test_reduction<double, 128><<<1,128,smem_size>>>(A_d, d, n, B_d);   CUERR
    cudaMemcpy(&resultd_128, B_d, sizeof(double), D2H);                 CUERR

    smem_size = d*SDIV(n, 64)*sizeof(double);
    double resultd_64 = 0;
    test_reduction<double, 64><<<1,64,smem_size>>>(A_d, d, n, B_d);     CUERR
    cudaMemcpy(&resultd_64, B_d, sizeof(double), D2H);                  CUERR

    smem_size = d*SDIV(n, 32)*sizeof(double);
    double resultd_32 = 0;
    test_reduction<double, 32><<<1,32,smem_size>>>(A_d, d, n, B_d);     CUERR
    cudaMemcpy(&resultd_32, B_d, sizeof(double), D2H);                  CUERR

    std::cout << "All values should be " << n*d << " for double results:\n";
    std::cout << "With 1024 threads: " << resultd_1024 << "\n";
    std::cout << "With 512 threads: " << resultd_512 << "\n";
    std::cout << "With 256 threads: " << resultd_256 << "\n";
    std::cout << "With 128 threads: " << resultd_128 << "\n";
    std::cout << "With 64 threads: " << resultd_64 << "\n";
    std::cout << "With 32 threads: " << resultd_32 << "\n";

    // smem_size = d*SDIV(n, 1024)*sizeof(float);
    // float resultf_1024 = 0;
    // test_reduction<float, 1024><<<1,1024,smem_size>>>(A_f, d, n, B_f);  CUERR
    
    // smem_size = d*SDIV(n, 512)*sizeof(float);
    // float resultf_512 = 0;
    // test_reduction<float, 512><<<1,512,smem_size>>>(A_f, d, n, B_f);    CUERR

    // smem_size = d*SDIV(n, 256)*sizeof(float);
    // float resultf_256 = 0;
    // test_reduction<float, 256><<<1,256,smem_size>>>(A_f, d, n, B_f);    CUERR

    // smem_size = d*SDIV(n, 128)*sizeof(float);
    // float resultf_128 = 0;
    // test_reduction<float, 128><<<1,128,smem_size>>>(A_f, d, n, B_f);    CUERR

    // smem_size = d*SDIV(n, 64)*sizeof(float);
    // float resultf_64 = 0;
    // test_reduction<float, 64><<<1,64,smem_size>>>(A_f, d, n, B_f);      CUERR

    // smem_size = d*SDIV(n, 32)*sizeof(float);
    // float resultf_32 = 0;
    // test_reduction<float, 32><<<1,32,smem_size>>>(A_f, d, n, B_f);      CUERR

    // std::cout << "All values should be " << n*d << " for float results:\n";
    // std::cout << "With 1024 threads: " << resultf_1024 << "\n";
    // std::cout << "With 512 threads: " << resultf_512 << "\n";
    // std::cout << "With 256 threads: " << resultf_256 << "\n";
    // std::cout << "With 128 threads: " << resultf_128 << "\n";
    // std::cout << "With 64 threads: " << resultf_64 << "\n";
    // std::cout << "With 32 threads: " << resultf_32 << "\n";
    return 0;
}