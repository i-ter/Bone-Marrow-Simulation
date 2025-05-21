/*
Compile the Metal kernel:
xcrun -sdk macosx metal -c matmul.metal -o matmul.air
xcrun -sdk macosx metallib matmul.air -o default.metallib

Compile the Objective-C wrapper:
clang++ -std=c++17 gpu.cpp  matmul_mps.mm -framework Metal -framework Foundation -framework MetalPerformanceShaders -ObjC++ -o mps_matmul_test

*/

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

// this my implementation of matmul in matmul.metal
extern "C" void mps_matmul(const float* A, const float* B, float* C, int M, int N, int K);
extern "C" void mps_optimised_matmul(const float* A, const float* B, float* C, int M, int N, int K);


void cpu_matmul(const float* A, const float* B, float* C, int M, int N, int K) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < K; ++k) {
                sum += A[i*K + k] * B[k*N + j];
            }
            C[i*N + j] = sum;
        }
    }
}

int main() {
    const int M = 2048, N = 2048, K = 2048;
    std::vector<float> A(M*K), B(K*N), C_cpu(M*N), C_gpu(M*N), C_optimised_gpu(M*N);


    // Fill matrices with random data
    for (auto& v : A) v = rand() / float(RAND_MAX);
    for (auto& v : B) v = rand() / float(RAND_MAX);
    
    std::cout<<A[0]<<" "<< A[10]<<" "<< A[100]<<" "<< A[1000]<<std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    cpu_matmul(A.data(), B.data(), C_cpu.data(), M, N, K);
    auto t2 = std::chrono::high_resolution_clock::now();
    double cpu_time = std::chrono::duration<double, std::milli>(t2-t1).count();
    std::cout << "CPU Time: " << cpu_time << " ms\n";

    t1 = std::chrono::high_resolution_clock::now();
    mps_matmul(A.data(), B.data(), C_gpu.data(), M, N, K);
    t2 = std::chrono::high_resolution_clock::now();
    double gpu_time = std::chrono::duration<double, std::milli>(t2-t1).count();
    std::cout << "GPU Time: " << gpu_time << " ms\n";
    
    // MPS highly optimized defualt version
    auto t3 = std::chrono::high_resolution_clock::now();
    mps_optimised_matmul(A.data(), B.data(), C_optimised_gpu.data(), M, N, K);
    auto t4 = std::chrono::high_resolution_clock::now();
    double mps_time = std::chrono::duration<double, std::milli>(t4-t3).count();
    std::cout << "MPS Optimized GPU Time: " << mps_time << " ms\n";
    
    std::cout<<C_optimised_gpu[0]<<" "<< C_optimised_gpu[10]<<" "<< C_optimised_gpu[100]<<" "<< C_optimised_gpu[1000]<<std::endl;

    // Check error
    float max_err = 0;
    for (int i = 0; i < M*N; ++i)
        max_err = std::max(max_err, std::abs(C_cpu[i] - C_gpu[i]));
    std::cout << "Max difference: " << max_err << "\n";    
    
    // Check error
    max_err = 0;
    for (int i = 0; i < M*N; ++i)
        max_err = std::max(max_err, std::abs(C_cpu[i] - C_optimised_gpu[i]));
    std::cout << "Max difference: " << max_err << "\n";

}

