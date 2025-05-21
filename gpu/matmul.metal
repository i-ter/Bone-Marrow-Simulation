#include <metal_stdlib>
using namespace metal;

kernel void matmul(
    device const float* A         [[ buffer(0) ]],
    device const float* B         [[ buffer(1) ]],
    device float* C               [[ buffer(2) ]],
    constant uint& M              [[ buffer(3) ]],
    constant uint& N              [[ buffer(4) ]],
    constant uint& K              [[ buffer(5) ]],
    uint2 gid                     [[ thread_position_in_grid ]]
) {
    uint row = gid.x;
    uint col = gid.y;
    if (row < M && col < N) {
        float sum = 0.0f;
        for (uint i = 0; i < K; ++i) {
            sum += A[row * K + i] * B[i * N + col];
        }
        C[row * N + col] = sum;
    }
}
