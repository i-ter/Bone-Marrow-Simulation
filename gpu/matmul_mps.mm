// matmul_mps.mm
// Objective-C wrapper for Metal to expose the kernel to c++

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

extern "C" void mps_matmul(const float* A, const float* B, float* C, int M, int N, int K) {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    id<MTLCommandQueue> queue = [device newCommandQueue];
    
    // Load the Metal library (assume .metallib is built and in the bundle)
    NSError* error = nil;
    NSString* path = [[NSBundle mainBundle] pathForResource:@"default" ofType:@"metallib"];
    id<MTLLibrary> lib = [device newLibraryWithFile:path error:&error];
    id<MTLFunction> func = [lib newFunctionWithName:@"matmul"];
    id<MTLComputePipelineState> pso = [device newComputePipelineStateWithFunction:func error:&error];
    
    int nA = M * K, nB = K * N, nC = M * N;
    id<MTLBuffer> bufA = [device newBufferWithBytes:A length:sizeof(float)*nA options:0];
    id<MTLBuffer> bufB = [device newBufferWithBytes:B length:sizeof(float)*nB options:0];
    id<MTLBuffer> bufC = [device newBufferWithLength:sizeof(float)*nC options:0];

    id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    [encoder setComputePipelineState:pso];
    [encoder setBuffer:bufA offset:0 atIndex:0];
    [encoder setBuffer:bufB offset:0 atIndex:1];
    [encoder setBuffer:bufC offset:0 atIndex:2];
    [encoder setBytes:&M length:sizeof(int) atIndex:3];
    [encoder setBytes:&N length:sizeof(int) atIndex:4];
    [encoder setBytes:&K length:sizeof(int) atIndex:5];

    MTLSize gridSize = MTLSizeMake(M, N, 1);
    NSUInteger threadGroupSize = pso.maxTotalThreadsPerThreadgroup;
    NSUInteger threadsPerThreadgroup = (threadGroupSize < 16) ? threadGroupSize : 16;
    MTLSize threadgroupSize = MTLSizeMake(threadsPerThreadgroup, threadsPerThreadgroup, 1);

    [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    memcpy(C, bufC.contents, sizeof(float) * nC);
}

extern "C" void mps_optimised_matmul(const float* A, const float* B, float* C, int M, int N, int K) {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    id<MTLCommandQueue> queue = [device newCommandQueue];

    // Create buffers
    int nA = M * K, nB = K * N, nC = M * N;
    id<MTLBuffer> bufA = [device newBufferWithBytes:A length:sizeof(float)*nA options:0];
    id<MTLBuffer> bufB = [device newBufferWithBytes:B length:sizeof(float)*nB options:0];
    id<MTLBuffer> bufC = [device newBufferWithLength:sizeof(float)*nC options:0];

    // Set up MPSMatrix descriptors
    MPSMatrixDescriptor *descA = [MPSMatrixDescriptor matrixDescriptorWithRows:M columns:K rowBytes:sizeof(float)*K dataType:MPSDataTypeFloat32];
    MPSMatrixDescriptor *descB = [MPSMatrixDescriptor matrixDescriptorWithRows:K columns:N rowBytes:sizeof(float)*N dataType:MPSDataTypeFloat32];
    MPSMatrixDescriptor *descC = [MPSMatrixDescriptor matrixDescriptorWithRows:M columns:N rowBytes:sizeof(float)*N dataType:MPSDataTypeFloat32];

    MPSMatrix *matA = [[MPSMatrix alloc] initWithBuffer:bufA descriptor:descA];
    MPSMatrix *matB = [[MPSMatrix alloc] initWithBuffer:bufB descriptor:descB];
    MPSMatrix *matC = [[MPSMatrix alloc] initWithBuffer:bufC descriptor:descC];

    // Set up the MPSMatrixMultiplication object (C = A*B)
    MPSMatrixMultiplication *mm = [[MPSMatrixMultiplication alloc]
        initWithDevice:device
        transposeLeft:NO
        transposeRight:NO
        resultRows:M
        resultColumns:N
        interiorColumns:K
        alpha:1.0
        beta:0.0];

    id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
    [mm encodeToCommandBuffer:commandBuffer
                  leftMatrix:matA
                 rightMatrix:matB
                resultMatrix:matC];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    memcpy(C, bufC.contents, sizeof(float) * nC);
}