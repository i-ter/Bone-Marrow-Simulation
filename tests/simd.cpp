#include <iostream>
#include <vector>
#include <numeric> // For std::iota and std::accumulate (though we'll write our own sum for timing)
#include <chrono>  // For timing
#include <iomanip> // For std::fixed and std::setprecision

#if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h> // For AVX/AVX2 intrinsics
#elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>  // For ARM NEON intrinsics
#endif


// To compile this program:
// On x86-64 (Intel/AMD with AVX2 support):
//   clang++ -O3 -mavx2 tests/simd.cpp -o tests/simd_test
// On ARM64 (like Apple Silicon M1/M2/M3 with NEON support):
//   clang++ -O3 tests/simd.cpp -o tests/simd_test (NEON is usually enabled by default)
// Then run with ./tests/simd_test

// Function to sum elements of a vector using a standard for loop
double sum_scalar(const std::vector<float>& arr) {
    double sum = 0.0;
    for (size_t i = 0; i < arr.size(); ++i) {
        sum += arr[i];
    }
    return sum;
}

// AVX instructions are not supported on Apple Silicon (ARM-based)
// Function to sum elements of a vector using AVX2 SIMD intrinsics
/*
double sum_simd_avx2(const std::vector<float>& arr) {
    double total_sum_double = 0.0;
    const size_t n = arr.size();
    if (n == 0) return 0.0;

    const float* data_ptr = arr.data();
    const size_t simd_width = 8; // 8 floats in a __m256 register (AVX)

    __m256 acc_vec = _mm256_setzero_ps(); // Accumulator for SIMD sums

    size_t i = 0;
    // Process chunks of 8 floats using AVX2
    for (; i + simd_width <= n; i += simd_width) {
        __m256 data_vec = _mm256_loadu_ps(data_ptr + i); // Unaligned load
        acc_vec = _mm256_add_ps(acc_vec, data_vec);
    }

    // Horizontal sum of the elements in acc_vec
    __m128 sum_128 = _mm_add_ps(_mm256_extractf128_ps(acc_vec, 0), _mm256_extractf128_ps(acc_vec, 1));
    __m128 temp_hadd = _mm_hadd_ps(sum_128, sum_128);
    temp_hadd = _mm_hadd_ps(temp_hadd, temp_hadd);
    
    float final_simd_sum_float = _mm_cvtss_f32(temp_hadd);
    total_sum_double = (double)final_simd_sum_float;

    // Add any remaining elements (tail) using scalar operations
    for (; i < n; ++i) {
        total_sum_double += data_ptr[i];
    }

    return total_sum_double;
}
*/

// Function to sum elements of a vector using ARM NEON SIMD intrinsics
double sum_simd(const std::vector<float>& arr) {
    const size_t n = arr.size();
    if (n == 0) return 0.0;

    double total_sum_double = 0.0; // Final accumulator in double precision
    const float* data_ptr = arr.data();
    const size_t simd_width = 4; // 4 floats in a float32x4_t

    // How many elements to process in SIMD single-precision before summing to double.
    // This should be a multiple of simd_width.
    // E.g., 1024 floats = 256 SIMD operations (256 * 4 = 1024).
    const size_t SIMD_FLUSH_BLOCK_SIZE = 1024; 
    
    size_t i = 0;
    float32x4_t current_block_acc_vec = vdupq_n_f32(0.0f); // Accumulator for current SIMD block (single precision)
    size_t elements_processed_in_current_block = 0;

    // Main SIMD loop processing blocks
    for (; i + simd_width <= n; i += simd_width) {
        float32x4_t data_vec = vld1q_f32(data_ptr + i);
        current_block_acc_vec = vaddq_f32(current_block_acc_vec, data_vec);
        elements_processed_in_current_block += simd_width;

        if (elements_processed_in_current_block >= SIMD_FLUSH_BLOCK_SIZE) {
            total_sum_double += (double)vaddvq_f32(current_block_acc_vec); // Sum block to double
            current_block_acc_vec = vdupq_n_f32(0.0f);                 // Reset block accumulator
            elements_processed_in_current_block = 0;
        }
    }

    // Add any sum remaining in the current_block_acc_vec (from the last, possibly partial, block)
    total_sum_double += (double)vaddvq_f32(current_block_acc_vec);
    
    // Add any remaining elements (tail) that were not processed by the SIMD loop.
    // 'i' is already at the correct position for the start of the tail.
    for (; i < n; ++i) {
        total_sum_double += data_ptr[i];
    }

    return total_sum_double;
}

int main() {
    const size_t N = 100000 * 10000;
    std::vector<float> data(N);

    // Initialize vector with some values (e.g., 0, 1, 2, ..., N-1)
    // std::iota(data.begin(), data.end(), 0.0f);
    for (size_t i = 0; i < N; ++i) {
        data[i] = (i % 2 == 0) ? 2.0f : 3.0f;
    }

    std::cout << "Processing " << N << " float elements." << std::endl;

    // Time scalar summation
    auto start_scalar = std::chrono::high_resolution_clock::now();
    double sum_s = sum_scalar(data);
    auto end_scalar = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_scalar = end_scalar - start_scalar;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Scalar sum: " << sum_s << ", Time: " << duration_scalar.count() << " ms" << std::endl;

    // Time SIMD summation
    auto start_simd = std::chrono::high_resolution_clock::now();
    double sum_v = sum_simd(data);
    auto end_simd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_simd = end_simd - start_simd;

    std::cout << "SIMD sum:   " << sum_v << ", Time: " << duration_simd.count() << " ms" << std::endl;
    std::cout << "Theoretical sum: " << N * (N - 1) / 2 << std::endl;

    if (duration_simd.count() > 0 && duration_scalar.count() > 0) {
        double speedup = duration_scalar.count() / duration_simd.count();
        std::cout << "Speedup: " << speedup << "x" << std::endl;
    }
    
    if (std::abs(sum_s - sum_v) > 1e-5 * std::abs(sum_s)) { // Check if sums are approximately equal
        std::cerr << "Warning: Sums do not match!" << std::endl;
    }

    return 0;
}
