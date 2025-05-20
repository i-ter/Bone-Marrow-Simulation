#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>


// compare running this file with
// clang++ -std=c++17 multithread_test.cpp -o test
// vs
// clang++ -std=c++17 -Xpreprocessor -fopenmp -lomp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib multithread_test.cpp -o test
// just compiling with openmp multithreading enables leads to 5-7x speadup


int main() {
    const int N = 100'000'000;
    std::vector<double> vec(N);

    auto timer = std::chrono::high_resolution_clock::now(); 
    #pragma omp parallel for

    for(int i=0; i<vec.size(); ++i) {
        vec[i] = sin(M_PI*float(i)/N);
    }
    auto end_time = std::chrono::high_resolution_clock::now();

    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - timer).count() << "ms" << '\n';

    return 0;
}