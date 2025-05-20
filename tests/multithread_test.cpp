#include <omp.h>
#include <iostream>

int main() {
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        #pragma omp critical
        std::cout << "Thread " << omp_get_thread_num() << " does i = " << i << std::endl;
    }
    return 0;
}
