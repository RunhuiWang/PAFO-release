#include <algorithm>
#include <atomic>
#include <memory.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <ctime>
#include <unistd.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

using namespace std;

int main() {
    
    int n = 260000000;
    int k = n/40;
    double* A = new double[n];

    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start);
    
    memset(A, 0, sizeof(double)*n);

    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("sequential memset took %f ms\n", elapsed/1000000);


    clock_gettime(CLOCK_REALTIME, &start);
    
    #pragma cilk grainsize=1
    cilk_for( int i = 0; i < 40; i ++) {
        memset(A+k*i, 0, sizeof(double)*k);    
    }

    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("parallel memset took %f ms\n", elapsed/1000000);
}
