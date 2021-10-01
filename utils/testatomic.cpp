#include<atomic>
#include<stdio.h>
#include<memory.h>
#include<time.h>
#include<cilk/cilk.h>
#include<cilk/cilk_api.h>
#include<numa.h>
#include<iostream>
#include<sstream>
using namespace std;

inline static double drand(){
	return rand()*1.0f/RAND_MAX;
}

atomic<double> ** numa_atomic;
double ** numa_array;
int nnodes;
static void setWorkers(int n) {
  __cilkrts_end_cilk();
  //__cilkrts_init();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}

void init_arrays(int n) {
    nnodes = numa_max_node() + 1;
    numa_atomic = new atomic<double>*[nnodes];
    numa_array  = new double*[nnodes];
    for (int i = 0; i < nnodes; i ++) {
        void* memory = numa_alloc_onnode(sizeof(atomic<double>)*n, i);
        numa_atomic[i] = (atomic<double>*) memory;

        //memory = numa_alloc_onnode(sizeof(double)*n, i);
        numa_array[i] = (double*) numa_atomic[i];

        memset(numa_atomic[i], 0, sizeof(atomic<double>)*n);
        //memset(numa_array, 0, sizeof(double)*n);
        for (int j = 0; j < n; j ++) {
            numa_atomic[i][j] = drand();
            //numa_array[i][j] = j;
        }
    }
}
int main() {
    int n = 2200000;
    init_arrays(n);
    atomic<double>* test = new atomic<double>[n];
    memset(test, 0, sizeof(atomic<double>)*n);
    printf("size of atomic double = %d\n", sizeof(atomic<double>));
    double* test1 = (double*) test;
    // heat up cilk plus runtime
    //for (int iter = 0; iter < 2; iter ++) {
    //    cilk_for (int i = 0; i < n; i ++) {
    //        for (int node = 0; node < nnodes; node++) {
    //            test1[i] += numa_array[node][i];
    //        }
    //    }
    //}
    setWorkers(40);
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start);
    for (int iter = 0; iter < 2; iter ++) {
        cilk_for (int i = 0; i < n/2; i ++) {
            for (int node = 0; node < nnodes; node++) {
                if (numa_array[node][i]> 0.2) {
                    test1[i] += numa_array[node][i];
                }
            }
            test1[i] += i * iter;
        }
        cilk_for (int i = n/2; i < n; i ++) {
            for (int node = 0; node < nnodes; node++) {
                if (numa_array[node][i]> 0.2) {
                    test1[i] += numa_array[node][i];
                }
            }
            test1[i] += i * iter;
        }
    }
    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("cilk forcalculated in %.3f miliseconds using %d workers.\n",
                elapsed/1000000, __cilkrts_get_nworkers());
    

    clock_gettime(CLOCK_REALTIME, &start);
    for (int iter = 0; iter < 2; iter ++) {
        for (int i = 0; i < n; i ++) {
            for (int node = 0; node < nnodes; node ++) {
                test[i] = test[i] + numa_atomic[node][i].load();
            }
            test1[i] += i * iter;
        }
    }
    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("for calculated in %.3f miliseconds using %d workers.\n",
                elapsed/1000000, __cilkrts_get_nworkers());
    
    for (int i = 0; i < 10; i ++) {
        printf("%f \n", test[i].load());
    }
    return 0;
}
