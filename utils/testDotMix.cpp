#include <cilk/cilk.h>
#include <time.h>
#include <sstream>
#include <iostream>
#include <cilkpub/dotmix.h>

static void setWorkers(int n) {
  __cilkrts_end_cilk();
  //__cilkrts_init();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}
int main() {
    int n = 214896320*5;
    cilkpub::DotMix rng(n);
    int x=0;

    //-----------------------1----------------------
    setWorkers(1);
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start);

    cilk_for (int i = 0; i < n; i ++) {
        rng.get();
    }

    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("rng executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());

    //-----------------------20-----------------------
    setWorkers(20);
    clock_gettime(CLOCK_REALTIME, &start);

    cilk_for (int i = 0; i < n; i ++) {
        rng.get();
    }

    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("rng executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());


    //-----------------------40-----------------------
    setWorkers(40);
    clock_gettime(CLOCK_REALTIME, &start);

    cilk_for (int i = 0; i < n; i ++) {
        rng.get();
    }

    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("rng executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());
}
