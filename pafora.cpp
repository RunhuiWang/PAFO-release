/*
    hybrid way ( bag and cilk_for scan)
    separate scan
    (No gorder, no set cover: commented)
*/
//#ifndef TBB_USE_GLIBCXX_VERSION
//#define TBB_USE_GLIBCXX_VERSION 11
//#endif
//#define FP_RESIDUE_DISTRIBUTION 1 // use different cores to sample 
                                  // different sources
                                  // each core maintain its residue and reserve
                                  // Will not run random walk when use this

#define TOPK_QUERY

#define _GLIBCXX_USE_CXX11_ABI 1
#define THRESHOLD 256    // bag threshold
//#define PARALLEL_INIT 1
#define SEPARATE 1
#ifndef FP_RESIDUE_DISTRIBUTION // only use bag when checking hot nodes
#define HYBRID 1
#define SOURCEREDUCE 1            // do not use reducer in multi-thread multi sources
#define RANDOMWALK 1
#define CACHE_SCHEDULE 1
#endif  // FP_RESIDUE_DISTRIBUTION
//#define MULTICORE 1           
#define RWINDEX 1
#define RWINDEXSCALE 3
#define RMAX_TYPE 0
#define NUMA_AWARE 1
#define INTEGER_RWUPDATE 1
#define UINT8_T_RWUPDATE 1
//#define TWO_PHASE_RWNUMA 1
#define NCOPY 4
#define NCORE 40
#ifndef NUMA_AWARE
#define LIGRA_NESTED 1
#endif
//#define GORDER 1
#define SOURCE_PLAW 1
//#define THRESARRAY 1
//#define NUMA_POP 1
//#define NUMA_TARGET 1
//#define USEPOUCH 1

//#define ITERATION_COUNT 1

//#define ITERATION_INFO 1
//#define BAG_PHASE_TIMER 1
//#define CHECK_INDEX_TARGETS 1
//#define DEBUG 1 
//#define TEST_PORTION_FP 1
//#define SPIN  0
//#define RACEINFO  1
//#define NON_ATOMIC_FP 1
//#define PHASE_TIMER 1
//#define NESTEDLOOP 1
//#ifndef FP_RESIDUE_DISTRIBUTION
//#endif
//#define FPWORKLOADCOUNT
//#define RWNESTEDLOOP 1
//#define RWPERMUTE 1
//#define DEBUGINFO 1
#define CHECK_TOPK 1
//#define FP_PERMUTE 1
//#define LOSSCHECK 1
//#define CACHEPADDING 0
//#define CACHEPAD_I32 1
#define CACHEPAD_U8 0
//#define TESTSINGLECORE 1
//#define FORWARD_PULL 1
//#define BINS_FOR_MC 1
//#define LOCAL_UPDATE_FP 1
//#define PULL 1
//#define LOCALTIMER 1
//#define BITSET 1
//#define SET_COVER 1
//#define SET_COVER 1
//#define SET_COVER_DUP 1
//#define CHECK_POP_NODES 1
#include <numa.h>
#include <atomic>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_list.h>
#include <climits>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <mutex>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <tr1/unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <sys/sysinfo.h>    // get cpu core number
#include <sched.h>


#ifdef USE_DOTMIX
#include <cilkpub/dotmix.h>
#endif

#ifdef USE_DOTMIX
cilkpub::DotMix global_rng(234909128);
#endif

#include "covergraph.h"
#include "bag.h"
//#include "rArray.h"
#include "cilkMax.h"
#include "pouch.h"
#include "pagerank.h"

#ifndef PROCESS_BLK_SERIALLY
#define PROCESS_BLK_SERIALLY 1
#endif

#ifndef PROCESS_EDGES_SERIALLY
#define PROCESS_EDGES_SERIALLY 0
#endif

#ifndef MUTEX
#define MUTEX 0
#endif

#ifndef USE_LOCK
//#define USE_LOCK 0
#endif

#define EDGE_THRESHOLD 128
bool run_topk=false;

inline bool file_exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}
// Compare and Swap
void atomic_add(atomic<double> &a, double b) {
    auto current = a.load();
    while (!a.compare_exchange_weak(current, current + b));
        ;
}
void atomic_add(atomic<int> &a, int b) {
    auto current = a.load();
    while (!a.compare_exchange_weak(current, current + b));
        ;
}
void atomic_add(atomic<uint16_t> &a, uint16_t b) {
    auto current = a.load();
    while (!a.compare_exchange_weak(current, current + b));
        ;
}
void atomic_add(atomic<uint8_t> &a, uint8_t b) {
    auto current = a.load();
    while (!a.compare_exchange_weak(current, current + b));
        ;
}

inline uint8_t Atomic_add(atomic<uint8_t> &a, uint8_t b) {
    while (1) {
        auto current = a.load();
        if (a.compare_exchange_weak(current, current + b)) {
            return current;
        } 
    }
}
cilk::reducer_opadd<int> race_num;  // use reducer to count race
// This version has return value that is used for discarding locks
inline double Atomic_add(atomic<double> &a, double b) {
    while (1) {
        auto current = a.load();
        if (a.compare_exchange_weak(current, current + b)) {
            return current;
        } 
        #ifdef RACEINFO
        else {
            race_num ++;
            //racecount ++;
        }
        #endif
    }
}
double Atomic_load_reset(atomic<double> &a) {
    while (1) {
        auto current = a.load();
        if (a.compare_exchange_weak(current, 0)) {
            return current;
        } 
    }
}
class MultiSources {
public:
    atomic<double>** residues;
    atomic<double>** reserves;
    double** push_counts;  // count nodes that push

    int num_core;
    atomic<double>* stats; // final results (times of each node being active)
    atomic<double>* final_counts; // final results (times of each node being active)
    atomic<long long int> total_num;

    MultiSources(int num_core) {
        this->num_core = num_core;
        residues = new atomic<double>*[num_core];
        reserves = new atomic<double>*[num_core];
        push_counts = new double*[num_core];
        for (int i = 0; i < num_core; i ++) {
            residues[i] = new atomic<double>[graph.n];
            reserves[i] = new atomic<double>[graph.n];
            memset(residues[i],  0, sizeof(atomic<double>)*graph.n);
            memset(reserves[i],  0, sizeof(atomic<double>)*graph.n);
            push_counts[i] = new double[graph.n];
            memset(push_counts[i],0,sizeof(double)*graph.n);
        }
        stats = new atomic<double>[graph.n];
        final_counts = new atomic<double> [graph.n];
        memset(stats, 0, sizeof(atomic<double>)*graph.n);
        memset(final_counts, 0, sizeof(atomic<double>) *graph.n);

        total_num.store(0);
    }

    // add result from core_id to the final results
    void merge_result_and_reset(int core_id) {
        // check result
        //double sum = 0;
        //double rsum= 0;
        //atomic<long long int> total_num_i;
        //total_num_i= 0;
        for (int i = 0; i < graph.n; i ++) {
            // check result
            //sum += residues[core_id][i] + reserves[core_id][i];
            //rsum+= residues[core_id][i];

            if (residues[core_id][i] > 0) {
                atomic_add(stats[i], 1);
                // check result
                //total_num_i += 1;
            }
            atomic_add(final_counts[i], push_counts[core_id][i]);
        }
        //total_num += total_num_i;
        //printf("%d sum of residues and reserves is %f\n", core_id ,sum);
        //printf("%d rsum is %f\n", core_id, rsum);
        memset(residues[core_id], 0, sizeof(atomic<double>)*graph.n);
        memset(reserves[core_id], 0, sizeof(atomic<double>)*graph.n);
        memset(push_counts[core_id],0,sizeof(double)*graph.n);
    }
};

class Result{
    public:
    int core_n;
    int length;
    int f_index;
    int r_index;
    double* f_results;
    double* r_results;
    double f_average;
    double r_average;
    double* rw_atomic;
    double* rw_setcover;
    int rwa_idx;
    int rws_idx;
    double rwa_average;
    double rws_average;
    // topk
    double* topk_results;
    double topk_average;
    int topk_index;
    // race count
    int race_idx;
    long* race_stats;
    long race_avg;
    vector<int> sources;

    Result(int n, int l) {
        core_n = n;
        length = l;
        f_results = new double[length];
        f_average =0.0;
        r_results = new double[length];
        r_average =0.0;
        f_index = 0;
        r_index = 0;

        rw_atomic = new double[length];
        rw_setcover = new double[length];
        rwa_idx = 0;
        rws_idx = 0;
        rwa_average = 0;
        rws_average = 0;
        
        // topk
        topk_results = new double[length];
        topk_index=0;
        topk_average = 0;
        // race
        race_idx=0;
        race_stats = new long[length];
        race_avg = 0;

        memset(f_results, 0, sizeof(double)*length);
        memset(r_results, 0, sizeof(double)*length);
        memset(topk_results, 0, sizeof(double)*length);
    }

    void set_sources(vector<int> & ss) {
        sources = ss;
    }

    void add_forward(double r) {
        f_results[f_index++] = r;    
    }
    void add_random(double r) {
        r_results[r_index++] = r;
    }
    void add_rw_setcover(double r) {
        rw_setcover[rws_idx++] = r;
    }
    void add_rw_atomic(double r) {
        rw_atomic[rwa_idx++] = r;
    }
    void add_topk(double r) {
        topk_results[topk_index++] = r;
    }
    void add_race_stat(long race_count) {
        race_stats[race_idx++] = race_count;
    }

    void print(int forward_p, int random_w) {
        for (int i = 0; i < length; i ++) {
            f_average += f_results[i];
            r_average += r_results[i];
            rws_average += rw_setcover[i]; 
            rwa_average += rw_atomic[i];
            race_avg += race_stats[i];

            topk_average += topk_results[i];
            //printf("%f %f\n", f_results[i], r_results[i]);
        }
        f_average /= length;
        r_average /= length;
        rws_average /= length;
        rwa_average /= length;
        race_avg /= length;

        topk_average /= length;

        if (forward_p)
            printf("Forward push executed in %f seconds with %d workers\n",
                    f_average/1000, core_n);
        #ifdef RANDOMWALK 
        if (random_w)
            printf("Random  walk executed in %f seconds with %d workers\n",
                    r_average/1000, core_n);

        //printf("setcover %f     atomic %f\n", rws_average, rwa_average);
        if (random_w && forward_p) {
            printf("Pafora executed in %f seconds with %d workers\n",
                    (f_average+r_average)/1000, core_n);
        }
        #endif
        #ifdef RACEINFO
        
        printf("PAFO race count is %ld\n", race_avg);
        #endif
        if (!run_topk) return;
        printf("topk queries executed in %f seconds with %d workers\n",
                    topk_average/1000, core_n);
    }
    void print_all_sources() {
        for (int i = 0; i < length; i ++) {
            printf("source = %d\n", sources[i]);
            printf("Forward push executed in %f seconds with %d workers\n",
                    f_results[i]/1000, core_n);
            printf("Random  walk executed in %f seconds with %d workers\n",
                    r_results[i]/1000, core_n);
        }
    }

    double print_source(int i) {
        printf("source = %d\n", sources[i]);
        printf("Forward push executed in %f seconds with %d workers\n",
                f_results[i]/1000, core_n);
        printf("Random  walk executed in %f seconds with %d workers\n",
                r_results[i]/1000, core_n);
        return f_results[i];
    }
};

map< int, vector< pair<int ,double> > > exact_topk_pprs;
string dataset;
int source_origin;
inline void load_exact_topk_ppr(){
    // string filename = dataset + "/livejournal.topk.pprs";
    string filename = dataset + "/orkut.topk.pprs";
    if(!file_exists_test(filename)){
        cout << "No exact topk ppr file" << filename << endl;;
        return;
    }
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> exact_topk_pprs;
    cout << "Exact topk pprs loaded" << endl;
}
class MyLock{
public:
    std::mutex lock;
    bool try_lock(){
        return lock.try_lock();
    }
    void unlock(){
        lock.unlock();
    }
};

//using namespace std::tr1;
using namespace std;
const string Green = "\033[0;32m";
const string Reset = "\033[0m";
const string Red = "\033[0;31m";

static string get_current_time_str() {
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
    std::string str(buffer);

    return str;
}



// variables
Graph graph;
double alpha;
PageRank pagerank;

int rounds = 100;
int nworker;    // number of cores
double source_portion=0;
double init_residue_value=0;
//const int GRAINSIZE = 2048;
//const int ALPHA = 0.2;
double rmax = 0.0;
double rsum = 1.0;
unsigned long long Omega = 0;

double* reserve_double;                // reserve array w.r.t s
atomic<double>* ppr_values;            // ppr array w.r.t s
//double* ppr_values;            // ppr array w.r.t s
//std::atomic_flag*    spinlocks;      // spinlocks to replace atomic double
atomic<double>* reserve_s;             // atomic reserve array w.r.t s
//atomic<double>* reserve_s_pop;       // atomic reserve array w.r.t s
//atomic<long long>* reserve_s_ll;     // atomic reserve array w.r.t s
//atomic<long>* reserve_s_l;           // atomic reserve array w.r.t s
//atomic<int>* reserve_s_int;          // atomic reserve array w.r.t s
//atomic<uint16_t>* reserve_s_u16;     // atomic reserve array w.r.t s
//atomic<uint8_t>* reserve_s_u8;       // atomic reserve array w.r.t s
// NUMA Aware reserve arrays
atomic<double>** reserve_numa;         // two dimensional array
//atomic<double>** reserve_numa_local; // two dimensional array local portion
//#ifdef INTEGER_RWUPDATE
atomic<uint8_t>** reserve_numa_u8;     // two dimensional array
atomic<int>** reserve_numa_i32;     // two dimensional array
//#endif
//double** reserve_numa_double;
int NUMA_NODES_SIZE = 1;
double* reserve_final;                 // final results

//double* double_incre;
double* residue_double;                // trivial array
//#ifndef FP_RESIDUE_DISTRIBUTION
#if SPIN
double* residue_s;
#else
atomic<double>* residue_s;             // residue array w.r.t s
#endif
//#endif
//atomic<double>* reserve_incre;       // reserves updated by random walk
cilk::reducer_opadd<double> s_residue;
//vector<int> first_enqueued;          // two marking arrays version
//vector<int> second_enqueued;
//bool* is_in_bag;                     // marking array
//#ifndef FP_RESIDUE_DISTRIBUTION
int source=0;
//#endif
int MAX_CORE_N = 100;
int dividingline=0;
int scan_type = 0;                    // 0 for set cover, 1 for sep_scan
//int forward_type = 0;               // 0 for bag, 1 for scan
int POUCH_SIZE= 8192;                 // for pouch


void permute_nodes(double);
#ifndef FP_RESIDUE_DISTRIBUTION
void init_random_walks(int);
#endif
void init_local_residues();
#ifdef FORWARD_PULL
void init_for_pull();
#endif
void set_up_bins_on_numa_nodes();
bool pair_sorter_large_first(pair<int, double> p1, pair<int, double> p2);
void topk_query_with_bound(int source, int k, Result* result);
// timer
struct timespec start, finish;
#ifdef PHASE_TIMER
struct timespec start_p1, finish_p1;
double total_time=0;
double total_time2=0;
int total_iter=0;
#endif
#ifdef BAG_PHASE_TIMER
struct timespec start_bag1, finish_bag1;
struct timespec start_bag2, finish_bag2;
double total_bag1_time=0;
double total_bag2_time=0;
#endif
#if defined(ITERATION_INFO)
double total_iter1=0;
double total_iter2=0;
double total_iter3=0;
#endif
#ifdef ITERATION_COUNT
double* iter_wl_counts;
int max_iter = 1000;
#endif
double g_elapsed=0.0;
// local timer
#ifdef LOCALTIMER
double local_elapsed=0.0;
#endif
// commandline arguments
int rwmode=0;
int with_rw_permutation=0;
double hybridsize=50.0;
#ifdef CHECK_POP_NODES
// Check popular nodes (most frequently hit)
double* global_pop_nodes=NULL;
double* ss_pop_nodes=NULL;
#endif
// GLOBAL DEBUGGING VARIABLES
int conflictCount = 0;
atomic<int> racecount(0);
int insertcount = 0;
#ifdef DEBUGINFO
cilk::reducer<MaxMonoid<int>> maxWork; 
long long* wlcount;
long long* wlcountotal;
long long totalMax = 0;
long long totalwl= 0;
int **real_locs;
int *real_locs_idx;
int real_loc_size=50;
#endif

// spinlock
/*
inline void spinlock(int i) {
    while (spinlocks[i].test_and_set(std::memory_order_acquire)) {
        
    }
}
inline void unlock(int i) {
    spinlocks[i].clear(std::memory_order_release);
}
*/

static void setWorkers(int n) {
  __cilkrts_end_cilk();
  //__cilkrts_init();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}
#ifndef FP_RESIDUE_DISTRIBUTION
// initialize residue_s array
void init_residue_s(int s) {
    race_num.set_value(0);
    racecount = 0;
    // may has been used for another source
    if (residue_s == NULL) {
        #if SPIN
        residue_s = new double[graph.n];
        spinlocks = new std::atomic_flag[graph.n];
        #else
        residue_s = new atomic<double>[graph.n];
        #endif

        //residue_double = new double[graph.n];
    }

    int o1 = graph.nodepointer[0][s];
    int o2 = graph.nodepointer[0][s+1];
    #ifdef DEBUGINFO
    printf("s = %d out degree (nodepointer[0][s]) is %d\n", s, o2-o1);
    printf("s = %d out degree (degree[0][s]) is %d\n", s, graph.degree[0][s]);
    #endif
    // r(s, v) <-- 0 for all v != s
#ifdef PARALLEL_INIT
    cilk_for(int i = 0; i < graph.n; i ++) {
        residue_s[i] = 0;
    }
#else
    #if SPIN
    memset(residue_s, 0, sizeof(double) * graph.n);
    #else
    memset(residue_s, 0, sizeof(atomic<double>) * graph.n);
    #endif
    //memset(residue_double, 0, sizeof(double)*graph.n);
#endif
    s_residue.set_value(0.0);
    // r(s, s) <-- 1
    #ifdef DEBUGINFO
    printf("s = %d\n", s);
    #endif

    if (s >= 0) {
        #if SPIN
        residue_s[s] = 1.0;
        #else
            #ifdef TEST_PORTION_FP
            if (init_residue_value == 0) {
                printf("please specify init residue value\n");
            }
        residue_s[s].store(source_portion);
            #else
        residue_s[s].store(1.0);
            #endif
        #endif
    }

    //residue_double[s] = 1;
    //init_local_residues();
    #ifdef FORWARD_PULL
    init_for_pull();
    #endif
}
#endif

// initialize reserve_s array
void init_s_reserve(int s) {

    // remove all the elements in source node's reserve vector
    if (reserve_s == NULL) {
        //reserve_s = new double[graph.n];
        if (reserve_s == NULL) {
            reserve_s = new atomic<double>[graph.n];
        }
        reserve_double = new double[graph.n];
        //NUMA aware reserve array
        #ifdef NUMA_AWARE
        NUMA_NODES_SIZE = numa_max_node()+1;
        int nworker = 80;
            #ifdef NCOPY
        NUMA_NODES_SIZE = NCOPY;
            #endif
        printf("NUMA_NODES_SIZE = %d\n", NUMA_NODES_SIZE);
        #ifndef INTEGER_RWUPDATE
        reserve_numa = new atomic<double>*[NUMA_NODES_SIZE];
        #endif
        //reserve_numa_double = new double*[NUMA_NODES_SIZE];
        #ifdef INTEGER_RWUPDATE
        reserve_numa_u8  = new atomic<uint8_t>*[NUMA_NODES_SIZE<<CACHEPAD_U8];
            #ifndef UINT8_T_RWUPDATE 
        reserve_numa_i32 = new atomic<int>*[NUMA_NODES_SIZE];
            #endif
        #endif
        for (int i = 0; i < NUMA_NODES_SIZE; i ++) {
            #ifndef INTEGER_RWUPDATE
            void* memory = numa_alloc_onnode(sizeof(atomic<double>) * graph.n, i%4);
            reserve_numa[i] = (atomic<double>*) memory;
            #endif
            //void* memory_double = numa_alloc_onnode(sizeof(double) * graph.n, i);
            //reserve_numa_double[i] = (double*) memory_double;
            #ifdef INTEGER_RWUPDATE
            void* memory_u8 = numa_alloc_onnode(sizeof(atomic<uint8_t>) * (graph.n), i);
            reserve_numa_u8[i]  = (atomic<uint8_t>*) memory_u8;
                #ifndef UINT8_T_RWUPDATE 
            void* memory_i32 = numa_alloc_onnode(sizeof(atomic<int>) * graph.n, i);
            reserve_numa_i32[i] = (atomic<int>*) memory_i32;
                #endif
            #endif
        }
        #endif
        //reserve_numa_local = new atomic<double>*[nworker];
        //for (int i = 0; i < nworker; i ++) {
        //    reserve_numa_local[i] = new atomic<double>[graph.target_divideline/5];
        //}
        reserve_final = new double[graph.n];
        memset(reserve_final, 0, sizeof(double)*graph.n);
    }
    // set to 0
#ifdef PARALLEL_INIT
    #pragma simd
    clk_for(int i = 0; i < graph.n; i ++) {
        ((double*)reserve_s)[i] = 0;
        ((double*)reserve_final)[i] = 0;
        ((double*)reserve_double)[i] = 0;
    }
#else
    memset(reserve_s, 0, sizeof(atomic<double>) * graph.n);
    //memset(reserve_s_pop, 0, sizeof(atomic<double>) * graph.n<<CACHEPADSIZE_D);
    //memset(reserve_s_ll, 0, sizeof(atomic<long long>) * graph.n);
    //memset(reserve_s_l, 0, sizeof(atomic<long>) * graph.n);
    //memset(reserve_s_int, 0, sizeof(atomic<int>) * graph.n<<CACHEPADSIZE_I);
    //memset(reserve_s_u16, 0, sizeof(atomic<uint16_t>)*graph.n);
    //memset(reserve_s_u8, 0, sizeof(atomic<uint8_t>)*graph.n);
    //memset(reserve_incre, 0, sizeof(atomic<double>) * graph.n);
    memset(reserve_double, 0, sizeof(double)*graph.n);
    if (reserve_final != NULL) delete[] reserve_final;
    reserve_final = new double[graph.n];
    memset(reserve_final, 0, sizeof(double) * graph.n);
#endif
    #ifdef NUMA_AWARE
    // NUMA Aware
    for (int i = 0 ; i < NUMA_NODES_SIZE; i ++) {
        #ifdef INTEGER_RWUPDATE
        memset(reserve_numa_u8[i], 0, sizeof(atomic<uint8_t>) * (graph.n<<CACHEPAD_U8));
            #ifndef UINT8_T_RWUPDATE 
        memset(reserve_numa_i32[i],0, sizeof(atomic<int>) * (graph.n));
            #endif
        #else
        memset(reserve_numa[i], 0, sizeof(atomic<double>) * graph.n);
        #endif
    }
    //for (int i = 0; i < nworker; i ++) {
    //    memset(reserve_numa_local[i], 0, 
    //            sizeof(atomic<double>) * graph.target_divideline/5);
    //}
    #endif
    
    #ifdef BINS_FOR_MC
    set_up_bins_on_numa_nodes();
    #endif
}

void delete_arrays() {
    #ifndef CACHE_SCHEDULE
    delete[] residue_s;
    delete[] reserve_s;
    #endif
    #ifdef NUMA_AWARE
    if (reserve_numa != NULL)
        for (int i = 0; i < NUMA_NODES_SIZE; i ++) {
            numa_free((void*)reserve_numa[i], sizeof(atomic<double>) * graph.n);
        }
    delete[] reserve_numa;
    #endif
    //delete[] reserve_final;
}

//////////////////////////////////Forward Push/////////////////////////////////
///////////////////// POUCH //////////////////////
#ifndef FP_RESIDUE_DISTRIBUTION
static inline
void put_in_pouch(Pouch*& pouch, int node, Bag<Pouch*>* bnext) {
    
    // First case: out degree is small enough
    long int load = graph.get_degree(node);
    int curload = pouch->get_load();

    if ( curload+load < POUCH_SIZE) {
        pouch->enqueue(node, load);
    } else if (curload + load == POUCH_SIZE) {
        pouch->enqueue(node, load);
        (*bnext).insert(pouch);
        pouch = new Pouch();
    } else {
        double  residue_i = Atomic_load_reset(residue_s[node]);
        // fill the gap in current pouch
        int gap = POUCH_SIZE - curload;
        int iter=0;
        #ifdef DEBUG 
        if (gap < 0) {
            printf("gap=%d curload=%d load=%ld\n", gap, curload, load);
            exit(1);
        }
        #endif
        pouch->set_subnode_2(node, 0, gap, residue_i);
        (*bnext).insert(pouch);
        pouch = new Pouch(); 
        
        // fill more pouches
        for (iter = gap; iter+POUCH_SIZE <= load;
                         iter += POUCH_SIZE ) {

            pouch->set_subnode_1(node,iter,iter+POUCH_SIZE, residue_i);
            (*bnext).insert(pouch);
            pouch = new Pouch();
        }
        // put the rest of the out neighbors in the pouch
        if (iter < load) { 
            pouch->set_subnode_1(node, iter, load, residue_i);
        }

    }
}

static inline void
propagate(int node, Bag_reducer<Pouch*> *next, Pouch*& pouch, 
          long start, long end , double residue){
    
    // get reducer's bag
    Bag<Pouch*>* bnext = &((*next).get_reference());
    double reserve = residue * graph.alpha;

    // update reserve of this node
    if ( start == 0) {
        
        #ifdef RACEINFO
        Atomic_add(reserve_s[node], reserve);
        #else
        atomic_add(reserve_s[node], reserve);
        #endif
        // CACHEPADDING
        //atomic_add(reserve_s[node<<CACHEPADSIZE_D], reserve);
    }
    
    // prepare for next iteration
    //long int num_out = graph.g[node].size();
    long num_out = graph.get_degree(node);
    //int num_out = graph.degree[0][node];
    double avg_push_residue = (residue - reserve)/num_out;

    //int offset = graph.r_edge_indices[node];
    long offset = graph.get_edges_start(node);
    //int offset = graph.nodepointer[0][node];
    start += offset;
    end += offset;
    for (long i = start; i < end; i++) {

        int neighbor = graph.edges[0][i];
        double preValue = Atomic_add(residue_s[neighbor], 
                                     avg_push_residue);
        // judge if this node should be inserted
        double residue_i = residue_s[neighbor].load();
        long long degree = graph.get_degree(neighbor);
        if (degree == 0) {
            residue_i =Atomic_load_reset(residue_s[neighbor]);
            double increase = residue_i * (1 - graph.alpha);
            //double preValue = Atomic_add(source,
            //                             increase);
            #ifdef SOURCEREDUCE
            s_residue += increase;
            #else
            #ifdef RACEINFO
            Atomic_add(reserve_s[node], reserve);
            #else
            atomic_add(residue_s[source], increase);
            #endif
            #endif
            //double threshold = graph.thres[source];
            //if (preValue < threshold &&
            //        residue_s[source].load() >= threshold) {
            //    // put neighbor in the pouch
            ////conflictCount++;
            //    put_in_pouch(pouch, source, bnext);
            //}
            #ifdef RACEINFO
            Atomic_add(reserve_s[neighbor], residue_i - increase);
            #else
            atomic_add(reserve_s[neighbor], residue_i - increase);
            #endif
            //CACHEPADDING
            //atomic_add(reserve_s[neighbor<<CACHEPADSIZE_D], residue_i*graph.alpha);
            
        } else if (residue_i >= degree*graph.rmax_scaled
                                 && preValue < degree*graph.rmax_scaled) {
            put_in_pouch(pouch, neighbor, bnext);
        }
        //}
    }
}

static inline void
pbfs_proc_Node(Pouch* const n[],
           int fillSize,
           Bag_reducer<Pouch*> *next)
{
    // TODO move newpouch to upper level function
    Pouch* newpouch = new Pouch();
    for (int j = 0; j < fillSize; ++j) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        const Pouch* pouch = n[j];

        // process the queue
        Queue* queue = pouch->queue;
    
        #ifdef DEBUGWL
        int size_count = queue->size;
        size_count += (pouch->subnode2.node == -1 ? 0 : 1);
        //maxWork->increase(size_count);
        int core_id = __cilkrts_get_worker_number(); 
        //wlcount[core_id] += size_count;
        #endif
        
        for (int i = 0; i < queue->size; i ++) {

            int node = queue->at(i);
            double residue = Atomic_load_reset(residue_s[node]);
            //propagate(node, next, newpouch, 
            //          0, graph.g[node].size(), residue);
            propagate(node, next, newpouch, 
                      0, graph.get_degree(node), residue);
            //propagate(node, next, newpouch, 
            //          0, graph.degree[0][node], residue);
            // DEBUG TODO need cache padding
            //wlcount[core_id] += graph.get_degree(node);
        }
        // process subnode 1
        Subnode subnode1 = pouch->subnode1;
        if (subnode1.node > -1) {
            propagate(subnode1.node, next, newpouch,
                      subnode1.start, subnode1.end, subnode1.residue);
            // #DEBUG TODO need cache padding
            #ifdef DEBUGWL
            wlcount[core_id] += subnode1.end - subnode1.start;
            #endif
        }
        // process subnode 2
        Subnode subnode2 = pouch->subnode2;
        if (subnode2.node > -1) {
            propagate(subnode2.node, next, newpouch,
                      subnode2.start, subnode2.end, subnode2.residue);
            // #DEBUG TODO need cache padding
            #ifdef DEBUGWL
            wlcount[core_id] += subnode2.end - subnode2.start;
            #endif
        }
        delete pouch;
    }
    if (newpouch->get_load() > 0) {
        Bag<Pouch*>* bnext = &((*next).get_reference());
        (*bnext).insert(newpouch);
//        printf("put newpouch in bag\n");
    }
}
static inline
void process_pennant(Pennant<Pouch*>* p, Bag_reducer<Pouch*>* next) {
    if (p->getLeft() != NULL)
        cilk_spawn process_pennant(p->getLeft(), next);

    if (p->getRight() != NULL)
        cilk_spawn process_pennant(p->getRight(), next);

    // Process the current element
    Pouch* const * n = p->getElements();
#if PROCESS_BLK_SERIALLY
    pbfs_proc_Node(n, BLK_SIZE, next);
#else
    #pragma cilk grainsize=1
    cilk_for (int i = 0; i < BLK_SIZE; i+=THRESHOLD) {
        // This is fine as long as THRESHOLD divides BLK_SIZE
        pbfs_proc_Node(n+i, THRESHOLD, next);
  }
#endif  // PROCESS_BLK_SERIALLY
    delete p;
    
}
static inline
void process_bag(Bag<Pouch*> *b,
                Bag_reducer<Pouch*>* next) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<Pouch*> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        cilk_spawn process_bag(b, next);
        process_pennant(p, next);

        cilk_sync;

    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        Pouch*const *n = b->getFilling();

#if PROCESS_BLK_SERIALLY
        pbfs_proc_Node(n, fillSize, next);
#else
        int extraFill = fillSize % THRESHOLD;
        cilk_spawn pbfs_proc_Node(n+fillSize-extraFill, extraFill, next);
        #pragma cilk grainsize = 1
        cilk_for (int i = 0; i < fillSize - extraFill; i += THRESHOLD) {
            pbfs_proc_Node(n+i, THRESHOLD, next);
        }
        cilk_sync;
#endif  // PROCESS_BLK_SERIALLY
    }
}

static inline void
pbfs_proc_Node(Pouch* const n[],
           int fillSize) {
    Pouch* newpouch = new Pouch();
    for (int j = 0; j < fillSize; ++j) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        const Pouch* pouch = n[j];
        // process subnode 2
        Subnode subnode2 = pouch->subnode2;
        if (subnode2.node > -1 && subnode2.start == 0) {
            int node = subnode2.node;
            double residual_i = subnode2.residue;
            atomic_add(residue_s[node], residual_i);
        }
        delete pouch;
    }
}
// put residual in pouch back
static inline
void process_pennant(Pennant<Pouch*>* p) {
    if (p->getLeft() != NULL)
        cilk_spawn process_pennant(p->getLeft());

    if (p->getRight() != NULL)
        cilk_spawn process_pennant(p->getRight());

    // Process the current element
    Pouch* const * n = p->getElements();
    pbfs_proc_Node(n, BLK_SIZE);
    delete p;
}
static inline
void put_residual_in_pouch_back(Bag<Pouch*> *b) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<Pouch*> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        cilk_spawn put_residual_in_pouch_back(b);
        process_pennant(p);

        cilk_sync;

    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        Pouch*const *n = b->getFilling();

        pbfs_proc_Node(n, fillSize);
    }
}
///////////////////END OF POUCH//////////////////////
#endif // FP_RESIDUE_DISTRIBUTION
/////////////////////BAG///////////////////////
#ifdef FP_RESIDUE_DISTRIBUTION

static inline void
pbfs_proc_Node(const int n[],
           int fillSize,
           Bag_reducer<int> *next, atomic<double>* residue_s,
           atomic<double>* reserve_s, int source, double* push_counts)
#else
static inline void
pbfs_proc_Node(const int n[],
           int fillSize,
           Bag_reducer<int> *next)
#endif
{
    #ifdef DEBUG
    //maxWork->increase(fillSize);
    //int core_id = __cilkrts_get_worker_number(); 
    //wlcount[core_id] += fillSize;
    //printf("fill size is %d\n", fillSize);
    #endif
    // get reducer's bag
    Bag<int>* bnext = &((*next).get_reference());
    for (int j = 0; j < fillSize; ++j) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        int node = n[j];

        #ifdef FP_RESIDUE_DISTRIBUTION
        push_counts[node] += 1;
        #endif

        #if SPIN
        spinlock(node);
        double residue = residue_s[node];
        residue_s[node] = 0;
        unlock(node);
        #else
        double residue = Atomic_load_reset(residue_s[node]);
        #endif
        
        // update reserve of this node
        double reserve = residue * graph.alpha;
        #ifdef RACEINFO
        Atomic_add(reserve_s[node], reserve);
        #else
        atomic_add(reserve_s[node], reserve);
        #endif
        // CACHEPADDING
        //atomic_add(reserve_s[node<<CACHEPADSIZE_D], reserve);
        //reserve_s[node] += reserve;
        // prepare for next iteration
        //long int num_out = graph.degree[0][node];
        long long num_out= graph.nodepointer[0][node+1] - graph.nodepointer[0][node];
        double avg_push_residue = (residue - reserve)/num_out;
        long long  start = graph.nodepointer[0][node];
        long long  end = graph.nodepointer[0][node+1];
        #ifdef NESTEDLOOP
        #pragma cilk grainsize = 8
        cilk_for (long long i = start; i < end; i++) {
        #else
        for (long long i = start; i < end; i++) {
        #endif
            // update the residues
            int neighbor = graph.edges[0][i];
//            #ifdef FP_RESIDUE_DISTRIBUTION
//            push_counts[neighbor] += 1;
//            #endif
            //int start_t = start;
            //int offset  = i;
            //int num_ou = num_out;
            #if SPIN
            spinlock(neighbor);
            double preValue = residue_s[neighbor];
            residue_s[neighbor] += avg_push_residue;
            unlock(neighbor);
            #else
            double preValue = Atomic_add(residue_s[neighbor], 
                                                       avg_push_residue);
            #endif
            // judge if this node should be inserted
            // double residue_i = residue_s[neighbor].load();
            double residue_i = residue_s[neighbor];
            long long degree = graph.get_degree(neighbor);
            if (degree == 0) {
                
                #if SPIN
                spinlock(neighbor);
                residue_i = residue_s[neighbor];
                residue_s[neighbor] = 0;
                unlock(neighbor);
                #else
                residue_i = Atomic_load_reset(residue_s[neighbor]);
                #endif

                double increase = residue_i * (1 - graph.alpha);
                //double preValue = Atomic_add(residue_s[source],
                //                             increase);
                #ifdef SOURCEREDUCE
                s_residue += increase;
                #else
                double preValue = Atomic_add(residue_s[source], increase);
                #ifdef THRESARRAY
                double threshold = graph.thres[source];
                #else
                long long degree_s = graph.get_degree(source);
                double threshold = degree_s * graph.rmax_scaled;
                #endif
                if (residue_i >= threshold && preValue < threshold) {
                    Bag<int>* bnext = &((*next).get_reference());
                    (*bnext).insert(source);
                }
                #endif
                #ifdef RACEINFO
                Atomic_add(reserve_s[neighbor], residue_i - increase);
                #else
                atomic_add(reserve_s[neighbor], residue_i - increase);
                #endif
                // CACHEPADDING
                //atomic_add(reserve_s[neighbor], residue_i - increase);
                //reserve_s[neighbor] += residue_i - increase;
            #ifdef THRESARRAY
            } else if (residue_i >= graph.thres[neighbor]
                                     && preValue < graph.thres[neighbor]) {
            #else
            } else if (residue_i >= graph.rmax_scaled * degree
                                     && preValue < graph.rmax_scaled * degree) {
            #endif
                Bag<int>* bnext = &((*next).get_reference());
                (*bnext).insert(neighbor);
            }
        }
    }
}

#ifdef FP_RESIDUE_DISTRIBUTION
void process_pennant(Pennant<int>* p, Bag_reducer<int>* next,
                    atomic<double>* residue_s, atomic<double>* reserve_s, int source, double* push_counts) {
    if (p->getLeft() != NULL)
        process_pennant(p->getLeft(), next, residue_s, reserve_s, source, push_counts);

    if (p->getRight() != NULL)
        process_pennant(p->getRight(), next, residue_s, reserve_s, source, push_counts);
    // Process the current element
    const int *n = p->getElements();
#if PROCESS_BLK_SERIALLY
    pbfs_proc_Node(n, BLK_SIZE, next, residue_s, reserve_s, source, push_counts);
#else
    #pragma cilk grainsize=1
    cilk_for (int i = 0; i < BLK_SIZE; i+=THRESHOLD) {
        // This is fine as long as THRESHOLD divides BLK_SIZE
        pbfs_proc_Node(n+i, THRESHOLD, next, residue_s, reserve_s, push_counts);
  }
#endif  // PROCESS_BLK_SERIALLY
    delete p;
    
}
#else
void process_pennant(Pennant<int>* p, Bag_reducer<int>* next) {
    // Process the current element
    if (p->getLeft() != NULL)
        cilk_spawn process_pennant(p->getLeft(), next);

    if (p->getRight() != NULL)
        cilk_spawn process_pennant(p->getRight(), next);
    // Process the current element
    const int *n = p->getElements();
#if PROCESS_BLK_SERIALLY
    pbfs_proc_Node(n, BLK_SIZE, next);
#else
    #pragma cilk grainsize=1
    cilk_for (int i = 0; i < BLK_SIZE; i+=THRESHOLD) {
        // This is fine as long as THRESHOLD divides BLK_SIZE
        pbfs_proc_Node(n+i, THRESHOLD, next);
  }
#endif  // PROCESS_BLK_SERIALLY
    delete p;
    
}
#endif

    
#ifdef FP_RESIDUE_DISTRIBUTION
void process_bag(Bag<int> *b,
                Bag_reducer<int>* next,
                atomic<double>* residue_s, atomic<double>* reserve_s, int source, double* push_counts) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<int> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        process_bag(b, next, residue_s, reserve_s, source, push_counts);
        process_pennant(p, next, residue_s, reserve_s, source, push_counts);
    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        const int *n = b->getFilling();

#if PROCESS_BLK_SERIALLY
        pbfs_proc_Node(n, fillSize, next, residue_s, reserve_s, source, push_counts);
#else
        int extraFill = fillSize % THRESHOLD;
        cilk_spawn pbfs_proc_Node(n+fillSize-extraFill, extraFill, next,
                                    residue_s, reserve_s, push_counts);
        #pragma cilk grainsize = 1
        cilk_for (int i = 0; i < fillSize - extraFill; i += THRESHOLD) {
            pbfs_proc_Node(n+i, THRESHOLD, next, residue_s, reserve_s, push_counts);
        }
        cilk_sync;
#endif  // PROCESS_BLK_SERIALLY
    }
}
#else

void process_bag(Bag<int> *b,
                Bag_reducer<int>* next) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<int> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        cilk_spawn process_bag(b, next);
        process_pennant(p, next);

        cilk_sync;

    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        const int *n = b->getFilling();

#if PROCESS_BLK_SERIALLY
        pbfs_proc_Node(n, fillSize, next);
#else
        int extraFill = fillSize % THRESHOLD;
        cilk_spawn pbfs_proc_Node(n+fillSize-extraFill, extraFill, next);
        #pragma cilk grainsize = 1
        cilk_for (int i = 0; i < fillSize - extraFill; i += THRESHOLD) {
            pbfs_proc_Node(n+i, THRESHOLD, next);
        }
        cilk_sync;
#endif  // PROCESS_BLK_SERIALLY
    }
}
#endif
#ifndef FP_RESIDUE_DISTRIBUTION
#ifdef USEPOUCH
int scan_forward(Bag_reducer<Pouch*>* next) {
#else
int scan_forward(Bag_reducer<int>* next) { 
#endif
    cilk::reducer_opadd<int> size(0);
    //#pragma cilk grainszie = 256
    cilk_for( int i = 0; i < graph.n; i ++) {
        int num_out = graph.get_degree(i);
        if (residue_s[i] > graph.rmax_scaled * num_out) {
            size += 1;
            int node = i;
            
            #if SPIN
            spinlock(node);
            double residue = residue_s[node];
            residue_s[node] = 0;
            unlock(node);
            #else
            double residue = Atomic_load_reset(residue_s[node]);
            #endif
             
            // update reserve of this node
            double reserve = residue * graph.alpha;
            #ifdef RACEINFO
            Atomic_add(reserve_s[node], reserve);
            #else
            atomic_add(reserve_s[node], reserve);
            #endif
            // CACHEPADDING
            //atomic_add(reserve_s[node<<CACHEPADSIZE_D], reserve);
            //reserve_s[node] += reserve;

            double avg_push_residue = (residue - reserve)/num_out;
            
            long long start = graph.get_edges_start(node);
            long long end = start + num_out;
            #ifdef NESTEDLOOP
            cilk_for (long i = start; i < end; i++) {
            #else
            for (long long i = start; i < end; i++) {
            #endif
                // propagate : update the residues
                int neighbor = graph.edges[0][i];

                #if SPIN
                spinlock(neighbor);
                residue_s[neighbor] += avg_push_residue;
                unlock(neighbor);
                #else
                #ifdef RACEINFO
                Atomic_add(residue_s[neighbor], avg_push_residue);
                #else
                atomic_add(residue_s[neighbor], avg_push_residue);
                #endif
                #endif

                double residue_i = residue_s[neighbor];
                
                if (graph.get_degree(neighbor) == 0) {
                    //atomic_add(residue_s[neighbor], -residue_i);
                    #if SPIN
                    spinlock(neighbor);
                    residue_i = residue_s[neighbor];
                    residue_s[neighbor] = 0;
                    unlock(neighbor);
                    #else
                    residue_i = Atomic_load_reset(residue_s[neighbor]);
                    #endif
                    
                    double increase = residue_i * (1 - graph.alpha);
                    #ifdef SOURCEREDUCE
                    s_residue += increase;
                    #else
                    #ifdef RACEINFO
                    Atomic_add(residue_s[source], increase);
                    #else
                    atomic_add(residue_s[source], increase);
                    #endif
                    #endif
                    //atomic_add(reserve_s[neighbor<<CACHEPADSIZE_D], residue_i - increase);
                    #ifdef RACEINFO
                    Atomic_add(reserve_s[neighbor], residue_i - increase);
                    #else
                    atomic_add(reserve_s[neighbor], residue_i - increase);
                    #endif
                    //reserve_s[neighbor] += residue_i - increase;
                    
                    continue;
                }
            }
            
        }
    }
    double source_residue = s_residue.get_value();
    s_residue.set_value(0.0);
    #if SPIN
    spinlock(source);
    residue_s[source] += source_residue;
    unlock(source);
    #else
    atomic_add(residue_s[source], source_residue);
    #endif

    int loadsize = size.get_value();
    if (loadsize < graph.n/hybridsize) {
        cilk_for( int i = 0; i < graph.n; i ++) {
            long long degree = graph.get_degree(i);
            if (residue_s[i] > graph.rmax_scaled * degree) {
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, i, bnext);
                #else
                Bag<int>* bnext = &((*next).get_reference());
                (*bnext).insert(i);
                #endif
            }
        }
    }
    return loadsize;
}
#endif

double * residues_c=NULL;
double * reserves_c=NULL;
//#ifndef FP_RESIDUE_DISTRIBUTION
int is_beginning = 1;
cilk::reducer_opadd<int> workload(0);
//#endif
#ifdef FP_RESIDUE_DISTRIBUTION
void inline propagate_from_a_node(int i, atomic<double>* residue_s,
                                         atomic<double>* reserve_s) {
#else
void inline propagate_from_a_node(int i) {
#endif
    #if FP_PERMUTE
    int node = graph.permute_map[i];
    long long degree = graph.get_degree(node);
    if (residue_s[node] > graph.rmax_scaled * degree){// && !graph.istrouble[node]) {
    #else
    long long degree = graph.get_degree(i);
    if (residue_s[i] > graph.rmax_scaled * degree){// && !graph.istrouble[i]) {
        #ifndef FP_RESIDUE_DISTRIBUTION
        workload += 1;
        #endif
        int node = i;
    #endif

        #ifdef CHECK_POP_NODES
        ss_pop_nodes[node] += 1;
        #endif

        #if SPIN
        spinlock(node);
        double residue = residue_s[node];
        residue_s[node] = 0;
        unlock(node);
        #else
        double residue = Atomic_load_reset(residue_s[node]);
        #endif
        // update reserve of this node
        double reserve = residue * graph.alpha;
        #ifdef RACEINFO
        Atomic_add(reserve_s[node], reserve);
        #else
        atomic_add(reserve_s[node], reserve);
        #endif
        // CACHEPADDING
        //atomic_add(reserve_s[node<<CACHEPADSIZE_D], reserve);
        //reserve_s[node] += reserve;

        int num_out= graph.nodepointer[0][node+1] - graph.nodepointer[0][node];
        double avg_push_residue = (residue - reserve)/num_out;
        
        long int offset = graph.get_edges_start(node);
        long long  start = offset;
        long long end = offset + num_out;
        //#pragma simd
        for (long i = start; i < end; i++) {
            // update the residues
            int neighbor = graph.edges[0][i];
            #if SPIN
            spinlock(neighbor);
            residue_s[neighbor] += avg_push_residue;
            unlock(neighbor);
            # else
            #ifdef RACEINFO
            Atomic_add(residue_s[neighbor], avg_push_residue);
            #else
            atomic_add(residue_s[neighbor], avg_push_residue);
            #endif
            #endif
            
//            #ifdef CHECK_POP_NODES
//            ss_pop_nodes[neighbor] += 1;
//            #endif
        }   
    }
}
// put residue_s into two non-atomic arrays
// scan the multiple set covers to do forwardpush in parallel
// propagate residue of nodes with no outgoing neighbors
#ifdef USEPOUCH
int sep_scan_forward(Bag_reducer<Pouch*>* next) {
#else
#ifdef FP_RESIDUE_DISTRIBUTION
int sep_scan_forward(Bag_reducer<int>* next, atomic<double>* residue_s,
                                             atomic<double>* reserve_s,
                                             int source) {
#else
int sep_scan_forward(Bag_reducer<int>* next) {
#endif
#endif
    
    #ifdef ITERATION_INFO
    total_iter2 += 1;
    #endif

    #ifdef FP_RESIDUE_DISTRIBUTION
    int workload=0;
    #else
    workload.set_value(0);
    #endif
    //cilk_for( int i = graph.topdegreeline; i < dividingline; i ++) {

    #ifdef LOCALTIMER
    // local timer
    struct timespec start_local, finish_local;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start_local);
    #endif 
    
    #ifdef FP_RESIDUE_DISTRIBUTION
    for( int i = 0; i < dividingline; i ++) {
    #else
    #pragma cilk grainsize=1024
    cilk_for( int i = 0; i < dividingline; i ++) {
    #endif
        //int core_id = sched_getcpu();
        #ifdef FP_RESIDUE_DISTRIBUTION
        //propagate_from_a_node(i, residue_s, reserve_s);
        #if FP_PERMUTE
        int node = graph.permute_map[i];
        long long degree = graph.get_degree(node);
        if (residue_s[node] > graph.rmax_scaled * degree){// && !graph.istrouble[node]) {
        #else
            #ifdef THRESARRAY
        if (residue_s[i] > graph.thres[i]) {
            #else
        long long degree = graph.get_degree(i);
        if (residue_s[i] > graph.rmax_scaled * degree){// && !graph.istrouble[i]) {
        //if (residue_s[i] > graph.thres[i]){// && !graph.istrouble[i]) {
            #endif
            workload += 1;
            int node = i;
        #endif

            #ifdef CHECK_POP_NODES
            ss_pop_nodes[node] += 1;
            #endif

            #if SPIN
            spinlock(node);
            double residue = residue_s[node];
            residue_s[node] = 0;
            unlock(node);
            #else
                #ifdef NON_ATOMIC_FP
            double residue = residue_double[node];
            residue_double[node] -= residue;
                #else
            double residue = Atomic_load_reset(residue_s[node]);
                #endif
            #endif
            // update reserve of this node
            double reserve = residue * graph.alpha;
            #ifdef RACEINFO
            Atomic_add(reserve_s[node], reserve);
            #else
            atomic_add(reserve_s[node], reserve);
            #endif
            // CACHEPADDING
            //atomic_add(reserve_s[node<<CACHEPADSIZE_D], reserve);
            //reserve_s[node] += reserve;

            int num_out= graph.nodepointer[0][node+1] - graph.nodepointer[0][node];
            double avg_push_residue = (residue - reserve)/num_out;
            
            long int offset = graph.get_edges_start(node);
            long long  start = offset;
            long long end = offset + num_out;
            //#pragma simd
            //if (num_out > 100000) {
            //    cilk_for(long i = start; i < end; i ++) {
            //        int neighbor = graph.edges[0][i];
            //        atomic_add(residue_s[neighbor], avg_push_residue);
            //    }
            //}
            //else {
                #ifdef NESTEDLOOP
                cilk_for (long i = start; i < end; i++) {
                #else
                for (long i = start; i < end; i++) {
                #endif
                    // update the residues
                    int neighbor = graph.edges[0][i];
                    #if SPIN
                    spinlock(neighbor);
                    residue_s[neighbor] += avg_push_residue;
                    unlock(neighbor);
                    # else
                    #ifdef RACEINFO
                    Atomic_add(residue_s[neighbor], avg_push_residue);
                    #else
                    atomic_add(residue_s[neighbor], avg_push_residue);
                    #endif
                    #endif
                    
//                    #ifdef CHECK_POP_NODES
//                    ss_pop_nodes[neighbor] += 1;
//                    #endif
                }   
            //}
        }
        #else
        propagate_from_a_node(i);
        #endif
        /* // two layer propagation    
        int node = i;
        long int start = graph.get_edges_start(node);
        long long end = graph.get_edges_start(node);
        for (int j = start; j < end; j ++) {
            propagate_from_a_node(j);
        }*/
        
    }
    
    #ifdef LOCALTIMER
    // local timer
    clock_gettime(CLOCK_REALTIME, &finish_local);
    elapsed = (finish_local.tv_sec - start_local.tv_sec)*1000000000;
    elapsed += (finish_local.tv_nsec - start_local.tv_nsec);
    local_elapsed += elapsed;
    //printf("Scan took %f ms\n", elapsed/1000000);
    #endif

    #ifdef FP_RESIDUE_DISTRIBUTION
    for( int i = dividingline; i < graph.n; i ++) {
    #else
    cilk_for( int i = dividingline; i < graph.n; i ++) {
    #endif
        double residue_i = residue_s[i];
        if (residue_i > 0) {
            residue_s[i] = 0;
            workload += 1;
            double increase = residue_i * (1 - graph.alpha);
            // update source residue
            #ifdef SOURCEREDUCE
            s_residue += increase;
            #else
            #ifdef RACEINFO
            Atomic_add(residue_s[source], increase);
            #else
            atomic_add(residue_s[source], increase);
            #endif
            #endif
            //reserve_s[i] += residue_i - increase;
            atomic_add(reserve_s[i],residue_i - increase);
            // CACHEPADDING
            //atomic_add(reserve_s[i<<CACHEPADSIZE_D],residue_i - increase);
        }
    }
    
    double source_residue = s_residue.get_value();
    s_residue.set_value(0.0);
    #if SPIN
    spinlock(source);
    residue_s[source] += source_residue;
    unlock(source);
    #else
    atomic_add(residue_s[source], source_residue);
    #endif
    //printf("scanned workload is %d\n", workoad.get_value());
    #ifdef FP_RESIDUE_DISTRIBUTION
    int loadsize = workload;
    #else
    int loadsize = workload.get_value();
    #endif
    if (loadsize < graph.n/hybridsize) {
        // put nodes in bag
        #ifdef FP_RESIDUE_DISTRIBUTION
        for (int i = 0; i < graph.n; i ++) {
        #else
        cilk_for (int i = 0; i < graph.n; i ++) {
        #endif
            #ifdef THRESARRAY
            if (residue_s[i] > graph.thres[i]) {
            #else
            int degree = graph.get_degree(i);
            if (residue_s[i] > graph.rmax_scaled * degree) {
            #endif
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, i, bnext);
                #else
                Bag<int>* bnext = &((*next).get_reference());
                (*bnext).insert(i);
                #endif
            }
        }
    }
    return loadsize;
}
#ifdef USEPOUCH
int sep_scan_forward_with_sched(Bag_reducer<Pouch*>* next) {
#else
int sep_scan_forward_with_sched(Bag_reducer<int>* next) {
#endif
    // BLOCK_SCHED
    graph.reset_block_masks();
    workload.set_value(0);

    #ifdef DEBUG
    printf("sep_scan_forward_with_sched is called\n");
    #endif


    #ifdef PHASE_TIMER
    total_iter ++;
    // phase 1
    clock_gettime(CLOCK_REALTIME, &start_p1);
    cilk_for (int i = 0; i < graph.phase1; i ++) {
        for (int j = 0; j < graph.NTHREADS; j ++) {
            int blk = graph.grain_blocks[i][j];
            int start = blk * graph.grain_size;
            int end = start + graph.grain_size;
            for (int k = start; k < end; k ++) {
                #ifdef THRESARRAY
                if (residue_s[k] > graph.thres[k]i) {
                #else
                long long degree = graph.get_degree(k);
                if (residue_s[k] > graph.rmax_scaled * degree) {
                #endif
                    workload += 1;
                    int node = k;
                    double residue = Atomic_load_reset(residue_s[node]);
                    // update reserve of this node
                    double reserve = residue * graph.alpha;
                    atomic_add(reserve_s[node], reserve);

                    int num_out= graph.nodepointer[0][node+1] - 
                                 graph.nodepointer[0][node];
                    double avg_push_residue = (residue - reserve)/num_out;
                    
                    long int offset = graph.get_edges_start(node);
                    long long  start = offset;
                    long long end = offset + num_out;
                    //#pragma simd
                    for (long i = start; i < end; i++) {
                        // update the residues
                        int neighbor = graph.edges[0][i];
                        atomic_add(residue_s[neighbor], avg_push_residue);
                    }   

                }
            }
        }
    }

    clock_gettime(CLOCK_REALTIME, &finish_p1);
    double elapsed = 0;
    elapsed = (finish_p1.tv_sec - start_p1.tv_sec)*1000000000;
    elapsed += (finish_p1.tv_nsec - start_p1.tv_nsec);
    //printf("elasped = %f\n", elapsed);
    total_time += elapsed/1000000000; // unit: second

    
    //phase 2
    clock_gettime(CLOCK_REALTIME, &start_p1);
    cilk_for (int j = graph.phase1*graph.NTHREADS; j < graph.NBLOCKS; j ++) {
        int blk = graph.grain_blocks[j/graph.NTHREADS][j%graph.NTHREADS];
        int start = blk * graph.grain_size;
        int end = start + graph.grain_size;
        for (int k = start; k < end; k ++) {
            #ifdef THRESARRAY
            if (residue_s[k] > graph.thres[k]) {
            #else
            long long degree = graph.get_degree(k);
            if (residue_s[k] > graph.rmax_scaled * degree) {
            #endif
                workload += 1;
                int node = k;
                double residue = Atomic_load_reset(residue_s[node]);
                // update reserve of this node
                double reserve = residue * graph.alpha;
                atomic_add(reserve_s[node], reserve);

                int num_out= graph.nodepointer[0][node+1] - 
                             graph.nodepointer[0][node];
                double avg_push_residue = (residue - reserve)/num_out;
                
                long int offset = graph.get_edges_start(node);
                long long  start = offset;
                long long end = offset + num_out;
                //#pragma simd
                for (long i = start; i < end; i++) {
                    // update the residues
                    int neighbor = graph.edges[0][i];
                    atomic_add(residue_s[neighbor], avg_push_residue);
                }   

            }
        }
    }

    clock_gettime(CLOCK_REALTIME, &finish_p1);
    elapsed = (finish_p1.tv_sec - start_p1.tv_sec)*1000000000;
    elapsed += (finish_p1.tv_nsec - start_p1.tv_nsec);
    //printf("elasped = %f\n", elapsed);
    total_time2 += elapsed/1000000000; // unit: second

    #else // PHASE_TIMER
    #pragma cilk grainsize=1
    cilk_for( int blk = 0; blk < graph.NBLOCKS; blk++) {
        int block_id = -1;
        
        // get core id
        int core_id = sched_getcpu()%NCORE;
        //int numa_node = numa_node_of_cpu( core_id);
        //if (numa_node != core_id % NCOPY) {
        //    printf("numa_node = %d, core_id = %d\n",numa_node, core_id);
        //}
        // get a block in its list first
        double idx = Atomic_add(graph.indices[core_id], 1);
        int idx_i = (int)lround(idx);
        if (idx_i < graph.nrounds[core_id]) {
            block_id = graph.grain_blocks[idx_i][core_id];
        } else {
        // then try to steal others block
            for (int tid = 0; tid < graph.NTHREADS; tid ++) {
                double idx = Atomic_add(graph.indices[tid], 1);
                int idx_i = (int)lround(idx);
                if (idx_i < graph.nrounds[tid]) {
                    block_id = graph.grain_blocks[idx_i][tid];
                    break;
                }
            }
        }
        //block_id = blk;
        if (block_id != -1) {
            int start = block_id * graph.grain_size; 
            int end = start + graph.grain_size;
            if (end >= graph.dividingline) end = graph.dividingline;

            for (int i = start; i < end; i ++) {
                #ifdef THRESARRAY
                if (residue_s[i] > graph.thres[i]) {
                #else
                long long degree = graph.get_degree(i);
                if (residue_s[i] > graph.rmax_scaled * degree){// && !graph.istrouble[i]) {
                #endif
                    workload += 1;
                    int node = i;
                    double residue = Atomic_load_reset(residue_s[node]);
                    // update reserve of this node
                    double reserve = residue * graph.alpha;
                    atomic_add(reserve_s[node], reserve);

                    int num_out= graph.nodepointer[0][node+1] - 
                                 graph.nodepointer[0][node];
                    double avg_push_residue = (residue - reserve)/num_out;
                    
                    long int offset = graph.get_edges_start(node);
                    long long  start = offset;
                    long long end = offset + num_out;
                    //#pragma simd
                    //if (block_id % 2 == 0) {
                        for (long i = start; i < end; i++) {
                            // update the residues
                            int neighbor = graph.edges[0][i];
                            atomic_add(residue_s[neighbor], avg_push_residue);
                        }   
                    //} else {
                    //    for (long i = end-1; i >= start; i--) {
                    //        // update the residues
                    //        int neighbor = graph.edges[0][i];
                    //        atomic_add(residue_s[neighbor], avg_push_residue);
                    //    }   
                    //}
                }
            }
        } else {
            //printf("block_id is -1, core id is %d\n",core_id);
        }
    }
    #endif

    cilk_for( int i = dividingline; i < graph.n; i ++) {
        double residue_i = residue_s[i];
        if (residue_i > 0) {
            residue_s[i] = 0;
            workload += 1;
            double increase = residue_i * (1 - graph.alpha);
            // update source residue
            #ifdef SOURCEREDUCE
            s_residue += increase;
            #else
            atomic_add(residue_s[source], increase);
            #endif
            atomic_add(reserve_s[i],residue_i - increase);
        }
    }
    
    double source_residue = s_residue.get_value();
    s_residue.set_value(0.0);
    atomic_add(residue_s[source], source_residue);
    int loadsize = workload.get_value();
    if (loadsize < graph.n/hybridsize) {
        // put nodes in bag
        cilk_for (int i = 0; i < graph.n; i ++) {
            #ifdef THRESARRAY
            if (residue_s[i] > graph.thres[i]) {
            #else
            long long degree = graph.get_degree(i);
            if (residue_s[i] > graph.rmax_scaled * degree) {
            #endif
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, i, bnext);
                #else
                Bag<int>* bnext = &((*next).get_reference());
                (*bnext).insert(i);
                #endif
            }
        }
    }
    return loadsize;
}
#ifdef FP_RESIDUE_DISTRIBUTION
void forward_push(int s, timespec& stoptime, atomic<double>* residue_s,
                                             atomic<double>* reserve_s,
                                             double* push_counts) {
    residue_s[s] = 1;
#else
void forward_push(int s, timespec& stoptime) {
#endif
    //printf("forward_push called %d\n", s);
    // check global popular nodes
    #ifdef CHECK_POP_NODES
    if (global_pop_nodes == NULL) {
        global_pop_nodes = new double[graph.n];
        ss_pop_nodes     = new double[graph.n];
    } else {
        memset(ss_pop_nodes, 0, sizeof(double)*graph.n);
    }
    #endif

    #ifdef LOCALTIMER
    // local timer
    local_elapsed = 0;
    #endif
    // DEBUG count size
    long long queuesize = 1;

    int source = s;
    int is_beginning = 1;

    //maxWork.get_value()->number = 0;
    int maxload = 0;
    // using two reducers to do the work
    #ifdef USEPOUCH
    Bag_reducer<Pouch*> *queue[2];
    Bag_reducer<Pouch*> b1;
    Bag_reducer<Pouch*> b2;
    queue[0] = &b1;
    queue[1] = &b2;
    #else
    Bag_reducer<int> *queue[2];
    Bag_reducer<int> b1;
    Bag_reducer<int> b2;
    queue[0] = &b1;
    queue[1] = &b2;
    #endif
    bool queuei = 1;

    // put the source node into one of the bags
    #ifdef USEPOUCH
    // put the source node into one of the bags
    Pouch* pouch_s = new Pouch();
    pouch_s->queue->enqueue(s);//, graph.get_degree(s));
    (*queue[1]).insert(pouch_s);
    #else
    (*queue[1]).insert(s);
    #endif
    // DEBUG
    vector<bool> isUpdated;
    int iter = 0;
    int scan_iter = 0;
    // TIMER
    g_elapsed = 0.0;
    int isdone = 0;
    int scan_size = 0;
    int forward_type = 0;
    // while there are nodes left in the bag
    int scan_iteration = 0;
    int scan_done = 0;
    int bag2_begin = 0;
    int scan_begin = 0;
    #ifdef BAG_PHASE_TIMER
    clock_gettime(CLOCK_REALTIME, &start_bag1);
    #endif
    //printf("FP start, source = %d\n",source);
    while ( !isdone ) {

        #ifdef DEBUG
        printf("iter = %d\n", iter);
        #endif

        #ifdef ITERATION_INFO
        if (scan_done == 0) total_iter1 += 1;
        else total_iter3 += 1;
        #endif

        #ifdef BAG_PHASE_TIMER
        if (scan_done == 1 && bag2_begin == 0) {
            bag2_begin = 1;
            clock_gettime(CLOCK_REALTIME, &start_bag2);
            printf("bag2 clock started\n");
        }
        if (scan_begin == 0 &&  (*queue[queuei]).numElements() > graph.n/hybridsize ) {
            scan_begin = 1;
            printf("bag1 clock stopped\n");
            clock_gettime(CLOCK_REALTIME, &finish_bag1);
        }
        #endif

        iter++;
        #ifdef HYBRID    
        if (scan_size > 0 || (*queue[queuei]).numElements() > graph.n/hybridsize) {
            #ifdef USEPOUCH
            if (forward_type == 0) {
                put_residual_in_pouch_back(&((*queue[queuei]).get_reference()));
            }
            #endif
            if(!(*queue[queuei]).isEmpty()) {
                (*queue[queuei]).clear();
            }
            forward_type = 1;
            #ifdef FP_RESIDUE_DISTRIBUTION
            scan_size = sep_scan_forward(queue[queuei], residue_s, reserve_s, s);
            scan_iter ++;
            #else
                #ifdef SEPARATE
                    #ifdef CACHE_SCHEDULE
            scan_size = sep_scan_forward_with_sched(queue[queuei]);
                    #else
            //scan_size = sep_scan_forward_with_sched(queue[queuei]);
            scan_size = sep_scan_forward(queue[queuei]);
                    #endif
                #else
            scan_size = scan_forward(queue[queuei]);
                #endif
            #endif
            scan_iteration += 1;
            #ifdef ITERATION_COUNT
            iter_wl_counts[iter] += scan_size;
            #endif
            if (scan_size < graph.n/hybridsize) {
                forward_type = 0;
                scan_size = 0;
                scan_done = 1;
                #ifdef DEBUG
                printf("scanning done. Last scanned size is %d\n", scan_size);
                printf("bag size is %d\n", (*queue[queuei]).numElements());
                #endif
            } else {
                #ifdef DEBUG
                printf("in scanning...scanned size is %d\n", scan_size);
                #endif
            }
            //continue;
        } else {
            isdone = 1; 
        }
        if (forward_type == 1)
            continue;
        #endif

        #ifdef DEBUG
        printf("iter = %d before clear outbag\n", iter);
        #endif
        if(!(*queue[!queuei]).isEmpty()) {
            (*queue[!queuei]).clear();
            
            #ifdef DEBUG
            printf("    iter = %d outbag cleared\n", iter);
            #endif
        }
        #ifdef DEBUG
        printf("iter = %d before process bag\n", iter);
        #endif
        #ifdef ITERATION_COUNT
        iter_wl_counts[iter] += (*queue[queuei]).numElements();
        #endif
        #ifdef FP_RESIDUE_DISTRIBUTION
        process_bag(&((*queue[queuei]).get_reference()), queue[!queuei], residue_s, reserve_s, source, push_counts);
        #else
        process_bag(&((*queue[queuei]).get_reference()), queue[!queuei]);
        #endif
        

        #ifdef DEBUG
        printf("iter = %d after process bag\n", iter);
        #endif

        double source_residue = s_residue.get_value();
        s_residue.set_value(0.0);
        long long degree_s = graph.get_degree(source);
        if (residue_s[source] >= graph.rmax_scaled * degree_s) {
            #if SPIN
            spinlock(source);
            residue_s[source] += source_residue;
            unlock(source);
            #else
            atomic_add(residue_s[source], source_residue);
            #endif
        } else {
            #if SPIN
            spinlock(source);
            residue_s[source] += source_residue;
            unlock(source);
            #else
            atomic_add(residue_s[source], source_residue);
            #endif
            long long degree_s = graph.get_degree(source);
            if (residue_s[source] > graph.rmax_scaled * degree_s) {
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag_reducer<Pouch*> *next = queue[!queuei];
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, source, bnext);
                (*bnext).insert(pouch);
                #else
                Bag_reducer<int> *next = queue[!queuei];
                Bag<int>* bnext = &((*next).get_reference());

                (*bnext).insert(source);
                #endif
            }
        }
        // prepare for next iteration
        queuei = !queuei;
        #ifdef DEBUG
        printf("bag size after process is %d\n", (*queue[queuei]).numElements());
        #endif
        isdone = (*queue[queuei]).isEmpty();
        
        #ifndef DEBUG
        continue;
        #endif
        queuesize += (*queue[queuei]).numElements();
        //maxload += maxWork.get_value()->number;
        //maxWork.get_value()->number = 0;

        
        #ifdef DEBUG
        totalwl = 0;
        long long iterMax = 0;
        int ncore = __cilkrts_get_nworkers();
        for (int i = 0; i < ncore; i ++) {
            //printf(" %lld ", wlcount[i]);
            totalwl += wlcount[i];
            wlcountotal[i] += wlcount[i];
            if (iterMax < wlcount[i]) 
                iterMax = wlcount[i];
            wlcount[i] = 0;
        }
        totalMax += iterMax;
        printf("iteration %d queue size is %d\n", iter, 
                (*queue[queuei]).numElements());//,maxload);
        printf("iter = %d\n", iter);
        int count = 0;
        double r_sum = 0.0;
        double u_sum = 0.0;
        
        for (int i = 0 ; i< graph.n; i += 1) {
            if (reserve_s[i] > 0.0) {
                count ++;
            }
            r_sum += reserve_s[i];
            u_sum += residue_s[i];
        }
        printf("reserve size is %d\n", count);

        printf("reserve sum is %.15f\n", r_sum);
        printf("residue sum is %.15f\n", u_sum);
        
        printf("total sum is %.15f\n", r_sum + u_sum);
        #endif        
    }
        
    #ifdef BAG_PHASE_TIMER
    clock_gettime(CLOCK_REALTIME, &finish_bag2);
    double elapsed = 0;
    elapsed = (finish_bag2.tv_sec - start_bag2.tv_sec)*1000000000;
    elapsed += (finish_bag2.tv_nsec - start_bag2.tv_nsec);
    elapsed = (finish_bag1.tv_sec - start_bag1.tv_sec)*1000000000;
    elapsed += (finish_bag1.tv_nsec - start_bag1.tv_nsec);
    //printf("elasped = %f\n", elapsed);
    total_bag2_time += elapsed/1000000000; // unit: second
    #endif

    #ifdef FP_RESIDUE_DISTRIBUTION
    //printf("%d scan iter is %d\n", s, scan_iter);
    //printf("%d iter is %d\n",s, iter);
    return;
    #endif
    // Stop the timer now
    clock_gettime(CLOCK_REALTIME, &stoptime);

    #ifdef CHECK_POP_NODES
    //double ss_sum_pop = 0;
    if (scan_iteration > 0) {
        for (int i = 0; i < graph.n; i ++) {
            //global_pop_nodes[i] += ss_pop_nodes[i]/scan_iteration;
            global_pop_nodes[i] += ss_pop_nodes[i];//scan_iteration;
            //ss_sum_pop += ss_pop_nodes[i]/scan_iteration;
        }
    }
    //printf("ss sum pop is %f\n", ss_sum_pop);
    #endif

    #ifdef LOCALTIMER
    // local timer
    printf("%f milliseconds\n", local_elapsed/1000000);
    #endif

    // get rsum
    rsum = 0;
    for (int i = 0 ; i < graph.n; i ++) {
        rsum += residue_s[i].load();
    }
    //cilk::reducer_opadd<double> residueSum(0.0);
    //cilk_for(int i = 0 ; i < graph.dividingline; i ++) {
    //    residueSum += residue_s[i].load();
    //}
    //rsum = residueSum.get_value();

    #ifdef FPWORKLOADCOUNT
    // get heavy work load
    int count = 0;
    for (int i = 0 ; i< graph.n; i += 1) {
        if (reserve_s[i] > 0.0) {
            count ++;
        }
    }
    printf("source is %d reserve size is %d portion = %f\n", graph.sep_reverse_map[source], count, double(count)/double(graph.n));

    double elapsed = 0;
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Forward push calculated in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());
    printf("\n");
    #endif

    #ifdef DEBUGINFO
    // calculate rsum
    
    printf("rsum before dividingline = %f\n", rsum);

    //residueSum.set_value(0.0);
    //cilk_for(int i = 0 ; i < graph.n; i ++) {
    //    residueSum += residue_s[i].load();
    //}
    //printf("rsum is %f\n", residueSum.get_value());
    //return;
    totalwl = 0;
    int ncore = __cilkrts_get_nworkers();
    for (int i = 0; i < ncore; i ++) {
        //printf(" %lld", wlcountotal[i]);
        totalwl += wlcountotal[i];
    }
    printf("Max workload line is %lld\n", totalMax);
    printf("total workload is %lld\n", totalwl);
    totalMax = 0;
    
    //delete[] is_in_bag;
    //printf("Source node = %d out degree = %ld\n",source, graph.get_degree(source));
    printf("iter = %d\n", iter);
    double r_sum = 0.0;
    double u_sum = 0.0;
    
    for (int i = 0 ; i< graph.n; i += 1) {
        if (is_beginning == 0) {
            r_sum += reserve_s[i];
            u_sum += residues_c[i];
        } else {
            r_sum += reserve_s[i];
            u_sum += residue_s[i];
        }
    }
    printf("%lld Nodes processed\n", queuesize);
    printf("Max workload is %d\n", maxload);

    printf("reserve sum is %.15f\n", r_sum);
    printf("residue sum is %.15f\n", u_sum);
    
    printf("total sum is %.15f\n", r_sum + u_sum);
    // output results

    rsum = u_sum;
    #endif
    #ifdef RACEINFO
    //printf("total race count is %d %d\n",race_num.get_value(), 
    //                        racecount.load());
    #endif
}

//&&&&&&&&&&&&&&&&&&&&&&&& LOCAL UPDATE FORWARD PUSH &&&&&&&&&&&&&&&&&&&&&&&&&&&

#ifdef LOCAL_UPDATE_FP
static const int BLOCK_SIZE = 1 << 10;
class LocalResidue{
public:
    double* residues;
    // block markings
    bool* flags;
    int flag_length;
    LocalResidue() {
        residues = new double[graph.n];
        flag_length = graph.n/BLOCK_SIZE + 1;
        flags    = new bool[flag_length];
    }
    void clear() {
        memset(residues, 0, sizeof(double)*graph.n);
        memset(flags, false, sizeof(bool)*flag_length);
    }
};

LocalResidue* local_residues = NULL;

void init_local_residues() {
    int nworkers = 80;
    if ( local_residues == NULL) {
        local_residues = new LocalResidue[nworkers];
    } else {
        // TODO try cilk for
        for (int i = 0; i < nworkers; i ++) {
            local_residues[i].clear();
        }
    }

}
static inline void
pbfs_proc_Node(const int n[],
           int fillSize) {

    #ifdef DEBUG
    //maxWork->increase(fillSize);
    //wlcount[core_id] += fillSize;
    #endif
    int core_id = __cilkrts_get_worker_number(); 
    double* local_residue_array = local_residues[core_id].residues;
    bool* local_flags = local_residues[core_id].flags;
    // get reducer's bag
    for (int j = 0; j < fillSize; ++j) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        int node = n[j];
        double residue = residue_double[node]; 
        residue_double[node] = 0;
        // update reserve of this node
        double reserve = residue * graph.alpha;
        reserve_double[node] += reserve;
        
        // prepare message lists
        //long int num_out = graph.degree[0][node];
        long long  start = graph.nodepointer[0][node];
        long long  end = graph.nodepointer[0][node+1];
        long long num_out= end - start;
        double avg_push_residue = (residue - reserve)/num_out;
        // process each neighbor
        for (long long i = start; i < end; i++) {
            // update to local residue array
            int neighbor = graph.edges[0][i];
            local_residue_array[neighbor] += avg_push_residue; 
            local_flags[neighbor/BLOCK_SIZE] = true;
        }
    }
}

void process_pennant(Pennant<int>* p) {
    if (p->getLeft() != NULL)
        cilk_spawn process_pennant(p->getLeft());

    if (p->getRight() != NULL)
        cilk_spawn process_pennant(p->getRight());

    // Process the current element
    const int *n = p->getElements();
#if PROCESS_BLK_SERIALLY
    pbfs_proc_Node(n, BLK_SIZE);
#else
    #pragma cilk grainsize=1
    cilk_for (int i = 0; i < BLK_SIZE; i+=THRESHOLD) {
        // This is fine as long as THRESHOLD divides BLK_SIZE
        pbfs_proc_Node(n+i, THRESHOLD);
  }
#endif  // PROCESS_BLK_SERIALLY
    delete p;
    
}
void process_bag(Bag<int> *b) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<int> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        cilk_spawn process_bag(b);
        process_pennant(p);

        cilk_sync;

    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        const int *n = b->getFilling();

#if PROCESS_BLK_SERIALLY
        pbfs_proc_Node(n, fillSize);
#else
        int extraFill = fillSize % THRESHOLD;
        cilk_spawn pbfs_proc_Node(n+fillSize-extraFill, extraFill, next);
        #pragma cilk grainsize = 1
        cilk_for (int i = 0; i < fillSize - extraFill; i += THRESHOLD) {
            pbfs_proc_Node(n+i, THRESHOLD);
        }
        cilk_sync;
#endif  // PROCESS_BLK_SERIALLY
    }
}

#ifndef FP_RESIDUE_DISTRIBUTION
void forward_push_local_update(int s, timespec& stoptime) {
    

    Bag_reducer<int> bag;
    bag.insert(s);
    bool is_done = false;
    int core_num = __cilkrts_get_nworkers();
    int block_num = graph.n/BLOCK_SIZE + 1;
    while(!is_done) {
        // process bag in parallel
        process_bag(&(bag.get_reference()));

        // merge result
        for (int i = 0; i < core_num; i ++) {
            // handle message list in each worker
            double* residues_i_core = local_residues[i].residues;
            bool*   flags = local_residues[i].flags;
            cilk_for(int j = 0; j < block_num; j ++) {
                int start = j * BLOCK_SIZE;
                if (flags[j]) {
                    for (int k = start; k < start+BLOCK_SIZE; k ++) {
                        residue_double[k] += residues_i_core[k];
                    }
                }
            }
        }
        #pragma cilk grainsize = 1
        cilk_for(int i = 0 ; i < core_num; i ++) {
            local_residues[i].clear();
        }
        // generate new bag
        bag.clear();
        cilk_for (int i = 0; i < graph.dividingline; i ++) {
            double residue = residue_double[i];
            long long degree = graph.get_degree(i);
            if (residue > graph.rmax_scaled * degree) {
                bag.insert(i);
            }
        }
        // process nodes without out degree
        cilk_for( int i = dividingline; i < graph.n; i ++) {
            double residue_i = residue_double[i];
            // TODO remove if-condition and see the effect
            if (residue_i > 0) {
                residue_double[i] = 0;
                double increase = residue_i * (1 - graph.alpha);
                // update source residue by reducer
                s_residue += increase;
                reserve_double[i] += residue_i - increase;
            }
        }
    
        double source_residue = s_residue.get_value();
        residue_double[source] += source_residue;
        s_residue.set_value(0.0);

        is_done = bag.isEmpty(); 
    }
    // Stop the timer now
    clock_gettime(CLOCK_REALTIME, &stoptime);

    // get rsum
    cilk::reducer_opadd<double> residueSum(0.0);
    cilk_for(int i = 0 ; i < graph.dividingline; i ++) {
        residueSum += residue_double[i];
    }
    rsum = residueSum.get_value();
    // set up reserve array for monte carlo
    cilk_for( int i = 0; i < graph.n; i ++) {
        reserve_s[i] = reserve_double[i];
    }
    #ifdef DEBUGINFO
    // calculate rsum
    cilk::reducer_opadd<double> reserveSum(0.0);
    
    cilk_for(int i = 0 ; i < graph.n; i ++) {
        reserveSum += reserve_double[i];
    }
    printf("rsum is %f\n", residueSum.get_value());
    //return;
    double r_sum = residueSum.get_value();
    double u_sum = reserveSum.get_value();
    printf("residue sum is %.15f\n", residueSum.get_value());
    printf("reserve sum is %.15f\n", u_sum);
    
    printf("total sum is %.15f\n", r_sum + u_sum);
    // output results
    //exit(0);
    #endif
}
#endif

#endif // LOCAO_UPDATE_FP
//===================== END OF MESSAGE LIST FORWARD PUSH ======================
//&&&&&&&&&&&&&&&&&&&&&&&& PULL FORWARD PUSH &&&&&&&&&&&&&&&&&&&&&&&&&&&
#ifdef FORWARD_PULL
typedef uint8_t BYTE;
typedef boost::dynamic_bitset<> BitSet;
int firsttime = 1;
#ifdef BITSET
BitSet* pull_active_local = NULL;
BitSet* push_active_local = NULL;
BitSet pull_active;
BitSet push_active;
#else
BYTE* pull_active = NULL;
BYTE* push_active = NULL;
#endif
double * residue_pull;
double * residue_push;
cilk::reducer_opadd<int> push_sum(0);
cilk::reducer_opadd<int> pull_sum(0);

void init_for_pull() {
    int nworkers = 80;
    if ( firsttime == 1) {
        firsttime = 0;
        #ifdef BITSET
        int nCore = MAX_CORE_N;
        pull_active_local = new BitSet[nCore];
        push_active_local = new BitSet[nCore];
        for (int i = 0; i < nCore; i ++) {
            pull_active_local[i] = BitSet(graph.n);
            push_active_local[i] = BitSet(graph.n);
        }
        pull_active = BitSet(graph.n);
        push_active = BitSet(graph.n);
        #else
        pull_active = new BYTE[graph.n];
        push_active = new BYTE[graph.n];
        #endif
        residue_pull = new double[graph.n];
        residue_push = new double[graph.n];
    } else {
        #ifdef BITSET
        int nCore = MAX_CORE_N;
        for (int i = 0; i < nCore; i ++) {
            pull_active_local[i] = BitSet(graph.n);
            push_active_local[i] = BitSet(graph.n);
        }
        //pull_active = BitSet(graph.n);
        //push_active = BitSet(graph.n);
        #else
        memset(pull_active, 0, graph.n);    // Both sizeof(BYTE)*graph.n and graph.n work
        memset(push_active, 0, sizeof(BYTE)*graph.n);
        #endif
        memset(residue_pull, 0, sizeof(double) * graph.n);
        memset(residue_push, 0, sizeof(double) * graph.n);
    }
}
// used for marking pull active
/*
static inline void
pbfs_proc_Node(const int n[],
           int fillSize, int number) {

    #ifdef DEBUG
    maxWork->increase(fillSize);
    //wlcount[core_id] += fillSize;
    #endif
    int core_id = __cilkrts_get_worker_number(); 
    double* local_residue_array = local_residues[core_id].residues;
    bool* local_flags = local_residues[core_id].flags;
    // get reducer's bag
    for (int j = 0; j < fillSize; ++j) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        int node = n[j];
        push_active[node] = 1;
        double residue = residue_push[node]; 
        // update reserve of this node
        double reserve = residue * graph.alpha;
        reserve_double[node] += reserve;
        
        // prepare message lists
        //long int num_out = graph.degree[0][node];
        long long  start = graph.nodepointer[0][node];
        long long  end = graph.nodepointer[0][node+1];
        long long num_out= end - start;
        double avg_push_residue = (residue - reserve)/num_out;
        // process each neighbor
        residue_push[node] = avg_push_residue;
        for (long long i = start; i < end; i++) {
            // update to byte array
            int neighbor = graph.edges[0][i];
            pull_active[neighbor] = 1;
        }
    }
}

void process_pennant(Pennant<int>* p, int number) {
    if (p->getLeft() != NULL)
        cilk_spawn process_pennant(p->getLeft(),1);

    if (p->getRight() != NULL)
        cilk_spawn process_pennant(p->getRight(),1);

    // Process the current element
    const int *n = p->getElements();
#if PROCESS_BLK_SERIALLY
    pbfs_proc_Node(n, BLK_SIZE,1);
#else
    #pragma cilk grainsize=1
    cilk_for (int i = 0; i < BLK_SIZE; i+=THRESHOLD) {
        // This is fine as long as THRESHOLD divides BLK_SIZE
        pbfs_proc_Node(n+i, THRESHOLD,1);
  }
#endif  // PROCESS_BLK_SERIALLY
    delete p;
    
}
void process_bag(Bag<int> *b, int number) {
    if (b->getFill() > 0) {
        // Split the bag and recurse
        Pennant<int> *p = NULL;

        b->split(&p); // Destructive split, decrements b->getFill()
        cilk_spawn process_bag(b,1);
        process_pennant(p,1);

        cilk_sync;

    } else {
       
        int fillSize = 0;
        fillSize = b->getFillingSize();
        const int *n = b->getFilling();

#if PROCESS_BLK_SERIALLY
        pbfs_proc_Node(n, fillSize,1);
#else
        int extraFill = fillSize % THRESHOLD;
        cilk_spawn pbfs_proc_Node(n+fillSize-extraFill, extraFill, next);
        #pragma cilk grainsize = 1
        cilk_for (int i = 0; i < fillSize - extraFill; i += THRESHOLD) {
            pbfs_proc_Node(n+i, THRESHOLD);
        }
        cilk_sync;
#endif  // PROCESS_BLK_SERIALLY
    }
}
*/
#ifndef FP_RESIDUE_DISTRIBUTION
int sep_scan_pull(Bag_reducer<Pouch>* next) {
    cilk::reducer_opadd<int> workload(0);
    cilk::reducer_opadd<double> s_residue(0.0);
    // process nodes without out degree
    cilk_for( int i = dividingline; i < graph.n; i ++) {
        double residue_i = residue_push[i];
        // TODO remove if-condition and see the effect
        if (residue_i > 0) {
            double increase = residue_i * (1 - graph.alpha);
            // update source residue by reducer
            s_residue += increase;
            reserve_double[i] += residue_i - increase;
            residue_push[i] = 0;
        }
    }
    double source_residue = s_residue.get_value();
    residue_push[source] += source_residue;
    s_residue.set_value(0.0);
    // push
    #ifdef LOCALTIMER
    // local timer
    struct timespec start_local, finish_local;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &start_local);
    #endif
    #pragma cilk grainsize=128*64
    cilk_for (int i = 0; i < graph.dividingline; i ++) {
        push_active[i] = 0;
        //int core_id = __cilkrts_get_worker_number(); 
        double residue = residue_push[i];
        long long degree = graph.get_degree(i);
        if (residue > graph.rmax_scaled * degree) {
            workload+=1;
            push_active[i] = 1;
            double reserve = residue * graph.alpha;
            reserve_double[i] += reserve;
            
            long long  start = graph.nodepointer[0][i];
            long long  end = graph.nodepointer[0][i+1];
            long long num_out = end - start;
            double avg_push_residue = (residue - reserve)/num_out;
            residue_push[i] = avg_push_residue;
            // process each neighbor
            #pragma simd
            for (long long j = start; j < end; j++) {
                // update to byte array
                int neighbor = graph.edges[0][j];
                #ifdef BITSET
                pull_active[neighbor] = 1;
                #else
                pull_active[neighbor] = 1;
                #endif
            }
            //push_sum += num_out;
        }
    }
    #ifdef LOCALTIMER
    // local timer
    clock_gettime(CLOCK_REALTIME, &finish_local);
    elapsed = (finish_local.tv_sec - start_local.tv_sec)*1000000000;
    elapsed += (finish_local.tv_nsec - start_local.tv_nsec);
    local_elapsed += elapsed;
    #endif
    // pull
    cilk_for(int i = 0; i < graph.n; i ++) {
        if (pull_active[i] == 1) {
            pull_active[i] = 0;
//            printf("node %d start pulling...\n", i);

            pull_sum += graph.gr_sep[i].size();
            // traverse its in neighbors
            for (int j = 0; j < graph.gr_sep[i].size(); j ++) { 
                // for each of its in neighbor
                int in_neighbor = graph.gr_sep[i][j];
                    
                if (push_active[in_neighbor]) {
                    residue_pull[i] += residue_push[in_neighbor];
//                    printf("%d pulls from %d\n", i, in_neighbor);
                }
            }
        }
    }
    
    cilk_for(int i = 0; i < graph.dividingline; i ++) {
        if (!push_active[i]) {
            residue_pull[i] += residue_push[i];
        }
        residue_push[i] = 0;
    }
    double* tmp = residue_push;
    residue_push = residue_pull;
    residue_pull = tmp;
    int loadsize = workload.get_value();
    if (loadsize < graph.n/hybridsize) {
        // put nodes in bag
        cilk_for( int i = dividingline; i < graph.n; i ++) {
            double residue_i = residue_push[i];
            if (residue_i > 0) {
                residue_push[i] = 0;
                workload += 1;
                double increase = residue_i * (1 - graph.alpha);
                // update source residue
                s_residue += increase;
                reserve_double[i] += residue_i - increase;
            }
        }
        double source_residue = s_residue.get_value();
        residue_push[source] += source_residue;
        s_residue.set_value(0.0);

        cilk_for (int i = 0; i < graph.n; i ++) {
            residue_s[i] = residue_push[i];
            reserve_s[i] = reserve_double[i];
            long long degree = graph.get_degree(i);
            if (residue_s[i] > graph.rmax_scaled * degree) {
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, i, bnext);
                #else
                Bag<int>* bnext = &((*next).get_reference());
                (*bnext).insert(source);
                #endif
            }
        }
    }
    return loadsize;
}
#endif

#ifdef SET_COVER
//  set cover
int set_cover_forward_push() {
    #ifdef DEBUGINFO
    printf("set cover scan called\n");
    #endif
    cilk::reducer_opadd<int> size(0);
    // get a non-atomic copy of reserve and residue values
    if (residues_c == NULL) {
        #ifdef DEBUGINFO
        printf("set cover first called\n");
        #endif
        residues_c = new double[graph.n];
    }
    if (is_beginning == 1) {
        cilk_for (int i = 0 ; i < graph.n; i ++) {
            residues_c[i] = residue_s[i];
        }
        is_beginning = 0;
    }
    //double * newresidues = new double[graph.n];
    //memset(newresidues, 0, sizeof(double) * graph.n);
    
    for (int i = 0; i < rounds; i ++) {
        // scan the vector
        cilk_for (int v = 0; v < graph.pr_multi_results[i].size(); v ++) {
            // get the node in set covers
            int node = graph.pr_multi_results[i][v];
            double residue = residues_c[node];
            long long degree = graph.get_degree(node);
            if (degree == 0 && residue > 0) {
                residues_c[node] = 0;
                //residues_c[node] -= residue;
                size += 1;
                double increase = residue * (1 - graph.alpha);
                // update source residue
                s_residue += increase;
                //reserve_s[i] += residue - increase;
                atomic_add(reserve_s[i],residue - increase);
                
            } else if (residue > graph.rmax_scaled * degree) {
                //residues_c[node] = 0;
                residues_c[node] -= residue;
                size += 1;
                double reserve = residue * graph.alpha;
                //reserve_s[node] += reserve;
                atomic_add(reserve_s[node], reserve);
                int num_out= graph.nodepointer[0][node+1] - graph.nodepointer[0][node];
                double avg_push_residue = (residue - reserve)/num_out;
                
                long long offset = graph.nodepointer[0][node];
                long long  start = offset;
                long long end = offset + num_out;
                #pragma simd
                for (long ii = start; ii < end; ii++) {
                //cilk_for (long ii = start; ii < end; ii++) {
                    // update the residues
                    int neighbor = graph.edges[0][ii];
                    //newresidues[neighbor] += avg_push_residue;
                    residues_c[neighbor] += avg_push_residue;
                }
            }
        }
    }
    //update residue_s array
    //cilk_for(int i = 0; i < graph.n; i ++) {
    //    residues_c[i] += newresidues[i];
    //}
    double source_residue = s_residue.get_value();
    s_residue.set_value(0.0);
//    Atomic_add(residue_s[source], source_residue);
    residues_c[source] += source_residue;
    //delete[] newresidues;
    if (size.get_value() < nworker * 32) {
        // set cover is done, switch to sep_scan
        scan_type = 1;
        // controls results stats
        is_beginning = 1;
    
        cilk_for(int i = 0 ; i < graph.n; i ++) {
            residue_s[i] = residues_c[i];
        }
        forward_type = 1;
    }
    double s1 = 0;
    double s2 = 0;
    for (int i = 0 ; i < graph.n; i ++) {
        s1 += residues_c[i];
        s2 += reserve_s[i];
    }
    #ifdef DEBUGINFO
    printf("residue sum = %.10f reserve sum = %.10f total %.10f\n",s1,s2,s1+s2);
    #endif
    return size.get_value();
}
#endif // SET_COVER
#ifndef FP_RESIDUE_DISTRIBUTION
void forward_push_pull_update(int s, timespec& stoptime) {
            
    //printf("pull\n");
    push_sum.set_value(0);
    pull_sum.set_value(0);
    #ifdef LOCALTIMER
    // local timer
    local_elapsed = 0;
    #endif
    
    int forward_type = 0;
    long long queuesize = 1;
    source = s;
    Bag_reducer<int> bag;
    bag.insert(s);
    bool is_done = false;
    cilk::reducer_opadd<int> done;
    int core_num = __cilkrts_get_nworkers();
    int block_num = graph.n/BLOCK_SIZE + 1;
    
    s_residue.set_value(0.0);
    //residue_push[s] = 1;

    Bag_reducer<int> *queue[2];
    Bag_reducer<int> b1;
    Bag_reducer<int> b2;
    queue[0] = &b1;
    queue[1] = &b2;
    bool queuei = 1;

    (*queue[1]).insert(s);
    int iter = 0;
    int first_scan = 0;
    int scan_size = 0;
    residue_s[s] = 1;
    while(!is_done) {
        #ifdef HYBRID    
        if (scan_size > 0 || (*queue[queuei]).numElements() > graph.n/hybridsize) {
            (*queue[queuei]).clear();
            forward_type = 1;
            if (first_scan == 0) {
                first_scan = 1;
                int grain_size = 1;
                cilk_for(int i = 0;i < graph.n; i += grain_size) {
                    //#pragma simd
                    //for (int j = i; j < i+grain_size; j++) {
                    //    residue_push[j] = residue_s[j];
                    //    reserve_double[j] = reserve_s[j];
                    //}
                    residue_push[i] = residue_s[i];
                    reserve_double[i] = reserve_s[i];
                }
            }
            //printf("push num = %d pull num = %d\n", push_sum.get_value(), pull_sum.get_value());
            push_sum.set_value(0);
            push_sum.set_value(0);
            iter++;
            scan_size = sep_scan_pull(queue[queuei]);

            if (scan_size < graph.n/hybridsize) {
                forward_type = 1;
                scan_size = 0;
            }
            //printf("scanned size is %d\n", scan_size);
            //continue;
        } else {
            is_done = 1; 
        }
        if (forward_type == 1)
            continue;
        #endif

        (*queue[!queuei]).clear();
        process_bag(&((*queue[queuei]).get_reference()), queue[!queuei]);

        double source_residue_l = s_residue.get_value();
        s_residue.set_value(0.0);
        long long degree_s = graph.get_degree(source);
        if (residue_s[source] >= graph.rmax_scaled * degree_s) {
            atomic_add(residue_s[source], source_residue_l);
        } else {
            atomic_add(residue_s[source], source_residue_l);
            long long degree_s = graph.get_degree(source);
            if (residue_s[source] > graph.rmax_scaled * degree_s) {
                #ifdef USEPOUCH
                Pouch* pouch = new Pouch();
                Bag_reducer<Pouch*> *next = queue[!queuei];
                Bag<Pouch*>* bnext = &((*next).get_reference());
                put_in_pouch(pouch, source, bnext);
                (*bnext).insert(pouch);
                #else
                Bag_reducer<int> *next = queue[!queuei];
                Bag<int>* bnext = &((*next).get_reference());

                (*bnext).insert(source);
                #endif
            }
        }
        // prepare for next iteration
        queuei = !queuei;
        is_done = (*queue[queuei]).isEmpty();
        iter++;
        continue;
    }
    // Stop the timer now
    clock_gettime(CLOCK_REALTIME, &stoptime);
    #ifdef LOCALTIMER
    // local timer
    printf("%f milliseconds\n", local_elapsed/1000000);
    #endif 

    // get rsum
    cilk::reducer_opadd<double> residueSum(0.0);
    cilk_for(int i = 0 ; i < graph.dividingline; i ++) {
        //residueSum += residue_push[i];
        residueSum += residue_s[i];
    }
    rsum = residueSum.get_value();
    // set up reserve array for monte carlo
//    cilk_for( int i = 0; i < graph.n; i ++) {
//       reserve_s[i] = reserve_double[i];
//    }
    #ifdef DEBUGINFO
    printf("iter = %d\n", iter);
    // calculate rsum
    cilk::reducer_opadd<double> reserveSum(0.0);
    
    cilk_for(int i = 0 ; i < graph.n; i ++) {
        reserveSum += reserve_double[i];
    }
    printf("rsum is %f\n", residueSum.get_value());
    //return;
    double r_sum = residueSum.get_value();
    double u_sum = reserveSum.get_value();
    printf("residue sum is %.15f\n", residueSum.get_value());
    printf("reserve sum is %.15f\n", u_sum);
    
    printf("total sum is %.15f\n", r_sum + u_sum);
    // output results
    exit(0);
    #endif
}
#endif


#endif  //FORWARD_PULL
//===================== END OF MESSAGE LIST FORWARD PUSH =======================
//////////////////////////////end of forward push//////////////////////////////


////////////////////////////RANDOM WALK starts here////////////////////////////
inline static unsigned long lrand() {
    return rand();
    // return sfmt_genrand_uint32(&sfmtSeed);
}

inline static double drand(){
	return rand()*1.0f/RAND_MAX;
    // return sfmt_genrand_real1(&sfmtSeed);
}

#ifndef FP_RESIDUE_DISTRIBUTION
// multithread random number generator
unsigned int SEED=1;
inline static unsigned long lrand_thd(int core_id) {
    return rand_r(&SEED);
}

inline static double drand_thd(int core_id){
    return ((double)lrand_thd(core_id)/(double)INT_MAX);
}

int nodes_dividingline;
void init_random_walks(int s) {
    // compress the node array
    //if (nodes == NULL)
    //    nodes = new int[graph.n];
    //size_nodes = graph.dividingline;
    //printf("in init random walks source = %d\n", source);

    #ifdef DEBUGINFO
    printf("rsum is %f\n", rsum);
    #endif
    #ifdef MC_BINS
    if (double_incre == NULL) {
        double_incre = new double[graph.n];
    }
    memset(double_incre, 0, sizeof(double)*graph.n);
    #endif // MC_BINS
}
void perform_int_rw_update(Result* result) {

    double omega = Omega * rsum;
    cilk::reducer_opadd<double> residue_source(0.0);
    #ifdef DEBUGINFO
    cilk::reducer_opadd<double> increaseSum(0.0);
    printf("omega is %f\n", omega);
    cilk::reducer_opadd<long int> rwcount(0);
    #endif
    //cilk_for (int i = 0; i < nodes_dividingline; i ++) {
    cilk_for (int i = 0; i < graph.dividingline; i ++) {

        //int node = nodes[i];
        int node = i;
        //int node = graph.permute_mapping_r[i];
        double residue_i = residue_s[node];
        if (residue_i > 0) {
            unsigned long num_rw = floor(residue_i * Omega);
            //rwcount += num_rw;
            
            long long int start_idx = graph.dest_pointers[node];
            long long j;
            for (j = 0; j < num_rw; j ++) {
                int t = graph.destinations[start_idx + j];
                //reserve_s_int[t] += 1;           
                //reserve_s_l[t] += 1;
            }
            
           if (ceil(residue_i * Omega) > floor(residue_i * Omega)) {
                double alpha_i = (Omega*residue_i - num_rw) / Omega;
                int t = graph.destinations[start_idx+j];
                atomic_add(reserve_s[t], alpha_i);
            }
            #ifdef DEBUGINFO
            increaseSum += residue_i;
            #endif
        }
    }
    cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
    //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
        
        //int node = nodes[i];
        int node = i;
        double residue_i = residue_s[node];
            //atomic_add(reserve_incre[source], increase);
        residue_source += residue_i;
        #ifdef DEBUGINFO
        increaseSum += residue_i;
        #endif
    }
    atomic_add(reserve_s[source], residue_source.get_value());
    
    // merge results
    return;
    cilk_for (int i = 0 ; i < graph.n; i ++) {
        //reserve_final[i] = reserve_s[i] + reserve_incre[i] +reserve_s_int[i]/Omega;
    }
    #ifdef DEBUGINFO
    printf("increased reserve sum is %.15f\n", increaseSum.get_value());
    printf("random walks count %ld\n", rwcount.get_value());
//    long long maxCount = 0;
    for (int i = 0 ; i < NCORE; i ++) {
        printf("%lld ",wlcount[i]);
    }
    printf("\n");
//    printf("Max workload is %lld\n", maxCount);
    #endif
    #ifdef RACEINFO
    printf("Racecount is %d %d\n", race_num.get_value(), racecount.load());
    #endif
}

//========================NUMA AWARE MONTE CARLO=============================
inline void mc_on_a_node(int node, double omega, int numa_node, int core_id) {
    //double residue_i = residue_s[node].load();
    double residue_i = residue_s[node];
    if (residue_i > 0) {

        
        unsigned long num_rw = floor(residue_i * Omega);
            #ifndef INTEGER_RWUPDATE
        double alpha_i = residue_i * omega / num_rw;        
        double increase = alpha_i / omega;
            #endif
        //long long int start_idx = graph.dest_pointers[i];
        //long long int start_idx = graph.gmap_dest_pointers[i];
        long long int start_idx = graph.dest_pointers[node];
        for (unsigned long long j = 0; j < num_rw; j ++) {
            int t = graph.destinations[start_idx + j];
            #ifdef INTEGER_RWUPDATE
                #ifdef UINT8_T_RWUPDATE 
            int prevalue=Atomic_add(reserve_numa_u8[numa_node][t], 1);
            if (prevalue == 255) {
                atomic_add(reserve_s[t], ((double)256.0)/Omega);
            }
                #else
            if (t < graph.target_divideline) {
                atomic_add(reserve_numa_i32[numa_node][t], 1);
            } else {
                atomic_add(reserve_numa_u8[numa_node][t], 1);
            }
                #endif
            #else
            atomic_add(reserve_numa[numa_node][t], increase);
            #endif
        }
        if (ceil(residue_i * Omega) > floor(residue_i * Omega)) {
            double alpha_i = (Omega*residue_i - num_rw) / Omega;
            int t = graph.destinations[start_idx+num_rw];
            atomic_add(reserve_s[t], alpha_i);
        }
    }

}
int GRAINSIZE_MC;
void my_cilk_for_mc(int start, int end) {
    if ( end-start < GRAINSIZE_MC) {
        
        double omega = Omega * rsum;
        //int core_id = __cilkrts_get_worker_number(); 
        int core_id = sched_getcpu();
        int numa_node = numa_node_of_cpu( sched_getcpu());
        for (int i = start; i < end; i ++) {
            #ifdef RWPERMUTE
            int node = graph.permute_mapping_r[i];
            #else
            int node = i;
            #endif
            //int node = graph.dest_gmap[i];
            mc_on_a_node(node, omega, numa_node, core_id);
        }
    } else {
        int mid = (start + end) / 2;
        cilk_spawn my_cilk_for_mc(start, mid);
        my_cilk_for_mc(mid, end);

        cilk_sync;
    }
}
//=======================END OF NUMA AWARE MONTE CARLO=========================

//&&&&&&&&&&&&&&&&&&&&&&&& MESSAGE LIST MONTE CARLO &&&&&&&&&&&&&&&&&&&&&&&&&&&
// 1. alloc bin arrays for each numa nodes
// 2. each numa node update messages to its own bin array
// 3. after messages are all sent, move to update stage
// 4. in update stage, workload balance is a problem

int BIN_OFFSET = 10;
int BIN_SIZE = 1 << BIN_OFFSET;
class Message{
public:
    int id;
    double update;
};

// use an atomic::int to achieve parallel queue
// the queue is divided into segments (message arrays)
// two segments are created initially, the second is a buffer. 
// when the index reaches the buffer, create a new segment
// segments are stored in a vector (only the pointer to the first message)

class Bin{
public:
    const static int ARRAY_SIZE_OFFSET = 16;
    atomic<int> index;
    int ARRAY_SIZE;
    vector<Message*> message_lists;    // store Message array pointers
    Message* messages_cur;

    Bin() {
        index.store(0);
        ARRAY_SIZE = 1<<ARRAY_SIZE_OFFSET;           // TODO tune this size
        message_lists = vector<Message*>();
        Message* message_array = new Message[ARRAY_SIZE];
        message_lists.push_back(message_array);
        messages_cur = message_array;
    }
    
    void create_buffer() {
        Message* message_buffer = new Message[ARRAY_SIZE];
        message_lists.push_back(message_buffer);
    }
 // need an array list
};
Bin* bins = NULL;
int bin_number = 0;
void set_up_bins_on_numa_nodes() {
    if (bins == NULL) {
        bin_number = graph.n / BIN_SIZE + 1;
        printf("bins number is %d\n", bin_number);
        bins = new Bin[bin_number];
        printf("bins setup\n");
    } else {
        cilk_for(int i = 0; i < bin_number; i ++) {
            Bin& bin = bins[i];
            for (int j = 0; j < bin.message_lists.size(); j ++) {
                Message* messages = bin.message_lists[j];
                memset(messages, 0, sizeof(Message) * bin.ARRAY_SIZE);
            }
        }
    }
}

#ifdef MC_BINS
void monte_carlo_by_bins() {
    
    race_num.set_value(0);
    racecount = 0;
    double omega = Omega * rsum;
    cilk::reducer_opadd<double> residue_source(0.0);
    #ifdef DEBUGINFO
    cilk::reducer_opadd<double> increaseSum(0.0);
    cilk::reducer_opadd<long int> rwcount(0);
    #endif
    
    // SCATTER
    cilk_for (int i = 0; i < graph.dividingline; i ++) {

        #ifdef RWPERMUTE
        int node = graph.permute_mapping_r[i];
        #else
        int node = i;
        #endif
        double residue_i = residue_s[node].load();
        if (residue_i > 0) {

            unsigned long num_rw = ceil(residue_i * Omega);
            double alpha_i = residue_i * omega / num_rw;        
            double increase = alpha_i / omega;
             
            long long int start_idx = graph.dest_pointers[node];
            #ifdef RWNESTEDLOOP
            cilk_for (unsigned long long j = 0; j < num_rw; j ++) {
            #else
            for (unsigned long long j = 0; j < num_rw; j ++) {
            #endif
                // get the destination
                int t = graph.destinations[start_idx + j];
                // get the corresponding bin
                int bin_id = t >> BIN_OFFSET;
                Bin& bin = bins[bin_id];
                // get the message array
                int queue_index = bin.index.fetch_add(1, std::memory_order_relaxed);
                int array_id = queue_index >> bin.ARRAY_SIZE_OFFSET;
                Message* messages = bin.message_lists[array_id];
                int message_idx = queue_index&(bin.ARRAY_SIZE-1);
                if (message_idx == 0) {
                    // create a new message buffer
                    Message* msgs = new Message[bin.ARRAY_SIZE];
                    bin.message_lists.push_back(msgs);
                }
                Message& message = messages[message_idx];
                message.id = t;
                message.update = increase;

                //double_incre[t] += increase;
            }
            #ifdef DEBUGINFO
            int core_id = __cilkrts_get_worker_number(); 
            wlcount[core_id] += num_rw;
            rwcount = rwcount + num_rw;
            #endif
        }
    }
    // test bins
    
    // GATHER
    
    for (int i = 0; i < bin_number; i ++) {
        Bin& bin = bins[i];
        int index = bin.index;
        int array_num = index / bin.ARRAY_SIZE;

        for (int j = 0; j < array_num; j ++) {
            Message* msgs = bin.message_lists[j];
            
            for (int k = 0; k < bin.ARRAY_SIZE; k ++) {
                int t = msgs[k].id;
                double update = msgs[k].update;
                double_incre[t] += update;
                #ifdef DEBUGINFO
                increaseSum += update;
                #endif
            }
        }
        Message* msgs = bin.message_lists[array_num];
        int length = index&(bin.ARRAY_SIZE-1);
        for (int k = 0; k < length; k ++) {
            int t = msgs[k].id;
            double update = msgs[k].update;
            double_incre[t] += update;
            #ifdef DEBUGINFO
            increaseSum += update;
            #endif
        }
    }
    
    cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
    //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
        
        //int node = nodes[i];
        int node = i;
        double residue_i = residue_s[node];
        //atomic_add(reserve_incre[source], increase);
        residue_source += residue_i;
        #ifdef DEBUGINFO
        increaseSum += residue_i;
        #endif
    }
    atomic_add(reserve_s[source], residue_source.get_value());

    #ifdef NUMA_AWARE // TODO check if reserve sum is 1
    cilk_for (int i = 0; i < graph.n; i ++) {
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            reserve_final[i] += reserve_numa[j][i];
        }
        reserve_final[i] += reserve_s[i];
    }
    #endif // NUMA_AWARE
    #ifdef DEBUGINFO
    for (int i = 0 ; i < NCORE; i ++) {
        printf("%lld ",wlcount[i]);
    }
    double reserve_sum = 0;
    for (int i = 0; i < graph.n; i ++) {
        reserve_sum += reserve_final[i];
    }
    printf("\n");
    printf("total reserve sum is %.15f\n", reserve_sum);
    printf("omega is %f\n", omega);
    printf("random walks count %ld\n", rwcount.get_value());
    long long maxCount = 0;
    printf("Max workload is %lld\n", maxCount);
    #endif
    #ifdef RACEINFO
    printf("Racecount is %d %d\n", race_num.get_value(), racecount.load());
    #endif

}
#endif //MC_BINS
//=======================END OF MESSAGE LIST MONTE CARLO=======================
// random walk with index
void perform_random_walk_update(Result* result) {

    #ifdef DEBUGINFO
    printf("perform_random_walk_update start rsum is %f\n", rsum);
    #endif
    #ifdef RACEINFO
    race_num.set_value(0);
    racecount = 0;
    #endif
    double omega = Omega * rsum;
    //cilk::reducer_opadd<double> residue_source(0.0);
    //if (graph.n != graph.dividingline) {
        cilk::reducer< cilk::op_add<double> > residue_source(0.0);
    //}
    #ifdef DEBUGINFO
    cilk::reducer_opadd<double> increaseSum(0.0);
    cilk::reducer_opadd<long int> rwcount(0);
    #endif
    //cilk_for (int i = 0; i < nodes_dividingline; i ++) {
#ifdef NUMA_AWARE
    #ifdef TWO_PHASE_RWNUMA
    int idx=0;
    cilk_for (; idx < graph.rw_dvdngline; idx += GRAINSIZE_MC) {
        int numa_node = sched_getcpu() % NCOPY;
        for (int j = idx; j < idx + GRAINSIZE_MC && j < graph.dividingline; j ++) {
            int node = j;
            double residue_i = residue_s[node];
            if (residue_i > 0) {
                unsigned long num_rw = ceil(residue_i * Omega);
                long long int start_idx = graph.dest_pointers[node];
                for (unsigned long long j = 0; j < num_rw; j ++) {
                    int t = graph.destinations[start_idx + j];
                    if (t < graph.target_divideline) {
                      reserve_numa_i32[numa_node][t] += 1;
                    } else {
                      reserve_numa_u8[numa_node][t<<CACHEPAD_U8] += 1;
                    }
                }
            }
        }
    }
    cilk_for (; idx < graph.dividingline; idx += GRAINSIZE_MC) {
        int numa_node = sched_getcpu() % NCOPY;
        for (int j = idx; j < idx + GRAINSIZE_MC && j < graph.dividingline; j ++) {
            int node = j;
            double residue_i = residue_s[node];
            if (residue_i > 0) {
                unsigned long num_rw = ceil(residue_i * Omega);
                long long int start_idx = graph.dest_pointers[node];
                for (unsigned long long j = 0; j < num_rw; j ++) {
                    int t = graph.destinations[start_idx + j];
                    reserve_numa_u8[numa_node][t<<CACHEPAD_U8] += 1;
                }
            }
        }
    }
    #else   // else TWO_PHASE_RWNUMA
    //my_cilk_for_mc(0, graph.dividingline);
    //vector<vector<int>> worker_map;
    //int nworker = __cilkrts_get_nworkers();
    //worker_map = vector<vector<int>> (80, vector<int>());
    #pragma cilk grainsize=1
    cilk_for(int i = 0; i < graph.dividingline ; i += GRAINSIZE_MC) {
    //for(int i = 0; i < graph.dividingline; i += GRAINSIZE_MC) {
        //int numa_node = sched_getcpu()%NCOPY;
        int core_id = sched_getcpu();
        //int core_id = 0;
        //int numa_node = numa_node_of_cpu( core_id);
        int numa_node = core_id % NCOPY;
        //int core_id = __cilkrts_get_worker_number(); 
        /*
        int core_id = __cilkrts_get_worker_number(); 
        worker_map[core_id].push_back(sched_getcpu());
        */
        int start = i;
        int end = i + GRAINSIZE_MC;
        if (end > graph.dividingline) end = graph.dividingline;
        for (int j = start; j < end; j ++) {
            mc_on_a_node(j, omega, numa_node, core_id);
        }   
    }
    #endif  // endif TWO_PHASE_RWNUMA
    /* 
    for (int i = 0; i < nworker; i ++) {
        int cpuid = worker_map[i][0];
        int numa_id = numa_node_of_cpu(cpuid);
        printf("worker %d is attached to # %d on numa node %d", i, cpuid,numa_id );
        for (int j = 0; j < worker_map[i].size(); j ++) {
            if ( cpuid != worker_map[i][j])
                printf(" %d", worker_map[i][j]);
        }
        printf("\n");
    }*/
#else   // NUMA_AWARE
    cilk_for (int i = 0; i < graph.dividingline; i ++) {

        #ifdef RWPERMUTE
        int node = graph.permute_mapping_r[i];
        #else
        int node = i;
        #endif
        //int node = graph.dest_gmap[i];
        double residue_i = residue_s[i].load();
        //double residue_i = residue_s[node];
        if (residue_i > 0) {

            unsigned long num_rw = floor(residue_i * Omega);
            double alpha_i = residue_i * omega / num_rw;        
            double increase = alpha_i / omega;
             
            //long long int start_idx = graph.dest_pointers[i];
            //long long int start_idx = graph.gmap_dest_pointers[i];
            long long int start_idx = graph.dest_pointers[node];
            for (unsigned long long j = 0; j < num_rw; j ++) {
                int t = graph.destinations[start_idx + j];
                //int t = graph.gmap_destinations[start_idx + j];
                //atomic_add(reserve_incre[t<<CACHEPADSIZE_D], increase);
                #ifdef INTEGER_RWUPDATE
                    #ifdef UINT8_T_RWUPDATE 
                int prevalue = Atomic_add(reserve_s_u8[t], 1);
                if (prevalue == 255) {
                    atomic_add(reserve_s[t], (double)256.0/Omega);
                }
                    #else
                if (t < graph.target_divideline) {
                    atomic_add(reserve_numa_i32[numa_node][t], 1);
                } else {
                    atomic_add(reserve_numa_u8[numa_node][t], 1);
                }
                    #endif
                #else
                atomic_add(reserve_s[t], increase);
                #endif
                #ifdef DEBUGINFO
                increaseSum += increase;
                #endif
            }
            #ifdef DEBUGINFO
            int core_id = __cilkrts_get_worker_number(); 
            wlcount[core_id] += num_rw;
            rwcount = rwcount + num_rw;
            #endif
            if (ceil(residue_i * Omega) > floor(residue_i * Omega)) {
                double alpha_i = (Omega*residue_i - num_rw) / Omega;
                int t = graph.destinations[start_idx+num_rw];
                atomic_add(reserve_s[t], alpha_i);
            }
        }
    }
#endif // NUMA_AWARE
    if (graph.n != graph.dividingline) {
        cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
        //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
            
            //int node = nodes[i];
            int node = i;
            double residue_i = residue_s[node];
            //atomic_add(reserve_incre[source], increase);
            *residue_source += residue_i;
            #ifdef DEBUGINFO
            increaseSum += residue_i;
            #endif
        }
        atomic_add(reserve_s[source], residue_source.get_value());
    }
    #ifdef NUMA_AWARE
    #ifdef INTEGER_RWUPDATE
    
        #ifdef UINT8_T_RWUPDATE 
    #pragma simd
    cilk_for (int i = 0; i < graph.n; i ++) {
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            ((double*)reserve_s)[i] += ((double)reserve_numa_u8[j][i])/Omega;
        }
    }
        #else
    #pragma simd
    cilk_for (int i = 0; i < graph.n; i ++) {
        int node = graph.target_map[i];
        int type = node < graph.target_divideline? 1:0;
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            reserve_final[i] += type ? ((int*)reserve_numa_i32[j])[node]/graph.Omega :
                                       ((uint8_t*)reserve_numa_u8[j])[node]/graph.Omega ;
        }
        reserve_final[i] += ((double*)reserve_s)[i];
    }
        #endif
    
    #else // INTEGER_RWUPDATE
    cilk::reducer_opadd<double> incrSum(0.0);
    #pragma simd
    cilk_for (int i = 0; i < graph.n; i ++) {
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            reserve_final[i] += ((double*)reserve_numa[j])[i];
        }
        incrSum += reserve_final[i];
    }
    cout << incrSum.get_value() << endl;;
    #pragma simd
    cilk_for (int i = 0; i < graph.n; i ++) {
        reserve_final[i] += ((double*)reserve_s)[i];
    }
    //printf("result merged\n");
    incrSum.set_value( 0);
    cilk_for (int i = 0; i < graph.n; i ++) {
        incrSum += ((double*)reserve_final)[i];
    }
    cout << incrSum.get_value() << endl;
    #endif

    #else // NUMA_AWARE
    cilk_for (int i = 0; i < graph.n; i ++) {
        reserve_final[i] += reserve_s[i];
    }
    #endif

    #ifdef CHECK_TOPK
    source = source_origin;
    map<int, double> exact_map;
    int size_e = exact_topk_pprs[source].size();
    int k = 500;
    //size_e = k;
    for(int i=0; i<size_e; i++){
        pair<int ,double>& p = exact_topk_pprs[source][i];
        if(p.second>0){
            exact_map.insert(p);
        }
    }
    vector<pair<int, double> > topk_pair;
    for (int i = 0; i < graph.n; i ++) {
        #ifdef UINT8_T_RWUPDATE
        topk_pair.push_back(make_pair(i, reserve_s[i].load())); 
        #else
        topk_pair.push_back(make_pair(i, reserve_final[i])); 
        #endif
    }
    sort(topk_pair.begin(), topk_pair.end(), pair_sorter_large_first);
    //cout << topk_pair[0].second << "    " << exact_topk_pprs[source][0].second << endl;

    double precision = 0.0;
    for (int i = 0; i < k; i ++) {
        if(exact_map.find(graph.sep_reverse_map[topk_pair[i].first])!=exact_map.end()){
            precision++;
        }
        //cout << "NO." << i << " pred:" << graph.sep_reverse_map[topk_pair[i].first] << ", " << topk_pair[i].second << endl << "\t exact:" << exact_topk_pprs[source][i].first << ", " << exact_topk_pprs[source][i].second << endl;
    }    
    //cout << precision / k << endl;
//    return;
//
//    FILE *topk= fopen("topk_p_4.txt", "w");
    double error_max = 0;
    double relative_max = 0;
    //topk_pair.resize(k);
    for(int i = 0; i < k; i ++) {
        int node = topk_pair[i].first;
        if (exact_map.find(graph.sep_reverse_map[node]) != exact_map.end()){
            double exact_ppr = exact_map[graph.sep_reverse_map[node]];
            if (exact_ppr < 1.0/graph.n) continue;
            double error = abs(exact_ppr - topk_pair[i].second);
            if (error > error_max) error_max = error;
            double relative = error/exact_ppr;
            if (relative > relative_max) relative_max = relative;
        }
    }
    cout << error_max << "      " << relative_max << endl;

    #endif

    return;
    
//    for (int i = 0 ; i < NCORE; i ++) {
//        printf("%lld ",wlcount[i]);
//    }
    double reserve_sum = 0;
    for (int i = 0; i < graph.n; i ++) {
        #ifdef UINT8_T_RWUPDATE
        reserve_sum += reserve_s[i];
        #else
        reserve_sum += reserve_final[i];
        #endif
    }
    printf("\n");
    printf("total reserve sum is %.15f\n", reserve_sum);
    #ifdef DEBUGINFO
    printf("omega is %f\n", omega);
    printf("random walks count %ld\n", rwcount.get_value());
    long long maxCount = 0;
    printf("Max workload is %lld\n", maxCount);
    #endif
    #ifdef RACEINFO
    printf("Racecount is %d %d\n", race_num.get_value(), racecount.load());
    #endif
}
void calculated_workload_random_walk_update(Result* result) {
    race_num.set_value(0);
    racecount = 0;
    double omega = Omega * rsum;
    cilk::reducer_opadd<double> residue_source(0.0);
    #ifdef DEBUGINFO
    cilk::reducer_opadd<double> increaseSum(0.0); 
    cilk::reducer_opadd<long int> rwcount(0);
    #endif
    //cilk_for (int i = 0; i < nodes_dividingline; i ++) {
    // calculate workload
    int core_num = __cilkrts_get_nworkers();
    double residue_thres = rsum / core_num;
    int* start_indices = new int[core_num+1];
    int idx=0;
    start_indices[idx++] = 0;
    double accum = 0;
    for (int i = 0 ; i < dividingline; i ++) {
        if (accum >= residue_thres) {
            accum = 0;
            start_indices[idx] = i+1;
            idx ++;
            if (idx == core_num) {
                start_indices[idx] = dividingline;
                break;
            }
        }
        accum += residue_s[i];
    }

    cilk_for (int k = 0; k < core_num; k ++) {
        for (int i = start_indices[k]; i < start_indices[k+1]; i++) {
            int node = i;
            //int node = graph.permute_mapping_r[i];
            double residue_i = residue_s[node];
            if (residue_i > 0) {

                unsigned long num_rw = ceil(residue_i * Omega);
                double alpha_i = residue_i * omega / num_rw;
                
                double increase = alpha_i / omega;
                 
                long long int start_idx = graph.dest_pointers[i];
                //long long int start_idx = graph.dest_pointers[node];
                for (unsigned long long j = 0; j < num_rw; j ++) {
                    int t = graph.destinations[start_idx + j];
                    atomic_add(reserve_s[t], increase);
                    //reserve_s_ll[t] += 1;
                    #ifdef DEBUGINFO
                    increaseSum += increase;
                    #endif
                }
                #ifdef DEBUGINFO
                int core_id = __cilkrts_get_worker_number(); 
                wlcount[core_id] += num_rw;
                rwcount = rwcount + num_rw;
                #endif
            }
        }
    }
    cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
    //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
        
        //int node = nodes[i];
        int node = i;
        double residue_i = residue_s[node];
        //atomic_add(reserve_incre[source], increase);
        residue_source += residue_i;
        #ifdef DEBUGINFO
        increaseSum += residue_i;
        #endif
    }
    atomic_add(reserve_s[source], residue_source.get_value());

    #ifdef DEBUGINFO
    for (int i = 0 ; i < NCORE; i ++) {
        printf("%lld ",wlcount[i]);
    }
    double reserve_sum = 0;
    for (int i = 0; i < graph.n; i ++) {
        reserve_sum += reserve_final[i];
    }
    printf("\n");
    printf("total reserve sum is %.15f\n", reserve_sum);
    printf("omega is %f\n", omega);
    printf("random walks count %ld\n", rwcount.get_value());
    long long maxCount = 0;
    printf("Max workload is %lld\n", maxCount);
    #endif
    #ifdef RACEINFO
    printf("Racecount is %d %d\n", race_num.get_value(), racecount.load());
    #endif
}
#ifdef RANDOMWALK
// online random walks
//cilkpub::DotMix global_rng(234909128);
inline int random_walk(int start){
    int cur = start;
    unsigned long k;
    while (true) {
        #ifdef USE_DOTMIX 
	    double random_num = (double)global_rng.get()/UINTMAX_MAX;
        #else
        double random_num = drand();
        #endif
        if ( random_num < graph.alpha) {
            return cur;
        }
		int degree = graph.get_degree(cur);
        if (degree){
            long long int start_idx = graph.nodepointer[0][cur];
            #ifdef USE_DOTMIX
            k = global_rng.get() % degree;
            #else
            k = rand() % degree;
            #endif
            cur = graph.edges[0][start_idx + k];
        }
        else{
            cur = start;
        }
    }
}
#else
inline int random_walk(int start) {
    return 0;
}
#endif
void online_random_walk_update(Result* result) {
    race_num.set_value(0);
    racecount = 0;
    double omega = Omega * rsum;
    cilk::reducer_opadd<double> residue_source(0.0);
    #ifdef DEBUGINFO
    cilk::reducer_opadd<double> increaseSum(0.0);
    cilk::reducer_opadd<long int> rwcount(0);
    #endif
    

    //cilk_for (int i = 0; i < nodes_dividingline; i ++) {
    cilk_for (int i = 0; i < graph.dividingline; i ++) {

        int node = i;
        //int node = graph.permute_mapping_r[i];
        //int node = graph.dest_gmap[i];
        double residue_i = residue_s[i].load();
        //double residue_i = residue_s[node];
        if (residue_i > 0) {

            unsigned long num_rw = ceil(residue_i * Omega);
            double alpha_i = residue_i * omega / num_rw;        
            double increase = alpha_i / omega;
             
            //long long int start_idx = graph.dest_pointers[i];
            //long long int start_idx = graph.gmap_dest_pointers[i];
            long long int start_idx = graph.dest_pointers[node];
            for (unsigned long long j = 0; j < num_rw; j ++) {
                //int t = graph.destinations[start_idx + j];
                int t = random_walk(i);
                //atomic_add(reserve_incre[t<<CACHEPADSIZE_D], increase);
                atomic_add(reserve_s[t], increase);
                //reserve_s_ll[t] += 1;
                //reserve_s_int[t<<CACHEPADSIZE_I] += 1;
                //double_incre[t] += increase;
                #ifdef DEBUGINFO
                increaseSum += increase;
                #endif
            }
            #ifdef DEBUGINFO
            int core_id = __cilkrts_get_worker_number(); 
            wlcount[core_id] += num_rw;
            rwcount = rwcount + num_rw;
            #endif
        }
    }
    cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
    //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
        
        //int node = nodes[i];
        int node = i;
        double residue_i = residue_s[node];
        //atomic_add(reserve_incre[source], increase);
        residue_source += residue_i;
        #ifdef DEBUGINFO
        increaseSum += residue_i;
        #endif
    }
    atomic_add(reserve_s[source], residue_source.get_value());

    #ifdef DEBUGINFO
    for (int i = 0 ; i < NCORE; i ++) {
        printf("%lld ",wlcount[i]);
    }
    printf("\n");
    printf("increased reserve sum is %.15f\n", increaseSum.get_value());
    printf("omega is %f\n", omega);
    printf("random walks count %ld\n", rwcount.get_value());
    long long maxCount = 0;
    printf("Max workload is %lld\n", maxCount);
    #endif
    #ifdef RACEINFO
    printf("Racecount is %d %d\n", race_num.get_value(), racecount.load());
    #endif
}
//////////////////////////////end of random walk///////////////////////////////
#endif // pair with the start of random walk

#ifndef FP_RESIDUE_DISTRIBUTION
// split the nodes with largest residue values into smaller pieces
/*
void redistribute(double maxMag) {
    cilk::reducer<cilk::op_list_append<int> > nodes_reducer;
    cilk_for(int i = 0; i < graph.n; i ++) {
        double residue = residue_s[i];
        if (residue > maxMag/100) {
            nodes_reducer->push_back(i);
        }
    }
    
    list<int> nodes_list = nodes_reducer.get_value();
    printf("node list size is %lu\n", nodes_list.size());
    int scale = 32*2;
    for (list<int>::iterator i = nodes_list.begin(); i!= nodes_list.end();i++){
        int node = *i;
        double residue = residue_s[node];
        
        //printf("residue is %f, out degree is %ld\n", residue, 
        //                                        graph.g[node].size());
        int scale_i = 1;
        if (residue > maxMag) scale_i = scale*128;
        else if (residue > maxMag/10) scale_i = scale*16;
        else scale_i = scale;
        residue /= scale_i;
        //residue_s[node].store(residue);
        // SPIN
        spinlock(node);
        residue_s[node] = residue;
        unlock(node);
        // this node is already in nodes so we need to append scale-1 times
        for (int j = 0; j < scale_i-1; j ++) {
            nodes[size_nodes++] = node;
        }
    }
}*/
#endif
//void generate_locations(int* locs, int size) {
//    int iter = 0;
//    int index = 0;
//    locs[index++] = 0;
//    for (int i = 0; i < size; i ++) {
//        int n = pow(2,i);
//        int denominator = pow(2,i+1);
//        for (int j = 0; j < n; j ++) {
//            int numerator = 2 * j + 1;
//            locs[index] = (int)(size_nodes * numerator / denominator) + 2;
//            //printf("%d %d\n", numerator, denominator);
//            index ++;
//            if (index == size) return;
//        }
//    }
//}
#ifndef FP_RESIDUE_DISTRIBUTION
// permute the nodes array
/*
void merge_sort(int* array, int begin, int end) {
    if (begin < end) {

        int mid = (begin + end) / 2;
        merge_sort(array, begin, mid);
        merge_sort(array, mid+1, end);

        int* result = new int[end - begin+1];
        int result_index = 0;
        int li = begin;
        int ri = mid+1;
        while (li <= mid && ri <= end) {
            int nodel = nodes[array[li]];
            int noder = nodes[array[ri]];
            if (residue_s[nodel] > residue_s[noder]) {
                result[result_index++] = array[li++];
            } else {
                result[result_index++] = array[ri++];
            }
        }
        if (li > mid) {
            while (ri <= end) {
                result[result_index++] = array[ri++];
            }
        } else {
            while (li <= mid) {
                result[result_index++] = array[li++];
            }
        }
        if (li != mid+1 || ri != end+1 || result_index != end-begin+1) {
            printf("merge sort is wrong!\n");
        }
        for (int i = 0; i < result_index; i ++) {
            array[begin+i] = result[i];
        }
    }
}*/
#endif
int *in_locations;
int *permute_map;
// relocate the nodes with largest residue values to critical locations
// Strategy: 

#ifndef FP_RESIDUE_DISTRIBUTION
/*
void permute_nodes(double maxMag) {

    cilk::reducer<cilk::op_list_append<int> > nodes_reducer;
    cilk_for(int i = 0; i < size_nodes; i ++) {
        double residue = residue_s[nodes[i]];
        if (residue > maxMag/50) {
            nodes_reducer->push_back(i);
        }
    }

    list<int> nodes_list = nodes_reducer.get_value();
    int list_size = nodes_list.size();
    int* locations = new int[list_size];
    int* nodes_array = new int[list_size];
    int j = 0;
    for (list<int>::iterator i = nodes_list.begin(); i!= nodes_list.end();i++, j ++){
        nodes_array[j] = *i;
    }
    merge_sort(nodes_array, 0, list_size);
    generate_locations(locations, list_size);
    // adjust locations
    for (int i = 0; i < list_size; i ++) {
        
    }

    printf("node list size is %lu source is%d\n", nodes_list.size(), source);
    int index = 0;
     
    int core_num = __cilkrts_get_nworkers();
    //int core_num = 8;
    int priority_limit = core_num * 300;
    
    for (int i = 0; i < priority_limit && i < list_size; i ++) {
        int niter = i / core_num;
        locations[i] = locations[i%core_num] + niter;
    }
    // The following for loop used to be above the for loop above
    in_locations = new int[size_nodes];
    for (int i = 0; i < list_size; i ++) {
        double residue = residue_s[nodes[locations[i]]];
        if (residue > maxMag/50) {
            in_locations[locations[i]] = 1;
        }
    }

    for (int i = 0; i < list_size && i < priority_limit; i ++) {
        
        int idx = nodes_array[i];

        int location = locations[index++];
        if (in_locations[location] == 1) {
            
            for (int j = i+1; j < list_size; j ++) {
                if (nodes_array[j] == location) {
                    nodes_array[j] = idx;
                }
            }
        } 
        int tmp = nodes[idx];
        nodes[idx] = nodes[location];
        nodes[location] = tmp;
    }
}*/
inline double residue_distribution_test() {
    int SPAN = 30;
    long* count = new long[SPAN];
    for (int i = 0; i < SPAN; i ++) {
        count[i] = 0;
    }
    int rcount = 0;

    double minmagnitude= 1.0/pow(10.0,SPAN);
    double maxMag = minmagnitude;
    printf("Residue Distribution: \n");
    for (int i = 0; i < graph.n; i ++) {
        double residue = residue_s[i];
        if (residue > 0) rcount ++;
        double bound = minmagnitude;
        int idx = 0;
        while (bound < 1) {
            if (residue < bound) {
                count[idx] ++;
                break;
            } else {
                idx++;
                bound *= 10;
                if (maxMag < bound) 
                    maxMag = bound;
            }
        }
    }
    for (int i = 0; i < SPAN; i ++) {
        printf("%ld ", count[i]);
    }
    printf("\n total = %d\n", rcount);
    return maxMag/10;
}


void print_timer(char*ss) {
    double elapsed = 0;
    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);

    printf("%s processed in %.3f miliseconds using %d workers.\n",
           ss, (elapsed)/1000000, __cilkrts_get_nworkers());
    printf("\n");
}
#endif //FP_RESIDUE_DISTRIBUTION
void start_timer() {

    clock_gettime(CLOCK_REALTIME, &start);
}
#ifdef FP_RESIDUE_DISTRIBUTION
void test1(int s, Result* result, atomic<double>* residue_s,
                                  atomic<double>* reserve_s,
                                  double* push_counts) {
#else
void test1(int s, Result* result) {
#endif
    #ifdef DEBUGINFO

    int nCore = MAX_CORE_N;
    if (wlcount == NULL)
        wlcount = new long long [nCore];
    for (int i = 0; i < nCore; i ++) {
        wlcount[i] = 0;
    }
    if (wlcountotal == NULL)
        wlcountotal = new long long [nCore];
    for (int i = 0; i < nCore; i ++) {
        wlcountotal[i] = 0;
    }
    #endif
    // prepare for random walk indexing
    // TEST START
    int nworkers = __cilkrts_get_nworkers();
    // maximize init efficiency
    #ifndef FP_RESIDUE_DISTRIBUTION
    setWorkers(80);
    #endif

    start_timer();
    // initialize residue 
    #ifndef FP_RESIDUE_DISTRIBUTION
    init_residue_s(s);
    #endif
    // initialize s reserve vector
    init_s_reserve(s);
    //print_timer("init reserve value");
//    struct timespec start, finish;
    double elapsed;

    #ifndef FP_RESIDUE_DISTRIBUTION
    setWorkers(nworkers);
    #endif

    clock_gettime(CLOCK_REALTIME, &start);
    
    #ifdef LOCAL_UPDATE_FP
    forward_push_local_update(s, finish);
    #else
        #ifdef PULL
        forward_push_pull_update(s, finish);
        #else
            #ifdef FP_RESIDUE_DISTRIBUTION
            forward_push(s, finish, residue_s, reserve_s, push_counts);
            #else
            forward_push(s, finish); 
            #endif
        #endif
    #endif
    //printf("after forward push source = %d\n", source);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    #ifdef DEBUGINFO
    printf("Forward push calculated in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());
    printf("\n");
    #endif
    result->add_forward(elapsed/1000000);

    result->add_race_stat(race_num.get_value());
    // END OF FORWARD PUSH 
    //#ifdef FPWORKLOADCOUNT
    //return;
    //__cilkscreen_disable_instrumentation();
    #ifdef FP_RESIDUE_DISTRIBUTION
    // RANDOM walk
    // return;
    return;
    
    #else // FP_RESIDUE_DISTRIBUTION

    init_random_walks(s);
    source = s;
    clock_gettime(CLOCK_REALTIME, &start);
    #ifdef RANDOMWALK
    if (rwmode == 0) {
        #ifdef RWINDEX
            #ifdef BINS_FOR_MC
            monte_carlo_by_bins(); 
            #else
            #ifndef FP_RESIDUE_DISTRIBUTION
            perform_random_walk_update(result);
            #endif
            #endif
        #else
        online_random_walk_update(result);
        #endif
    } else if (rwmode == 1){
        #ifndef FP_RESIDUE_DISTRIBUTION
        perform_int_rw_update(result);
        #endif
    } else {
        #ifndef FP_RESIDUE_DISTRIBUTION
        calculated_workload_random_walk_update(result);
        #endif
    }
    #endif // RANDOMWALK
    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);

    #ifdef DEBUGINFO
    printf("Random walk calculated in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    printf("\n");
    #endif
    result->add_random(elapsed/1000000);

    #ifdef LOSSCHECK
    // incre sum
    double sum = 0;
    for (int i = 0; i < dividingline; i ++) {
        sum += double_incre[i];
    }
    printf("%f\n", sum/rsum);
    #endif
    #ifdef DEBUGINFO
    printf("------------------------END OF TEST-------------------------\n");
    #endif

    #endif // FP_RESIDUE_DISTRIBUTION
}

// test performance of parallel random number generator
void testParallelRand(int n) {

    printf("------------START OF PARALLEL RAND TEST---------------\n\n");
    struct timespec start, finish;
    double elapsed;

    int workerN = n;
    long reps = 100000000;
    setWorkers(workerN);
    clock_gettime(CLOCK_REALTIME, &start);

    cilk_for(int j = 0; j < workerN; j ++) {
        int core_id = __cilkrts_get_worker_number(); 
        for (int i = 0; i < reps; i ++) {
            //drand_thd(core_id);
            //drand();
            global_rng.get();
        } 
    }

    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Parallel drand executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());

    setWorkers(1);
    clock_gettime(CLOCK_REALTIME, &start);
        for (int i = 0; i < workerN * reps; i ++) {
            //drand();
            global_rng.get();
        }

    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("serial drand executed in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    printf("------------END OF PARALLEL RAND TEST---------------\n\n");
    
}

void testAtomicInt() {

    printf("------------------BEGIN OF ATOMIC ADD TEST-----------------\n\n");
    struct timespec start, finish;
    double elapsed;

    int testsize = 1000000000;
    clock_gettime(CLOCK_REALTIME, &start);
    int* ordin = new int[testsize];
    for (int i = 0; i < testsize; i ++) {
        ordin[i] += 1; 
    } 

    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Add executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());

    int corenum = __cilkrts_get_nworkers();
    int size = 1000000;
    int n = corenum * size;
    clock_gettime(CLOCK_REALTIME, &start);
    int* testintarr = new int[n];
    // RANDOM walk
    atomic<int>* ad = new atomic<int>[testsize];
    for (int i = 0; i < testsize; i ++) {
        ad[i].fetch_add(1, std::memory_order_relaxed);
    } 
    clock_gettime(CLOCK_REALTIME, &finish);


    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Atomic add executed in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    // atomic int
    clock_gettime(CLOCK_REALTIME, &start);
    // RANDOM walk
    memset(testintarr, 0 , sizeof(double)*n); 
    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Atomic int add executed in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    printf("------------------END OF ATOMIC ADD TEST--------------------\n\n");
    
}
void testSomething() {

    printf("------------------BEGIN OF ATOMIC ADD TEST-----------------\n\n");
    struct timespec start, finish;
    double elapsed;
    
    setWorkers(40);
    int testsize = 1000000;
    clock_gettime(CLOCK_REALTIME, &start);

    //do something here
    cilk_for (int i = 0; i < testsize; i ++) {
        __cilkrts_get_worker_number();
        //sched_getcpu();
    }

    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Add executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());
}
// test performance difference between atomic add and regular add
void testAtomic() {

    printf("------------------BEGIN OF ATOMIC ADD TEST-----------------\n\n");
    struct timespec start, finish;
    double elapsed;

    int testsize = 1000000000;
    clock_gettime(CLOCK_REALTIME, &start);
    double* ordin = new double[testsize];
    for (int i = 0; i < testsize; i ++) {
        ordin[i] += 1; 
    } 

    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Add executed in %.3f miliseconds using %d workers.\n",
           elapsed/1000000, __cilkrts_get_nworkers());

    int corenum = __cilkrts_get_nworkers();
    int size = 1000000;
    int n = corenum * size;
    clock_gettime(CLOCK_REALTIME, &start);
    double* testdoublearr = new double[n];
    // RANDOM walk
    atomic<double>* ad = new atomic<double>[testsize];
    for (int i = 0; i < testsize; i ++) {
        atomic_add(ad[i], 1);
    } 
    clock_gettime(CLOCK_REALTIME, &finish);


    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Atomic add executed in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    // atomic int
    clock_gettime(CLOCK_REALTIME, &start);
    // RANDOM walk
    memset(testdoublearr, 0 , sizeof(double)*n); 
    clock_gettime(CLOCK_REALTIME, &finish);

    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    printf("Atomic int add executed in %.3f miliseconds using %d workers.\n",
           (elapsed)/1000000, __cilkrts_get_nworkers());
    printf("------------------END OF ATOMIC ADD TEST--------------------\n\n");
    
}
#ifdef FP_RESIDUE_DISTRIBUTION
void scaletest(int n, int source, Result* result, atomic<double>* residue_s,
                                                  atomic<double>* reserve_s,
                                                  double* push_counts) {
#else
void scaletest(int n, int source, Result* result) {
#endif
    #ifdef FP_RESIDUE_DISTRIBUTION
    test1(source, result, residue_s, reserve_s, push_counts);
    #else
    setWorkers(n);
    if (run_topk) {
        topk_query_with_bound(source, 500, result);
    } else
        test1(source, result);
    #endif
    return;
}

void load_ss_query(vector<int>& queries, string loc){
    string filename = loc+"/ssquery.txt";
    ifstream queryfile(filename);
    int v;
    while(queryfile>>v){
        queries.push_back(v);
    }
}

bool pair_sorter_large_first(pair<int, double> p1, pair<int, double> p2) {
    if (p1.second > p2.second)
        return true;
    else if (p1.second < p2.second)
        return false;
    else if (p1.first < p2.first)
        return true;
    else
        return false;
}
void write_forwardpush_stats(MultiSources* msource, string graph_path) {

    vector<pair<int, double> > rc_pair; // residue count pairs
    for (int i = 0; i< graph.n; i ++) {
        rc_pair.push_back(make_pair(i, msource->stats[i].load()));
    }
    // printf some results
    for (int i = 0; i < 20; i ++) {
//        printf("%f ", rc_pair[i].second);
    }
    // sort
    sort( rc_pair.begin(), rc_pair.end(), pair_sorter_large_first);
    // print some results
    printf("Top 20 hits\n");
    for (int i = 0; i < 20; i ++) {
        printf("node = %d, counts = %f, out degree = %ld, in degree = %ld\n",
                rc_pair[i].first, rc_pair[i].second,
                graph.get_degree(rc_pair[i].first),
                graph.get_in_degree(rc_pair[i].first));
                //graph.gr[rc_pair[i].first].size());
    }
    // create a pagerank map
    //map<int, double> pagerank_map;
    //for (int i = 0; i < graph.n; i ++) {
    //    pagerank_map.insert(graph.pagerank_pair[i]);
    //}

    FILE *stats = fopen((graph_path+"/stats.txt").c_str(), "w");
    fprintf(stats, "#%s\n", graph_path.c_str());
    for (int i = 0; i < graph.n; i++) {
        int node = rc_pair[i].first;
        fprintf(stats, "%d %d %f %ld %ld\n", i,  node, rc_pair[i].second, 
                            //pagerank_map[rc_pair[i].first],
                            graph.get_degree(node), graph.get_in_degree(node));
    }
    fprintf(stats, "e\n");
    fflush(stats);
}
void write_forwardpush_push_counts(MultiSources* msource, string graph_path) {

    vector<pair<int, double> > rc_pair; // residue count pairs
    for (int i = 0; i< graph.n; i ++) {
        rc_pair.push_back(make_pair(i, msource->final_counts[i].load()));
    }
    // printf some results
    for (int i = 0; i < 20; i ++) {
//        printf("%f ", rc_pair[i].second);
    }
    // sort
    sort( rc_pair.begin(), rc_pair.end(), pair_sorter_large_first);
    // print some results
    printf("Top 20 pushes\n");
    for (int i = 0; i < 20; i ++) {
        printf("node = %d, counts = %f, out degree = %ld, in degree = %ld\n",
                rc_pair[i].first, rc_pair[i].second,
                graph.get_degree(rc_pair[i].first),
                graph.get_in_degree(rc_pair[i].first));
                //graph.gr[rc_pair[i].first].size());
    }
    // create a pagerank map
    //map<int, double> pagerank_map;
    //for (int i = 0; i < graph.n; i ++) {
    //    pagerank_map.insert(graph.pagerank_pair[i]);
    //}

    FILE *stats = fopen((graph_path+"/hot_nodes_stats.txt").c_str(), "w");
    fprintf(stats, "#%s\n", graph_path.c_str());
    for (int i = 0; i < graph.n; i++) {
        int node = rc_pair[i].first;
        fprintf(stats, "%d %d %f %ld %ld\n", i,  node, rc_pair[i].second, 
                            //pagerank_map[rc_pair[i].first],
                            graph.get_degree(node), graph.get_in_degree(node));
    }
    fprintf(stats, "e\n");
    fflush(stats);
}
//==================TOP K =================

double threshold_topk;
double zero_ppr_upper_bound;
iMap<double> ppr;
double rmax_scale = 1;
void pafora_topk_version(int s, double rsum, double rmax, int round) {
    
    // prepare for forward push
    source = s;
    int is_beginning = 1;
    Bag_reducer<int> *queue[2];
    Bag_reducer<int> b1;
    Bag_reducer<int> b2;
    queue[0] = &b1;
    queue[1] = &b2;
    bool queuei = 1;

    // put the source node into one of the bags
    (*queue[1]).insert(s);
    int iter = 0;
    int scan_iter = 0;
    int isdone = 0;
    int scan_size = 0;
    int forward_type = 0;
    int scan_iteration = 0;
    int scan_done = 0;
    int bag2_begin = 0;
    int scan_begin = 0;

    if (round > 0) forward_type = 1;
    // run forward push
    while( !isdone ) {
        
        iter++;
        if (forward_type > 0 || (*queue[queuei]).numElements() > graph.n/hybridsize) {
            
            if(!(*queue[queuei]).isEmpty()) {
                (*queue[queuei]).clear();
            }
            forward_type = 1;
            //scan_size = sep_scan_forward_with_sched(queue[queuei]);
            scan_size = sep_scan_forward(queue[queuei]);
            scan_iteration += 1;
            if (scan_size < graph.n/hybridsize) {
                forward_type = 0;
                scan_size = 0;
                scan_done = 1;
            }
            //continue;
        } else {
            isdone = 1; 
        }
        if (forward_type == 1)
            continue;

        if(!(*queue[!queuei]).isEmpty()) {
            (*queue[!queuei]).clear();
        }
        process_bag(&((*queue[queuei]).get_reference()), queue[!queuei]);
        
        double source_residue = s_residue.get_value();
        s_residue.set_value(0.0);
        long long degree_s = graph.get_degree(source);
        if (residue_s[source] >= graph.rmax_scaled * degree_s) {
            atomic_add(residue_s[source], source_residue);
        } else {
            atomic_add(residue_s[source], source_residue);
            long long degree_s = graph.get_degree(source);
            if (residue_s[source] > graph.rmax_scaled * degree_s) {
                Bag_reducer<int> *next = queue[!queuei];
                Bag<int>* bnext = &((*next).get_reference());

                (*bnext).insert(source);
            }
        }
        // prepare for next iteration
        queuei = !queuei;
        isdone = (*queue[queuei]).isEmpty();
    }
    // random walks (integer arrays initialized in init_s_reserve())
    // init ppr values (used for collecting random walk destinations)
    if (ppr_values == NULL) ppr_values = new atomic<double>[graph.n];
    memset(ppr_values, 0, sizeof(atomic<double>) * graph.n);
    // DEBUG
    //cilk::reducer< cilk::op_add<long int> > pcount_sum(0);
    //cilk_for (int i = 0; i < graph.n; i ++) {
    //    int node = i;
    //    for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
    //        *pcount_sum += reserve_numa_i32[j][node] +
    //                        reserve_numa_u8[j][node];
    //    }
    //}
    //printf("pcount sum = %ld, Omega = %llu, %f\n\n", 
    //    pcount_sum.get_value(), Omega,(double)pcount_sum.get_value()/Omega);
    // END OF DEBUG
    double omega = Omega * rsum;
    cilk::reducer< cilk::op_add<double> > residue_source(0.0);
    cilk::reducer< cilk::op_add<long> > rw_count(0);
    #pragma cilk grainsize=1
    cilk_for(int i = 0; i < graph.dividingline ; i += GRAINSIZE_MC) {
    //for(int i = 0; i < graph.dividingline; i += GRAINSIZE_MC) {
        int core_id = sched_getcpu();
        int numa_node = core_id % NCOPY;
        int start = i;
        int end = i + GRAINSIZE_MC;
        if (end > graph.dividingline) end = graph.dividingline;
        for (int j = start; j < end; j ++) {
            //mc_on_a_node(j, omega, numa_node, core_id);
            int node = j;
            double residue_i = residue_s[node];
            if (residue_i > 0) {
                unsigned long num_rw = floor(residue_i * Omega);
                long long int start_idx = graph.dest_pointers[node];
                for (unsigned long long j = 0; j < num_rw; j ++) {
                    int t = graph.destinations[start_idx + j];
                    int old = Atomic_add(reserve_numa_u8[numa_node][t], 1);
                    if (old == 255) {
                        atomic_add(ppr_values[i], 256.0/Omega);
                    }
                    *rw_count += 1;
                }
                if (ceil(residue_i * Omega) > floor(residue_i * Omega)) {
                     double alpha_i = (Omega*residue_i - num_rw) / Omega;
                     int t = graph.destinations[start_idx+num_rw];
                     atomic_add(ppr_values[t], alpha_i);
                }
            }
        }   
    }

    // handle nodes with no outgoing edges
    if (graph.n != graph.dividingline) {
        cilk_for (int i = graph.dividingline; i < graph.n; i ++) {
        //cilk_for (int i = nodes_dividingline; i < size_nodes; i ++) {
            
            //int node = nodes[i];
            int node = i;
            double residue_i = residue_s[node];
            //atomic_add(reserve_incre[source], increase);
            *residue_source += residue_i;
        }
        atomic_add(reserve_s[source], residue_source.get_value());
    }
    
    // merge results and compute ppr values

    #pragma simd
    cilk_for (int i = 0; i < graph.n; i ++) {
        int node = graph.target_map[i];
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            ((double*)ppr_values)[i] += 
                (double)(((uint8_t*)reserve_numa_u8[j])[node])/graph.Omega;
                //(((int*)reserve_numa_i32[j])[node]) = 0;
                //(((uint8_t*)reserve_numa_u8[j])[node]) = 0;
        }
        ((double*)ppr_values)[i] += ((double*)reserve_s)[i];
    }
    
    return; 
    
    cilk::reducer< cilk::op_add<double> > ppr_sum(0.0);
    cilk_for(int i = 0; i < graph.n; i ++) {
        *ppr_sum += ppr_values[i];
    }
    printf("ppr sum = %.10f\n", ppr_sum.get_value());

    cilk::reducer< cilk::op_add<double> > rr_sum(0.0);
    cilk_for(int i = 0; i < graph.n; i ++) {
        *rr_sum += reserve_s[i]+residue_s[i];
    }
    printf("rr sum = %.10f\n", rr_sum.get_value());

    cilk::reducer< cilk::op_add<double> > r_sum(0.0);
    cilk_for(int i = 0; i < graph.n; i ++) {
        *r_sum += residue_s[i];
    }
    printf("rsum = %.10f\n", r_sum.get_value());

    cilk::reducer< cilk::op_add<double> > p_sum(0.0);
    cilk_for(int i = 0; i < graph.n; i ++) {
        *p_sum += reserve_s[i];
    }
    printf("psum = %.10f\n", p_sum.get_value());

    cilk::reducer< cilk::op_add<long int> > count_sum(0);
    cilk_for (int i = 0; i < graph.n; i ++) {
        int node = i;
        for (int j = 0; j < NUMA_NODES_SIZE; j ++) {
            *count_sum += reserve_numa_u8[j][node];
        }
    }
    printf("num_rw sum is %ld\n", rw_count.get_value());
    printf("count sum = %ld, Omega = %llu, %f\n\n", 
        count_sum.get_value(), Omega, (double)count_sum.get_value()/Omega);
}

inline bool cmp(double x, double y) {
    return x > y;
}


// obtain the top-k ppr values from ppr_values
double kth_ppr(int k) {
    
    struct timespec kstart, kfinish;
    clock_gettime(CLOCK_REALTIME, &kstart);

    vector<double> temp_ppr;
    temp_ppr.clear();
    bool use_reducer = true;
    if (!use_reducer) {
        int length = graph.n;
        temp_ppr.resize(length);
        for (int i = 0; i < length; i ++) {    
            temp_ppr[i] = ppr_values[i];
        }
    } else {
        cilk::reducer< cilk::op_list_append<double> > nodes_reducer;
        cilk_for (int i = 0; i < graph.n; i ++) {
            if (ppr_values[i] > 0) {
                nodes_reducer->push_back(ppr_values[i]);
            }
        }
        const list<double> &nodes = nodes_reducer.get_value();
        for (list<double>::const_iterator i = nodes.begin(); i != nodes.end(); i ++) {
            temp_ppr.push_back(*i);
        }
    }

    nth_element(temp_ppr.begin(), temp_ppr.begin()+k-1, temp_ppr.end(), greater<double>());
    double kth = temp_ppr[k-1];
    //for (int i = 0; i < length; i ++) {
    //    if (temp_ppr[i] >= 0.2)
    //        cout << i <<" " <<temp_ppr[i] << endl;
    //}
    double elapsed = 0;
    clock_gettime(CLOCK_REALTIME, &kfinish);
    elapsed = (kfinish.tv_sec - kstart.tv_sec)*1000000000;
    elapsed += (kfinish.tv_nsec - kstart.tv_nsec);
    printf("find kth among %ld processed in %.3f miliseconds using %d workers.\n",temp_ppr.size(),
            (elapsed)/1000000, __cilkrts_get_nworkers());

    return kth;
}

bool if_stop(double delta, int k) {
    double kth = kth_ppr(k);
    if (kth >= 2.0 * delta) {
        printf("k = %d, kth = %.5e\n",k, kth);
        return true;
    }

    if (delta > threshold_topk) return false;

}
inline void set_fora_parameters(int n, long long m, double epsilon, 
                                double pfail, double delta){
    rmax = epsilon*sqrt(delta/3/m/log(2/pfail));
    rmax *= rmax_scale;
    graph.rmax_scaled = rmax;
    // rmax *= config.multithread_param;
    graph.Omega = Omega = (2+epsilon)*log(2/pfail)/delta/epsilon/epsilon;
}
void init_rw_counts() {
    
    for (int i = 0 ; i < NUMA_NODES_SIZE; i ++) {
        memset(reserve_numa_u8[i], 0, sizeof(uint8_t) * graph.n);
    }
}
void init_arrays_for_pafo(int source) {

    // remove all the elements in source node's reserve vector
    if (reserve_s == NULL) {
        reserve_s = new atomic<double>[graph.n];
        residue_s = new atomic<double>[graph.n];
        //NUMA aware reserve array
        NUMA_NODES_SIZE = numa_max_node()+1;
        int nworker = 80;
        printf("NUMA_NODES_SIZE = %d\n", NUMA_NODES_SIZE);
        reserve_numa_u8  = new atomic<uint8_t>*[NUMA_NODES_SIZE<<CACHEPAD_U8];
        for (int i = 0; i < NUMA_NODES_SIZE; i ++) {
            //void* memory_double = numa_alloc_onnode(sizeof(double) * graph.n, i);
            //reserve_numa_double[i] = (double*) memory_double;
            void* memory_u8 = numa_alloc_onnode(sizeof(atomic<uint8_t>) * (graph.n), i);
            reserve_numa_u8[i]  = (atomic<uint8_t>*) memory_u8;
        }
    }

    // set to 0
    //#pragma simd
    //cilk_for(int i = 0; i < graph.n; i ++) {
    //    ((double*)reserve_s)[i] = 0;
    //    ((double*)residue_s)[i] = 0;
    //}
    memset(reserve_s, 0, sizeof(atomic<double>) * graph.n);
    memset(residue_s, 0, sizeof(atomic<double>) * graph.n);
    // NUMA Aware
    for (int i = 0 ; i < NUMA_NODES_SIZE; i ++) {
        //#pragma simd
        //cilk_for(int j = 0; j < graph.n; j ++){
        //    ((uint16_t*)reserve_numa_u16[i])[j]= 0;
        //    ((int*)reserve_numa_i32[i])[j]= 0;
        //    ((uint8_t*)reserve_numa_u8[i])[j] = 0;
        //}
        memset(reserve_numa_u8[i], 0, sizeof(uint8_t) * (graph.n<<CACHEPAD_U8));
    }
    residue_s[source].store( 1); 
    s_residue.set_value(0.0);

}
void init_array_for_residues(int source) {
    int s = source;
    // may has been used for another source
    if (residue_s == NULL) {
        residue_s = new atomic<double>[graph.n];
    }

    // r(s, v) <-- 0 for all v != s
#ifdef PARALLEL_INIT
    cilk_for(int i = 0; i < graph.n; i ++) {
        residue_s[i] = 0;
    }
#else
    memset(residue_s, 0, sizeof(atomic<double>) * graph.n);
    //memset(residue_double, 0, sizeof(double)*graph.n);
#endif
    s_residue.set_value(0.0);
    // r(s, s) <-- 1
        #ifdef TEST_PORTION_FP
        if (init_residue_value == 0) {
            printf("please specify init residue value\n");
        }
    residue_s[s].store(source_portion);
        #else
    residue_s[s].store(1.0);
        #endif

    //residue_double[s] = 1;
    //init_local_residues();
    #ifdef FORWARD_PULL
    init_for_pull();
    #endif
    
}
void topk_query_with_bound(int source, int k, Result* result) {
    printf("-----source = %d\n", source);
    start_timer();
    //clock_gettime(CLOCK_REALTIME, &start);
    const static double min_delta = 1.0 / graph.n;
    const static double init_delta = 1.0 / 4;
    double ppr_decay_alpha = 0.77;
    double epsilon = 0.5;
    threshold_topk = (1.0-ppr_decay_alpha)/pow(500, ppr_decay_alpha) / pow(graph.n, 1-ppr_decay_alpha);

    //const static double new_pfail = 1.0 / graph.n / graph.n/log(graph.n/k);
    const static double new_pfail = 1.0 / graph.n / graph.n/log(graph.n);

    const static double lowest_delta_max = epsilon*sqrt(min_delta/3/graph.m/log(2/new_pfail));

    double rsum = 1.0;


            struct timespec kstart, kfinish;
            clock_gettime(CLOCK_REALTIME, &kstart);
    init_arrays_for_pafo(source);
    //init_array_for_residues(source);
    //init_s_reserve(source);
    zero_ppr_upper_bound = 1.0;

            double elapsed = 0;
            clock_gettime(CLOCK_REALTIME, &kfinish);
            elapsed = (kfinish.tv_sec - kstart.tv_sec)*1000000000;
            elapsed += (kfinish.tv_nsec - kstart.tv_nsec);
            printf("initialization processed in %.3f miliseconds using %d workers.\n",
                    (elapsed)/1000000, __cilkrts_get_nworkers());

            clock_gettime(CLOCK_REALTIME, &kstart);
    // for delta: try value from 1/4 to 1/n
    int round = 0;
    // iMap
    //upper_bounds.reset_one_values();
    //lower_bounds.reset_zero_values();
    
    //double delta = init_delta;
    double delta = threshold_topk;
    while ( delta >= min_delta ) {
        set_fora_parameters(graph.n, graph.m, epsilon, new_pfail, delta);
        // call PAFO : update reserve and residue array
        // Note: estimated ppr values should be calculated
        //       in each round
            struct timespec kstart, kfinish;
            clock_gettime(CLOCK_REALTIME, &kstart);

        pafora_topk_version(source, rsum, rmax, round++);

            clock_gettime(CLOCK_REALTIME, &kfinish);
            elapsed = (kfinish.tv_sec - kstart.tv_sec)*1000000000;
            elapsed += (kfinish.tv_nsec - kstart.tv_nsec);
            printf("    pafo processed in %.3f miliseconds using %d workers.\n",
                    (elapsed)/1000000, __cilkrts_get_nworkers());

            clock_gettime(CLOCK_REALTIME, &kstart);

        init_rw_counts();

            clock_gettime(CLOCK_REALTIME, &kfinish);
            elapsed = (kfinish.tv_sec - kstart.tv_sec)*1000000000;
            elapsed += (kfinish.tv_nsec - kstart.tv_nsec);
            printf("        init_rw_counts() processed in %.3f miliseconds using %d workers.\n",
                    (elapsed)/1000000, __cilkrts_get_nworkers());

        if (if_stop(delta, k) || delta <= min_delta) {
            break;
        } else {
            delta = max (min_delta, delta / 2.0);
        }
    }
    //printf("source= %.10f %d rounds\n", ppr_values[source].load(), round);

    print_timer("topk");
    printf("delta = %.5e\npfail= %.5e\n", delta, new_pfail);
    printf("Omega= %.5e\n", (double)Omega);
    printf("rmax = %.5e, round = %d\n",graph.rmax_scaled, round);

    //double elapsed = 0;
    clock_gettime(CLOCK_REALTIME, &finish);
    elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
    elapsed += (finish.tv_nsec - start.tv_nsec);
    result->add_topk(elapsed/1000000);

    //exit(0);
}


int main(int argc, char *argv[]) {

    cout << Green << "--------------start------------" << get_current_time_str() << Reset << endl;
    #ifdef ITERATION_COUNT
    iter_wl_counts = new double[max_iter];
    for (int i = 0; i < max_iter; i ++) {
        iter_wl_counts[i] = 0;
    }
    #endif
    //testParallelRand(40);
    //return 0;
    //testAtomic();
    //testAtomicInt();
    //testSomething();
    //return 0;
    nworker = -1;
    rounds = 20;
    string defaultdata = "../../../uqswan22/datasets/twitter";
    string datagraph = defaultdata;;
    const char* num_of_workers;
    int pouch_size_cmd = 0;
    double rwscale = 0;
    for(int i=0; i<argc; i++){
        if(string(argv[i])=="--worker"){
            num_of_workers = argv[i+1];
            __cilkrts_set_param("nworkers", num_of_workers);
            nworker = atoi(num_of_workers);
        }
        if(string(argv[i]) == "--dataset"){
            datagraph = string(argv[i+1]);   
            dataset = datagraph;
        }
        if(string(argv[i]) == "--nsetcover"){
            rounds = atoi(argv[i+1]);   
        }
        if(string(argv[i]) == "--rwmode"){
            rwmode = atoi(argv[i+1]);   
        }
        if(string(argv[i]) == "--rwpermute"){
            with_rw_permutation = atoi(argv[i+1]);   
        }
        if(string(argv[i]) == "--hybridsize"){
            hybridsize = atoi(argv[i+1]);   
        }
        if(string(argv[i]) == "--pouchsize"){
            pouch_size_cmd = atoi(argv[i+1]);   
        }
        if(string(argv[i]) == "--scale"){
            rwscale = atof(argv[i+1]);   
            rmax_scale = rwscale;
        }
        if(string(argv[i]) == "--topk"){
            run_topk = true;
        }
        if(string(argv[i]) == "--portion"){
            int p = atoi(argv[i+1]);    // this is used to test different 
            source_portion = 1/p;       //init residue  
        }
    }
    
	vector<int> query;
	load_ss_query(query, datagraph); 
    printf("query file has %ld sources\n", query.size());


    //numa_set_interleave_mask(numa_all_nodes_ptr);
    // load the graph
    printf("%s \n",datagraph.c_str());
    graph = Graph(datagraph);
    alpha = graph.alpha;
    //numa_set_localalloc();
    #ifdef FPWORKLOADCOUNT
    int query_num = 100;
    #else
    int query_num = 20;
    #endif
    // check popular nodes
    #ifdef CHECK_POP_NODES
    query_num = 1000;
    string query_file_str = datagraph+"/ssquery.txt";
    #else
    query_num = 20;
    #ifdef ITERATION_COUNT
    query_num = 1;
    #endif
    string query_file_str = datagraph+"/heavyquery.txt";
    //query_num = 6;
    //string query_file_str = datagraph+"/slowquery.txt";
    #endif

    #ifdef FP_RESIDUE_DISTRIBUTION
    query_num = 1000;
    query_file_str = datagraph+"/ssquery_2000.txt";
    #endif

    #ifdef SOURCE_PLAW
    query_num = 150;
    query_file_str = datagraph+"/ssquery_plaw.txt";
    #endif
    
    vector<int> sources;
    if(file_exists_test(query_file_str)){
         ifstream query_file(query_file_str.c_str());
         int v;
         while(query_file>>v) {
             //int node = graph.sep_mapping[v];
             int node = v;
             if (graph.get_degree(node) > 0) {
                sources.push_back(v);
             }
         }
         
    }else{
        printf("=========query file %s not found!\n", query_file_str.c_str());
        //exit(0);
        srand (time(NULL));
        ofstream query_file(query_file_str.c_str());
        #ifdef SOURCE_PLAW
        for (int i = 0; i < graph.n; i ++) {
            double rn = drand();
            if (rn <= 1.0*graph.get_degree(i)*query_num/graph.m) {
                sources.push_back(i);
            }
        }
        #else
        for(int i=0; i<query_num; i++){
            int v = rand()%graph.n;
            int mapped_v = graph.sep_mapping[v]; 
            while(graph.get_degree(mapped_v)==0){
                v = rand()%graph.n;
                mapped_v = graph.sep_mapping[v];
            }
            //printf("mapped_v = %d\n", mapped_v);
            sources.push_back(v);
        }
        #endif
            
        for (int i = 0; i < sources.size(); i ++)
            query_file<<sources[i]<<endl;
        printf("=========query file %s created!\n", query_file_str.c_str());
    }

    printf("valid source number is % ld\n", sources.size());
    #ifdef SOURCE_PLAW
    query_num = 50;
    query_file_str = datagraph+"/ssquery_plaw.txt";
    #endif

    // INIT residue and reserve arrays
    init_s_reserve(-1);
    init_residue_s(-1);

    graph.scale = rwscale;
    #ifdef PAGERANK_ORDER
    pagerank = PageRank();
    pagerank.calculate_exact_pagerank(graph);

    graph.order_by_pagerank(pagerank.pg_values);
    #endif
    load_exact_topk_ppr();
    double n = graph.n;
    double epsilon = 0.5;
    double delta = 1 / n;
    double pfail = 1 / n;
    long long m = graph.m;
    // test topk parameters
    double ppr_decay_alpha = 0.77;
    const static double min_delta = 1.0 / graph.n;
    threshold_topk = (1.0-ppr_decay_alpha)/pow(500, ppr_decay_alpha) / pow(graph.n, 1-ppr_decay_alpha);

    const static double new_pfail = 1.0 / graph.n / graph.n/log(graph.n/500);

    const static double lowest_delta_max = epsilon*sqrt(min_delta/3/graph.m/log(2/new_pfail));
    if (run_topk) {
      pfail = new_pfail;
      delta = threshold_topk;
    }
    // end of test
    rmax = epsilon * sqrt(delta / 3 / (double)m/log(2/pfail));
    
    Omega = (2 + epsilon)*log(2/pfail)/delta/epsilon/epsilon;
    printf("rmax = %.5e\n", rmax*rmax_scale);
    printf("Omega = %llu\n", Omega);
    printf("Omega * rmax = %f\n", Omega*rmax);

    #ifdef NUMA_TARGET 
    // load popularity
    graph.load_popularity();
    #endif
    //graph.compute_threshold(rmax, 0);
    // must separate verices before finding set covers
    #ifdef SEPARATE
    bool sep = 1;
    #else
    bool sep = 0;
    #endif
    dividingline = graph.separate_vertices(sep);
    // For forward push schedule
    printf("Cache schedule called\n");
#ifdef CACHE_SCHEDULE

    #ifdef BLK_LOCALITY_FIRST
    graph.generate_block_order_locality_first();
    #else
    //graph.generate_block_order_by_cache_list();
    graph.permute_blocks_by_contention_locality_pop();
    #endif // BLK_LOCALITY_FIRST
#endif

    printf("Cache schedule ends\n");
    // SET COVER
    int rounds = 250;
    #ifdef SET_COVER
    printf("set cover rounds = %d\n", rounds);
    #endif
    graph.rounds = rounds;
    #ifdef SET_COVER
    graph.gen_multi_cover(rounds);
    graph.check_set_covers();
    #endif
    #ifdef SET_COVER_DUP
    graph.gen_multi_cover_dup(rounds);
    graph.check_set_covers_dup();
    #endif
    #ifdef SET_COVER_PR
    graph.gen_pagerank_multi_cover(rounds);
    #endif
    // END OF SET COVER

    #ifdef CACHE_SCHEDULE
    graph.numa_aware_load();    // this function has to be called after hotsetschedule called
    residue_s = graph.residues;
    reserve_s = graph.reserves;
    #endif
    // Forward push threshold
    if (rwscale == 0) {
        //graph.compute_threshold(rmax, sep, Omega, RMAX_TYPE, (double)RWINDEXSCALE*Omega*rmax);
        graph.compute_threshold(rmax, sep, Omega, RMAX_TYPE, (double)RWINDEXSCALE);
    } else {
        //graph.compute_threshold(rmax, sep, Omega, RMAX_TYPE,(double)rwscale*Omega*rmax);
        graph.compute_threshold(rmax, sep, Omega, RMAX_TYPE,(double)rwscale);
    }
    // permute for random walk scheduling
    #ifdef RWPERMUTE
    graph.permute_nodes_for_rw();
    #endif
    // generate random walk destination 
    #ifdef RANDOMWALK
    //GRAINSIZE_MC = graph.dividingline / (128) / __cilkrts_get_nworkers();
    int nscale = graph.n / 1000000;
    nscale /= 8;
    if (nscale == 0) nscale = 1;
    GRAINSIZE_MC = nscale*19200*3/(2*graph.m/graph.n);
    //GRAINSIZE_MC = 19200*3/(2*graph.m/graph.n);
    printf("rw grain size is %d\n", GRAINSIZE_MC);
    #ifndef FP_RESIDUE_DISTRIBUTION
        #ifdef RWINDEX
    graph.generate_random_walks();  // should be called after separation
        #endif
    // the following will change destination index, can be used for Gorder
    #ifdef GORDER
    //graph.generate_permuted_random_walk_index();
    #endif
    printf("random walk destinations ready\n");
    #endif
    #endif
    // generate rw set covers
    //graph.generate_multiple_rw_set_covers(10);
    //set pouch size
    #ifdef USEPOUCH
    nscale = graph.n / (1<<18);
    POUCH_SIZE = (graph.m/graph.n)*16; // 16 is INTIIAL_P_SIZE in pouch.h
    //POUCH_SIZE = 768;
    if (pouch_size_cmd != 0) {
        POUCH_SIZE = pouch_size_cmd;
    }
    printf("pouch size is %d \n", POUCH_SIZE);
    #endif
    cout<<"ready for query"<<endl;
    usleep(1000000);
    printf("query start\n");
    // test SSPPR
    //test1(s);
    //test1(s);
    //return 0;
    printf("query number is %d\n", query_num);
    vector<int> num_core;
    if (nworker != -1) {
        
        Result *result = new Result(nworker, query_num);
        //for(int j=0; j<query_num; j++){
        //    //scaletest(nworker, sources[j], result);
        //    #ifndef FP_RESIDUE_DISTRIBUTION
        //    if (sep) {
        //        int node = sources[j];
        //        #ifdef GORDER
        //        node = graph.gorder[node];
        //        #endif
        //        scaletest(nworker, graph.sep_mapping[node], result);
        //    } else {
        //        scaletest(nworker, sources[j], result);
        //    }
        //    #endif
        //}
            for(int j=0; j<query_num; j++){
                source_origin = sources[j];
                if (sep) {
                    int node = sources[j];
                    #ifdef GORDER
                    node = graph.gorder[node];
                    #endif
                    #ifdef MULTICORE
                    //printf("source = %d\n", graph.sep_mapping[node]);
                    #endif

                    #ifndef FP_RESIDUE_DISTRIBUTION
                    node = graph.sep_mapping[node];
                    scaletest(nworker, node, result);
                    #endif
                } else {
                    #ifdef MULTICORE
                    printf("source = %d\n", graph.sep_mapping[sources[j]]);
                    #endif
                    #ifndef FP_RESIDUE_DISTRIBUTION
                    scaletest(nworker, sources[j], result);
                    #endif
                }
            }
        result->print(1,1);
        #ifdef ITERATION_INFO
        printf("ITERATION average bag 1 phase : %.1f\n", total_iter1);
        printf("ITERATION average scan  phase : %.1f\n", total_iter2);
        printf("ITERATION average bag 2 phase : %.1f\n", total_iter3);
        #endif
        #ifdef BAG_PHASE_TIMER
        printf("average bag 2 phase : %f\n", total_bag2_time/query_num);
        #endif
    } else {
    #ifndef FP_RESIDUE_DISTRIBUTION
        num_core.push_back(1);
        num_core.push_back(4);
        num_core.push_back(8);
        num_core.push_back(12);
        num_core.push_back(16);
        num_core.push_back(20);
        num_core.push_back(24);
        num_core.push_back(28);
        num_core.push_back(32);
        num_core.push_back(36);
        num_core.push_back(40);
//        num_core.push_back(80);
//        num_core.push_back(80);
//        num_core.push_back(1);
//        num_core.push_back(1);
    #else // FP_RESIDUE_DISTRIBUTION
        num_core.push_back(1);
    #endif //FP_RESIDUE_DISTRIBUTION
        Result** results = new Result*[num_core.size()];

        int cpu_core_number = get_nprocs_conf()/2;
        #ifdef FP_RESIDUE_DISTRIBUTION
        MultiSources* msource = new MultiSources(cpu_core_number);
        double* total_count = new double[graph.n];
        #endif
        printf("%d cpu_cores\n", cpu_core_number);

        for(int i=0; i<num_core.size(); i++){
            //int core= num_core[num_core.size()-1-i];
            int core= num_core[i];
            results[i] = new Result(core, query_num);
            results[i]->set_sources(sources);
            #ifdef FP_RESIDUE_DISTRIBUTION
            // timer
            struct timespec start, finish;
            double elapsed;
            clock_gettime(CLOCK_REALTIME, &start);

            setWorkers(cpu_core_number);
            //setWorkers(1);
            cilk::reducer_opadd<int> ss_workload(0);
            cilk::reducer_opadd<int> fp_count(0);
            int k;
            printf("query_num = %d, cpu_core_number = %d\n", query_num, cpu_core_number);
            for (k = 0; k < query_num/cpu_core_number; k ++) {
                #pragma cilk grainsize=1
                cilk_for(int l = 0; l < cpu_core_number; l ++) {
                //for(int l = 0; l < cpu_core_number; l ++) {
                    int node = sources[k*cpu_core_number+l];
                    scaletest(core, graph.sep_mapping[node], 
                              results[i],
                              msource->residues[l],
                              msource->reserves[l],
                              msource->push_counts[l]);
                    msource->merge_result_and_reset(l);
                    ss_workload++;
                }
                printf("k = %d\n", k);
            }
            #pragma cilk grainsize=1
            cilk_for(int l = 0;l<query_num%cpu_core_number;l++) {
            //for(int l = 0;l<query_num%cpu_core_number;l++) {
                int node = sources[k*cpu_core_number+l];
                scaletest(core, graph.sep_mapping[node], 
                          results[i],
                          msource->residues[l],
                          msource->reserves[l],
                          msource->push_counts[l]);
                msource->merge_result_and_reset(l);
                ss_workload++;
            }
            
            printf("%d sources processed\n", ss_workload.get_value());
            printf("total counts is %lld\n", msource->total_num.load());
            printf("total counts is %d\n", fp_count.get_value());
            // timer
            clock_gettime(CLOCK_REALTIME, &finish);
            elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
            elapsed += (finish.tv_nsec - start.tv_nsec);
            printf("Forward push calculated in %.3f miliseconds using %d workers.\n",
                   elapsed/1000000, __cilkrts_get_nworkers());
            printf("\n");
            // sort stats result for forward push and write to file
            write_forwardpush_stats(msource, datagraph);
            write_forwardpush_push_counts(msource, datagraph);
            #else
            for(int j=0; j<query_num; j++){
                //scaletest(core, sources[j], results[i]);
                if (sep) {
                    //int node = graph.gorder[sources[j]]; // Gorder
                    int node = sources[j];
                    #ifdef GORDER
                    node = graph.gorder[node];
                    #endif
                    //#ifdef MULTICORE
                    //printf("source = %d\n", graph.sep_mapping[node]);
                    //#endif
                    //printf("source = %d core = %d, i = %d, j = %d, query.size = %d %d\n", node, core, i,j, query.size(), query_num);
                    scaletest(core, graph.sep_mapping[node], results[i]);
                } else {
                    #ifdef MULTICORE
                    printf("source = %d\n", sources[j]);
                    #endif
                    //scaletest(core, graph.gorder[sources[j]], results[i]); // Gorder
                    #ifndef FP_RESIDUE_DISTRIBUTION
                    scaletest(core, sources[j], results[i]);
                    #endif
                }
            }
            results[i]->print(1, 1);
            #endif

            #ifdef CHECK_POP_NODES
            break;
            #endif
        }
        #ifndef FP_RESIDUE_DISTRIBUTION
        for (int i = 0; i < query_num ; i++) {
            //double f1 = results[0]->print_source(i);
            //double f2 = results[1]->print_source(i);
            //printf("=======fp speed up = %f\n", f1/f2);
        }
        for (int i = 0; i < num_core.size(); i ++) {
            //results[i]->print(1, 1);
            #ifdef PHASE_TIMER
            printf("average cilk for iteration is %d\n", total_iter/query_num);
            printf("average phase1 time: %f\n", total_time/query_num);
            printf("average phase2-3 time: %f\n", total_time2/query_num);
            #endif
            
            #ifdef BAG_PHASE_TIMER
            printf("average bag 2 phase : %f\n", total_bag2_time/query_num);
            #endif
            
            #ifdef ITERATION_INFO
            printf("ITERATION averate bag 1 phase : %.1f\n", total_iter1/query_num);
            printf("ITERATION averate scan  phase : %.1f\n", total_iter2/query_num);
            printf("ITERATION averate bag 2 phase : %.1f\n", total_iter3/query_num);
            #endif
        }
        #endif
        // delete resutls
        for (int i = 0; i < num_core.size(); i ++) {
            delete results[i];
        }
        delete[] results;
    }
    delete_arrays();
    graph.free_memory();

    #ifdef ITERATION_COUNT
    for (int i = 0; i < 50; i ++) {
        printf("%0.15f nodes processed in iter %d \n", iter_wl_counts[i]/(graph.n*query_num), i);
    }
    for (int i = 0; i < 50; i ++) {
        printf("%.0f nodes processed in iter %d \n", iter_wl_counts[i], i);
    }
    #endif
    
    if (rwmode == 0) {
        printf("random walk by atomic<double>\n");
    } else {
        printf("random walk by atomic<int>\n");
    }
    printf("hybrid size is %f\n", hybridsize);
    printf("==================================\n");
    
    //for (int i = 0; i < query_num; i ++) {
    //    printf("%d ", graph.sep_mapping[sources[i]]);
    //}
    printf("\n");
    cout << Red << "--------------stop------------" << get_current_time_str() << Reset << endl << endl << endl;

}
