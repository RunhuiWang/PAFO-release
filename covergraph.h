// ATTENTION! 64 is for twitter only
//#define TWITTER 1
//#define FRIENDSTER 1

//#define FRIENDSTER_SCALE 25

#ifdef TWITTER
#define GRAINSIZE_SCALE 8*8
#elif defined FRIENDSTER
#define GRAINSIZE_SCALE 16*8
#else
//#define GRAINSIZE_SCALE 8
#define GRAINSIZE_SCALE 8*2
#endif
#ifndef __GRAPH_H__
#define __GRAPH_H__

#ifndef USE_DOTMIX
#define USE_DOTMIX 1
#endif

#include <atomic>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include <memory.h>
//#include <numaif.h>
//#include <linux/mempolicy.h>
//#include <numa.h>
//#define OLD_VERSION_COMPARISON
#include <sched.h>
//#include <linux/kernel.h>

//#define FAST_IO 1
#define FAST_INDEX 1
//MMAP
#include <stdint.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
//#include <thread>

#ifdef USE_DOTMIX
#include <cilkpub/dotmix.h>
#endif

#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
//#define CHECKING

// parallel graph IO
//#define PARALLEL_LOAD_GRAPH

#define PARALLEL_LOAD_IDX
#include "paLoad.h"
#include "cilkMin.h"

//#define BLK_LOCALITY_FIRST 1
#define CACHELINE_SIZE 64
#define SET_COVER_NROUND 100

#ifdef USE_DOTMIX
    #ifdef DS_TEST
extern cilkpub::DotMix global_rng;    
    #else
cilkpub::DotMix global_rng(234909128);
    #endif
#endif

/*the input edges (u,v) must be sorted by increasing order of u (tiers broken arbitarily if u is the same)!!!!
Otherwise, it will be incorrect*/
using namespace std;

extern double hybridsize;

class Graph {
public:
    int rounds = SET_COVER_NROUND;
    //assume that # of nodes fits into int, and # of edges fits into long long
    double rmax;
    double rmax_scaled;
    double scale=1;
    double Omega;
    int thres_type;
    vector<double> thres;
    //vector<double> r_thres;
    #ifdef PHASE_TIMER
    int phase1;
    #endif
    vector<vector<int> > g;
    vector<vector<int> > gr;
    vector<int> rw_idx;
    vector< pair<unsigned long long, unsigned long> >rw_idx_info;

#ifdef FAST_IO
#endif // FAST_IO
    atomic<long int>* in_degree;
    int* edges[2]={NULL, NULL};
    long* nodepointer[2]= {NULL, NULL};
    int* degree[2]{NULL, NULL};
    int d_in_max;
    vector<int> pagerank_order;
    vector<pair<int,double> > pagerank_pair;
    vector<int> pr_reverse_map;

    bool* istrouble = NULL;    // If a node's out degree is too large, then this node is a trouble
    vector<int> degree_order;   // degree map
    vector<pair<int, long int> > degree_pair;
    vector<int> degree_reverse_map;
    int* random_walk_map = NULL;

    string data_folder;
    double dbar;
    vector<double> global_ppr;
    // the tele ratio for random walk
    double alpha;
    vector<int> gorder;
    vector<int> reverse_map;
    // seperated forward push
    int* sep_mapping=NULL;   // seperate vertices without out outneighbors
    vector<int> sep_reverse_map;
    vector<vector<int> > g_sep;
    vector<vector<int> > gr_sep;
    long long int dividingline;
    // set cover
    int* cover_map=NULL;
    int* cover_pos=NULL;     // reserve of cover_map
    vector<int> set_cover;
    int DEGREE_THRES = 100;
    double threshold_scale;

    inline int get(int direction, int node, int pos){
        long nodepos = nodepointer[direction][node]+pos;
        return edges[direction][nodepos];
    }

    static bool cmp(const pair<int, double> &t1, const pair<int, double> &t2) {
        return t1.second > t2.second;
    }

    int n;
    long long m;
    long long scanned_edges;

    inline long int get_degree(int n) {
        return nodepointer[0][n+1] - nodepointer[0][n];
    }
    inline long int get_in_degree(int n) {
        return in_degree[n];
    }
    inline long int get_edges_start(int n) {
        return nodepointer[0][n];
    }
    Graph(string data_folder) {
        //INFO("sub constructor");
        this->data_folder = data_folder;
        //this->alpha = ALPHA_DEFAULT;
        this->alpha = 0.2;
        //init_graph_handling_dead_end();
        init_graph();
        thres.assign(n, 0.0);
        //r_thres.assign(n, 0.0);
        cout << "init graph n: " << this->n << " m: " << this->m << endl;
        cout<< "Avg degree: "<<(m/n)<<"  log(n): "<<(log(n))<<endl;
        sync();

        std::ofstream ofs("/proc/sys/vm/drop_caches");
        ofs << "1" << std::endl;
    }
    Graph(){}
    
    void free_memory(){
        #ifdef PARALLEL_LOAD_GRAPH
        pagraph.del();
        #endif
        if (nodepointer[0]  != NULL) delete[] nodepointer[0];
        if (nodepointer[1]  != NULL) delete[] nodepointer[1];
        if (degree[0]       != NULL) delete[] degree[0];
        if (degree[1]       != NULL) delete[] degree[1];
        if (istrouble       != NULL) delete[] istrouble;
        #ifndef PARALLEL_LOAD_GRAPH
        if (edges[0]        != NULL) delete[] edges[0];
        #endif
        if (edges[1]        != NULL) delete[] edges[1];
        if (sep_mapping     != NULL) delete[] sep_mapping;
        if (cover_map       != NULL) delete[] cover_map;
        if (cover_pos       != NULL) delete[] cover_pos;
        if (is_in_sets      != NULL) delete[] is_in_sets;
        if (is_node_covered != NULL) delete[] is_node_covered;
        if (is_edge_covered != NULL) delete[] is_edge_covered;
        if (dest_pointers   != NULL) delete[] dest_pointers;
        if (destinations    != NULL) delete[] destinations;
        if (is_occupied     != NULL) delete[] is_occupied;
        if (permute_mapping != NULL) delete[] permute_mapping;
        if (permute_mapping_r   != NULL) delete[] permute_mapping_r;
        if (gmap_destinations   != NULL) delete[] gmap_destinations;
        if (gmap_dest_pointers  != NULL) delete[] gmap_dest_pointers;
        if (random_map      != NULL) delete[] random_map;
        if (random_map_r    != NULL) delete[] random_map_r;
        if (target_map      != NULL) delete[] target_map;
        if (target_map_r    != NULL) delete[] target_map_r;
        if (grain_blocks    != NULL) {
            for (int i = 0; i < NROUNDS; i ++) 
                delete[] grain_blocks[i];
            delete[] grain_blocks;
        }
        if (indices         != NULL) delete[] indices;
        if (block_selected  != NULL) delete[] block_selected;
        if (nrounds         != NULL) delete[] nrounds;
    }

    void compute_threshold(double rmax, bool sep, double Omega, int type, double scale) {
        this->rmax = rmax;
        this->rmax_scaled = rmax*scale;
        this->Omega = Omega;
        this->thres_type = type;
        threshold_scale = scale;
        printf("threshold scale is %f\n", threshold_scale);
        printf("scaled rmax is %.15f\n", rmax_scaled);
        cilk_for (int i = 0; i < n; i ++) {
            if (sep)
                thres[i] = rmax * get_degree(i);
            else
                thres[i] = rmax * get_degree(i);
            //thres[i] = rmax * r_degrees[i];
            // Type 1 means to use D(v)/Omega as threshold
            if (type == 1) {
                long int size = get_degree(i);
                thres[i] = size/Omega;// * log(size);
            }
            thres[i] *= scale;
        }
    }

    void init_nm() {
        //string attribute_file = data_folder + FILESEP + "attribute.txt";
        string attribute_file = data_folder + "/attribute.txt";
        assert_file_exist("attribute file", attribute_file);
        ifstream attr(attribute_file);
        string line1, line2;
        char c;
        while (true) {
            attr >> c;
            if (c == '=') break;
        }
        attr >> n;
        while (true) {
            attr >> c;
            if (c == '=') break;
        }
        attr >> m;
        dbar = m*1.0/n;

    }

    #ifdef PARALLEL_LOAD_GRAPH
    paGraph<int> pagraph;
    #endif
    #ifdef PARALLEL_LOAD_IDX
    paGraph<int> paidx;
    #endif

    void map_file_to_memory(string filename) {
        struct stat sb;
        int fd = open(filename.c_str(), O_RDONLY);
        // get the size in bytes of the file
        fstat (fd, &sb);
 
        // map the file in a memory area
        void *pointer = mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
 
        char* p = (char*) pointer;
        // print 3 char of the file to demostrate it is loaded ;)
        printf("first 3 chars of the file: %c %c %c\n", p[0], p[1], p[2]);
 
        close(fd);
    }
    void deserialize_graph() {
        
        string object_file= data_folder + "/graph_obj.txt";
        string object_file_r= data_folder + "/graph_obj_r.txt";
        if (exists_test( object_file )) {
        //numa_set_interleave_mask(numa_all_nodes_ptr);
            std::ifstream ifs(object_file);
            boost::archive::binary_iarchive ia(ifs);
            ia >> g;
        //numa_set_localalloc();
            // length of nodepointer is n+1
            nodepointer[0][0] = 0;
            for (int i = 1; i <= n; i ++) {
                nodepointer[0][i] = nodepointer[0][i-1] +
                                    g[i-1].size();
                if (i < 10) {
                    cout << nodepointer[0][i] << endl;
                }
            }
            edges[0] = new int[m];
            cilk_for(int i = 0; i < n; i ++) {
                long long start = nodepointer[0][i];
                long length = get_degree(i);
                for (long j = 0; j < length; j ++) {
                    int node = g[i][j];
                    edges[0][start+j] = node;
                }
                degree[0][i] = g[i].size();
            }
            printf("m = %ld, nodepointer[0][n] = %ld\n",m,nodepointer[0][n]);
            // prepare nodepointer array
            cilk_for(int i = 0; i < n; i ++) {
                g[i].clear();
                vector<int>().swap(g[i]);
            }
            g.clear();
            ifs.clear();
            // init in_degree array
            in_degree = new atomic<long int>[n];
            memset(in_degree, 0, sizeof(atomic<long int>)*n);
            cilk_for(long int i = 0; i < m; i ++) {
            //    in_degree[edges[0][i]] += 1;
            }
            return;
            std::ifstream ifs_r(object_file_r);
            boost::archive::binary_iarchive ia_r(ifs_r);
            ia_r >> gr;
        }
    }

    inline string get_idx_file_name(){
        string file_name;
        if(scale==1)
            file_name = data_folder+"/randwalks.idx";
        else
            file_name = data_folder+"/randwalks."+to_string(scale)+".idx";
        
        return file_name;
    }
    
    inline string get_idx_info_name(){
        string file_name;
        if(scale==1)
            file_name = data_folder+"/randwalks.info";
        else
            file_name = data_folder+"/randwalks."+to_string(scale)+".info";
        return file_name;   
    }
    void deserialize_rw_idx() {
        
        string file_name = get_idx_info_name();
        assert_file_exist("index file", file_name);
        std::ifstream info_ifs(file_name);
        boost::archive::binary_iarchive info_ia(info_ifs);
        info_ia >> rw_idx_info;
        printf("index info size %ld\n", rw_idx_info.size());
        
        file_name = get_idx_file_name();
        assert_file_exist("index file", file_name);
        printf("loading random walk index...\n");
        std::ifstream ifs(file_name);
        boost::archive::binary_iarchive ia(ifs);
        ia >> rw_idx;
        printf("index size %ld\n", rw_idx.size());
        cilk_for(int i = 0; i < n; i ++) {
            dest_pointers[i] = rw_idx_info[i].first;
        }
        dest_pointers[n] = rw_idx_info[n-1].first + rw_idx_info[n-1].second;
        cout << dest_pointers[n] << endl;
        destinations = new int[dest_pointers[n]];
        cilk_for (long long i = 0; i < dest_pointers[n]; i ++) {
            destinations[i] = rw_idx[i];
        }
        for (int i = 0; i < n; i += 10000000) {
            for (long long j = dest_pointers[i]; j < dest_pointers[i]+10; j ++) {
                cout << destinations[j] << " ";
            }
            cout << endl;
        }
        rw_idx.clear();
        vector<int>().swap(rw_idx);
        
        rw_idx_info.clear();
        vector< pair<unsigned long long, unsigned long> >().swap(rw_idx_info);
        info_ifs.clear();
        //exit(0); 
    }
    void init_graph() {

        init_nm();
        nodepointer[0]=  new long[n+1];
        nodepointer[1]= new long[n+1];
        degree[0] = new int [n];
        degree[1] = new int[n];
        istrouble = new bool[n];
        memset(nodepointer[0], 0, sizeof(long) * (n+1));
        memset(nodepointer[1], 0, sizeof(long) * (n+1));
        memset(istrouble, false, sizeof(bool) * n);
        memset(degree[0], 0, sizeof(int)*n);
        memset(degree[1], 0, sizeof(int)*n);
        scanned_edges=0;

        g_sep = vector<vector<int>> (n, vector<int>());
        gr_sep = vector<vector<int>> (n, vector<int>());
        
      #ifdef FAST_IO
        deserialize_graph();
        return;
      #endif
    #ifdef PARALLEL_LOAD_GRAPH
        // load adjgraph in parallel
        std::string string_path = (data_folder+"/graph_adj.txt");
        //map_file_to_memory(string_path);
        //char* file_path = new char[string_path.size()];
        //std::copy(string_path.begin(), string_path.end(), file_path);
        pagraph = readGraphFromFile<int>(string_path.c_str());
        n = pagraph.n;
        m = pagraph.m;
        edges[0] = pagraph.edges;
        printf("adj graph read! \n");
        in_degree = new atomic<long int>[n];
        memset(in_degree, 0, sizeof(atomic<long int>)*n);
        cilk_for(long int i = 0; i < m; i ++) {
            in_degree[edges[0][i]] += 1;
        }
    #else
        cout << "-=-=-=-=-=-=read graph by scan" << endl;
        //numa_set_interleave_mask(numa_all_nodes_ptr);
        edges[0] = new int[m];
        //edges[1] = new int[m];
        //memset(edges[0], 0, sizeof(int)*m);
        //memset(edges[1], 0, sizeof(int)*m);
        //bool graph_split = false;
        #ifdef GORDER
        bool graph_split = false;
        #else
        bool graph_split = true;
        #endif
        if (graph_split) {
            cout << "-=-=-=-=-=-= graph split" << endl;
            FILE *split_info = fopen((data_folder+"/graph_split_info.txt").c_str(), "r");
            int npiece = 0;
            fscanf(split_info, "%lld", &npiece);
            long long* edge_starts = new long long [npiece];
            for (int i = 0; i < npiece; i ++) {
                fscanf(split_info, "%lld", &edge_starts[i]);
                cout << edge_starts[i] << " ";
            }
            cout << endl;
            fclose(split_info);
            setWorkers(npiece);
            #pragma cilk grainsize=1
            cilk_for(int i = 0; i < npiece; i ++) {
                cpu_set_t  mask;
                CPU_ZERO(&mask);
                CPU_SET(i, &mask);
                sched_setaffinity(0, sizeof(mask), &mask);
                cout << sched_getcpu() << endl;

                string file_name = data_folder+"/graph_part";
                stringstream ss;
                ss << i;
                string n_part = ss.str();
                file_name = file_name + n_part + ".txt";
                printf("%s\n", file_name.c_str());
                ss.clear();

                FILE *fin= fopen(file_name.c_str(), "r");
                int t1, t2;
                long long offset = edge_starts[i];
                while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
                    if (t1 == t2) {
                        m--; continue;
                    }
                    degree[0][t1] ++;
                    degree[1][t2] = 0;
                    edges[0][offset++] = t2;
                }
                fclose(fin);
                cout << offset << endl;
                CPU_ZERO(&mask);
                for (int i = 0; i < sysconf(_SC_NPROCESSORS_ONLN); i ++)
                    CPU_SET(i, &mask);
                sched_setaffinity(0, sizeof(mask), &mask);
                
            }
            setWorkers(sysconf(_SC_NPROCESSORS_ONLN));
//            // check correctness
//            string graph_file = data_folder + "/graph.txt";
//            assert_file_exist("graph file", graph_file);
//            FILE *fin = fopen(graph_file.c_str(), "r");
//            int t1, t2;
//            long long index = 0;
//            while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
//                if(t1 == t2) {
//                    m--;
//                    continue;
//                }
//                if (edges[0][index++] != t2) {
//                    cout << edges[0][index-1] << " " << t2 << endl;
//                    exit(0);
//                }
//            }
//            cout << "split graph is consistent with grahp.txt" << endl;
        } else {
		    #ifdef GORDER
            string graph_file = data_folder + "/graph_Gorder.txt";
		    #else
            string graph_file = data_folder + "/graph.txt";
		    #endif
            cout << graph_file << endl;
            assert_file_exist("graph file", graph_file);
            FILE *fin = fopen(graph_file.c_str(), "r");
            int t1, t2;
            int last_node;
            while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
                if(t1 == t2) {
                    m--;
                    continue;
                }
                degree[0][t1]++;
                degree[1][t2]++;
                edges[0][scanned_edges]= t2;
                scanned_edges++;
            }
        }
        in_degree = new atomic<long int>[n];
        memset(in_degree, 0, sizeof(atomic<long int>)*n);
        cilk_for(long int i = 0; i < m; i ++) {
            in_degree[edges[0][i]] += 1;
        }
        //numa_set_localalloc();
    #endif // PARALLEL_LOAD_GRAPH
        #ifdef PARALLEL_LOAD_GRAPH
        cilk_for(int i = 0; i < n; i ++) {
            degree[0][i] = pagraph.V[i].degree;
        }
        d_in_max = 0;
        //for (int i = 0; i < 10; i ++) {
        //    printf("degree[0][i] = %ld\n", degree[0][i]);
        //}
        #endif
        nodepointer[0][0] = 0;
        for(int i=1; i<n; i++){
            nodepointer[0][i]=nodepointer[0][i-1]+degree[0][i-1];
            //nodepointer[1][i]=nodepointer[1][i-1]+degree[1][i-1];
            if(degree[1][i]>d_in_max) d_in_max = degree[1][i];
        }
        
        // to avoid boundary check
        nodepointer[0][n] = m;
        printf("scanned edges =  %lld m = %lld\n", scanned_edges, m);
        
    #ifdef CHECKING
        // check adj graph
        for (int i = 0; i < n; i ++) {
            int adj_degree = pagraph.V[i].degree;
            if (adj_degree != degree[0][i]) {
                printf("i = %d, ajd_degree = %d, degree = %d\n",
                        i, adj_degree, degree[0][i]);
                exit(0);
            }
        }
        printf("degree checked\n");
        // check edges
        for (long long i = 0; i < m; i ++) {
            int node = pagraph.edges[i];
            if (edges[0][i] != node) {
                printf("i = %d, edges_adj = %d, edges= %d\n",
                        i, node, degree[0][i]);
                exit(0);
            }
        }
        printf("edges checked\n");
    #endif
        /*
        vector<int> temp_degree_count = vector<int>(n,0);
        for(int i=0; i<n; i++){
            for(int j=0; j<degree[0][i]; j++){
                 int u =i;
                 long long pos= nodepointer[0][i]+j;
                 int v = edges[0][pos];
                 long long reverse_pos = nodepointer[1][v]+temp_degree_count[v];
                 edges[1][reverse_pos] = u;
                 temp_degree_count[v]++; 
            }
        }
        */
		#ifdef GORDER
        load_gorder_mapping();
		#endif
        load_degree_mapping();

        return;
    }

    void split_graph(int npieces) {
        long long* splitting_points = new long long [npieces];
        memset(splitting_points, 0, sizeof(long long) * npieces);
        double frag_size = (double)m / (double) npieces;
        for (int i = 0; i < npieces; i ++) {
            splitting_points[i] = frag_size * i;
                cout << splitting_points[i] << endl; 
        }
        long * begins = new long[npieces+1];
        memset(begins, 0, sizeof(long)*(npieces+1));
        int index = 1;
        cout << "n = " << n << endl;
        FILE *split_info = fopen((data_folder+"/graph_split_info.txt").c_str(), "w");
        fprintf(split_info, "%lld\n", npieces);
        fprintf(split_info, "%lld\n", nodepointer[0][0]);
        for (int i = 0; i < n; i ++) {
            if (nodepointer[0][i] <= splitting_points[index] 
                && nodepointer[0][i+1] > splitting_points[index] ) {
                cout << nodepointer[0][i] << endl; 
                fprintf(split_info, "%lld\n", nodepointer[0][i]);
                begins[index++] = i;
            }
        }
        fflush(split_info);
        begins[npieces] = n;
        cout << "splitting points" << endl;
        for (int i = 0; i < npieces; i ++) {
            cout << begins[i] << endl;
        }
        
        cout << "--------------start------------" << get_current_time_str() << endl;
//        setWorkers(npieces);
        //#pragma cilk grainsize=1
        //cilk_for(int i = 0; i < npieces; i ++) {
        for(int i = 0; i < npieces; i ++) {
            
            FILE *stats = fopen((data_folder+"/graph_part"+to_string(i)+".txt").c_str(), "w");
            for (int node = begins[i]; node < begins[i+1]; node++) {
                for (long long offset = nodepointer[0][node]; 
                            offset < nodepointer[0][node+1]; offset ++) {
                    int t = edges[0][offset];
                    fprintf(stats, "%d %d\n", node,t);
                                    
                }
            }
            fflush(stats);
        }
        
        cout << "--------------stop------------" << get_current_time_str() << endl << endl << endl;
    }
    void snapToAdj() {
        
        // generate adj graph
        string adj_file_str = data_folder + "/graph_adj.txt";
        if(!exists_test(adj_file_str)){
            cout << "generating adjacent list" << endl;

            cout << "--------------start------------" << get_current_time_str() << endl;
            FILE *stats = fopen((data_folder+"/graph_adj.txt").c_str(), "w");
            fprintf(stats, "AdjacencyGraph\n");
            fprintf(stats, "%d\n", n);
            fprintf(stats, "%lld\n", m);
            long long offset = 0;
            for(int i=0; i<n; i++){
                fprintf(stats, "%lld\n", offset);
                offset += degree[0][i];
            }
            for (long long i = 0; i < m; i ++) {
                if (i % 100000000 == 0) {
                    cout << i <<" "<< edges[0][i] << endl;
                    fflush(stats);
                    //fsync(fileno(stats));
                }
                fprintf(stats, "%d\n", edges[0][i]);
            }
            fflush(stats);
            //fsync(fileno(stats));
            //fclose(stats);
            cout << "--------------stop------------" << get_current_time_str() << endl << endl << endl;
        }
    }
    void load_degree_mapping() {
        string degree_order_file_str = data_folder + "/order.txt";
        //config.top_heavy = int(log(n)/log(2));
        degree_order.reserve(n);
        
        for(int i=0; i<n; i++)
            degree_order.push_back(0);
        if(!exists_test(degree_order_file_str)){
            ofstream degree_file(degree_order_file_str);
            for(int i=0; i<n; i++){
                degree_pair.push_back(make_pair(i, degree[0][i]));
            }
            sort(degree_pair.begin(),degree_pair.end(), pair_sorter_large_first);
            for(int i=0; i<n; i++){
                int v = degree_pair[i].first;
                double degree = degree_pair[i].second;
                degree_order[v] = i;
                degree_file<<v<<" "<<degree<<endl;
            }
        }else{
            ifstream degree_file(degree_order_file_str);
            int v;
            double degree;
            int line=0;
            while(degree_file>>v){
                degree_order[v] = line;
                degree_file>>degree;
                degree_pair.push_back(make_pair(v,degree));
                line++;
            }
        }
        degree_reverse_map = vector<int>(n, 0);
        for (int i = 0; i < n; i ++) {
            degree_reverse_map[degree_order[i]] = i;
        }

        printf("graph degree order read\n");
        DEGREE_THRES = m/n * 3;
        return;
    }

    void load_gorder_mapping() {
        //string degree_order_file_str = data_folder + FILESEP + "order.txt";
        string gorder_file_str = data_folder + "/graph_Gmap.txt";
        //config.top_heavy = int(log(n)/log(2));
        gorder.reserve(n);
        
        for(int i=0; i<n; i++)
            gorder.push_back(0);
        if(!exists_test(gorder_file_str)){
            printf("graph_Gmap.txt not found\n");
            exit(0);
        }else{
            ifstream gorder_file(gorder_file_str);
            int v;
            int line = 0;
            while(gorder_file>>v){
                //degree_order[v] = line;
                gorder[line] = v;
                line++;
            }
        }
        reverse_map = vector<int>(n, 0);
        for (int i = 0; i < n; i ++) {
            reverse_map[gorder[i]] = i;
        }

        printf("graph_Gmap read\n");
        
    }

    double get_avg_degree() const {
        return double(m) / double(n);
    }

	bool exists_test(const std::string &name) {
	    ifstream f(name.c_str());
	    if (f.good()) {
	        f.close();
	        return true;
	    }
	    else {
	        f.close();
	        return false;
	    }
	}
    void assert_file_exist(string desc, string name) {

    	if (!exists_test(name)) {
    	    cerr << desc << " " << name << " not find " << endl;
    	    exit(1);
    	}
	} 

    // put vertices without out neighbors at the end of the edgemap
    int separate_vertices(int sep) {

        sep_mapping = new int[n];
        memset(sep_mapping, 0, sizeof(int)*n);
        if (sep == 0) {
            // no need to separate
            dividingline = n;
            #pragma simd
            for (int i = 0; i < n ; i ++) {
                sep_mapping[i] = i;
            }
            return n;
        }
        
        printf("separation begin\n");
        int sidx = 0;
        long int edgecount = 0;
        for (int i = 0; i < n; i ++) {
            if ( degree[0][i] > 0) {
                sep_mapping[i] = sidx++;
            }
            edgecount += degree[0][i];
            if (degree[0][i] < 0) {
                printf("degree[0][%d] = %d\n", i, degree[0][i]);
            }
        }
        printf("edge count is %ld\n", edgecount);
        dividingline = sidx;
        for (int i = 0; i < n; i ++) {
            if ( degree[0][i] == 0) {// && is_moved[i] == 0) {
                sep_mapping[i] = sidx++;
//                tmp_nodepointer[sidx-1] = scanned_edges;
            }
        }
 
        printf("dividingline = %lld\n", dividingline);
        // check index
        if (sidx != n) {
            printf("Separating wrong! sidx = %d n = %d\n",sidx, n);
            exit(0);
        }
        // get reverse map
        sep_reverse_map = vector<int>(n, 0);
        for (int i = 0; i < n; i ++) {
            sep_reverse_map[sep_mapping[i]] = i;
        }
        if (dividingline == n) {
            printf("dividingline == n\n");
            //return n;
        }

        // update degree[0] array thres[] array and edges[0]
        // copy results to temporal arrays
        int* tmp = new int[n];
        int* tmp_r = new int[n];
        long* tmp_np = new long[n+1];
        long* tmp_np_r = new long[n+1];
        //vector<double> tmp_thres = vector<double>(n, 0);
        for (int i = 0; i < n; i ++) {
            
            tmp[sep_mapping[i]] = degree[0][i];
            //tmp_thres[sep_mapping[i]] = thres[i];

            tmp_r[sep_mapping[i]] = degree[1][i];
            
            if (degree[0][i] == 0) {
                tmp_np[sep_mapping[i]] = m;
            } else {
                tmp_np[sep_mapping[i]] = nodepointer[0][i];
            }
            tmp_np_r[sep_mapping[i]] = nodepointer[1][i];
        }
        // write results 
        for (int i = 0; i < n; i ++) {
            degree[0][i] = tmp[i];
            degree[1][i] = tmp_r[i];
            //degree[0][i] = 0;
            //nodepointer[0][i] = tmp_nodepointer[i];
            nodepointer[0][i] = tmp_np[i];
            //nodepointer[1][i] = tmp_np_r[i];
            //thres[i] = tmp_thres[i];
        }
         
        nodepointer[0][n] = m;
        // the edge map does not have to relocate
        cilk_for (long int i = 0; i < m; i ++) {
            edges[0][i] = sep_mapping[edges[0][i]];
        }
        cout << "separation finished" << endl;        

    #ifdef CHECKING
        // check edges
        for (int i = 0; i < n; i ++) {
            int oldnode = i;
            int newnode = sep_mapping[i];
            long offset = nodepointer[0][newnode];
            long length = nodepointer[0][newnode+1] - offset;
            #ifdef OLD_VERSION_COMPARISON 
            if (length != g[i].size() ) {
                printf("out degree is not right %d %d %ld %ld\n", i,newnode, length, g[i].size());
                exit(0);
            }
            for (long j = 0; j < g[i].size(); j ++) {
                int neighbor_o = g[i][j];
                int neighbor_n = edges[0][offset + j];
                if (sep_mapping[neighbor_o] != neighbor_n) {
                    printf("separation failed i=%d j=%ld o=%d n=%d map[o] = %d\n", i, j, neighbor_o, neighbor_n, sep_mapping[neighbor_o]);
                    printf("oldnode %d newnode %d\n", i, newnode);
                    for (int ii = 0 ; ii < g[i].size(); ii ++) {
                        printf("%d (map=%d) %d\n",g[i][ii], sep_mapping[g[i][ii]], edges[0][offset+ii]);
                    }
                    exit(0);
                }
            }
            #endif
        }

        for (int i = dividingline; i <= n; i ++) {
            if (nodepointer[0][i] != m) {
                printf("failed\n");
                exit(0);
            }
        }
    #endif // CHECKING

        //cout << "checking in separate() done" << endl;
        sort_edges();
        cout << "sorting in separate() done" << endl;
        // relocate adjacent lists so that they are in right postion
        // just move the nodepointers
        delete []tmp;
        delete []tmp_r;
        delete []tmp_np;
        delete []tmp_np_r;
        return dividingline;
    }

    // NOTICE: this function has been checked
    void sort_edges() {
//        setWorkers(40);
        #pragma cilk grainsize = 8192
        cilk_for (int i = 0; i < n; i ++) {
            long long start = nodepointer[0][i];
            long long end = nodepointer[0][i+1];
            merge_sort(edges[0], start, end-1);
        }
    }

    void merge_sort(int* array, long long begin, long long end) {
        if (begin < end) {
    
            long long mid = (begin + end) / 2;
            merge_sort(array, begin, mid);
            merge_sort(array, mid+1, end);
    
            int* result = new int[end - begin+1];
            long long result_index = 0;
            long long li = begin;
            long long ri = mid+1;
            while (li <= mid && ri <= end) {
                int nodel = array[li];
                int noder = array[ri];
                if (nodel < noder) {
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
            for (long long i = 0; i < result_index; i ++) {
                array[begin+i] = result[i];
            }
            delete[] result;
        }
    }

    // Please perform separate_vertices() before calling this function
    // This operation is based on seperation
    // Randomly change the position of nodes before dividingline
    // SHARED by all set cover functions
    void gen_cover_map() {
        if (cover_map == NULL) {
            cover_map = new int[n];
        }
        memset(cover_map, 0 , sizeof(int) * n);
        if (cover_pos == NULL) {
            cover_pos = new int[n];
        }
        memset(cover_pos, 0, sizeof(int) * n);
        int* pos = new int[n];
        // only check vertices with outgoing neighbors
        for (int i = 0; i < n; i ++) {
            pos[i] = i;
            //pos[n-1-i] = degree_order[i];
            //pos[i] = degree_order[i];
        }
        // random permute
        /*
        for (int i = dividingline-1; i >= 0; i --) {
            // generate a random position
            long r = rand();
            int index = r % dividingline;
            
            if (index >= dividingline || index < 0) {
                printf("random index out of bound\n");
                exit(0);
            }
            // swap content
            int tmp = pos[i];
            pos[i] = pos[index];
            pos[index] = tmp;
        }
        */
        // get mapping and reverse mapping
        for (int i = 0; i < n; i ++) {
            cover_map[pos[i]] = i;
            cover_pos[i] = pos[i];
        }
        // check cover map (1 to 1 mapping)
        int* checked = new int[n];
        memset(checked, 0 , sizeof(int)*n);
        for (int i = 0; i< n; i ++) {
            checked[pos[i]] += 1;
        }

        // NOTICE: no duplication in map
        for (int i = 0 ; i< n; i ++) {
            if (checked[i] != 1) {
                printf("cover map check failed\n");
                exit(0);
            }
        }
        delete[] pos;
        delete[] checked;
    }

    vector<int> result_single;
    bool* is_in_sets = NULL;
    bool * is_node_covered = NULL;
    bool * is_edge_covered = NULL;
#ifdef SET_COVER
    void gen_set_cover(vector<int>& result_set, bool is_multi) {
        printf("set cover begin\n");
        if (is_in_sets == NULL) {
            is_in_sets  = new bool[n];
            memset(is_in_sets, 0, sizeof(bool) * n);
        }
        if (is_node_covered == NULL) {
            is_node_covered = new bool[n];
            memset(is_node_covered, 0, sizeof(bool) * n);
        }
        if (is_edge_covered == NULL) {
            is_edge_covered = new bool[m];
            memset(is_edge_covered, 0, sizeof(bool) * m);
        }
        gen_cover_map();
        //reset cover map
        // mark the vertex that has been checked (either in set or removed)
        int * checked = new int[n];
        memset(checked, 0, sizeof(int) * n);
        result_set = vector<int>();
        for (int i = 0; i < dividingline; i ++) {
            // get node 
            int node = cover_pos[i];
            if (checked[i] == 0 && is_in_sets[node] == 0) {
                // this vertex is not checked yet
                
                // mark its outneighbors' in neighbor
                int start = nodepointer[0][node];
                int end = nodepointer[0][node+1];
                //for (int j = 0; j < g_sep[node].size(); j ++) {
                for (int j = start; j < end; j ++) {
                    //int neighbor = g_sep[node][j];
                    int neighbor = edges[0][j];
                    for (int k = 0; k < gr_sep[neighbor].size(); k ++) {
                        int conflict = gr_sep[neighbor][k];
                        conflict = cover_map[conflict];
                        checked[conflict] = 1;
                    }
                }
                result_set.push_back(node);
                //checked[node] = 1;
                if (is_multi)
                    is_in_sets[node] = 1;
                // mark its out neighbors and edges
                for (int j = start; j < end; j ++) {
                    is_node_covered[edges[0][j]] = 1;
                    is_edge_covered[j] = 1;
                }
            }
        }

        // check workload
        int sum = 0;
        for (int v : result_set) {
            sum += g_sep[v].size();
        }
        
        printf("Set cover size is %ld\n", result_set.size());
        printf("Set cover workload is %d\n", sum);
        // chech coverage
        int node_covered = 0;
        int edge_covered = 0;
        for (int i = 0; i < n; i ++) {
            if (is_node_covered[i])
                node_covered ++;
        }
        for (int i = 0; i < m; i ++) {
            if (is_edge_covered[i] == 1) 
                edge_covered ++;
        }
        printf("%d nodes covered\n", node_covered);
        printf("%d edges covered\n", edge_covered);
        delete[] checked;
    }
    
    vector< vector<int> > multi_results;
    void gen_multi_cover(int nset) {
        printf("generate %d set covers", nset);
        multi_results = vector<vector<int>> (nset, vector<int>());
        int sum = 0;
        int size = 0;
        for (int i = 0; i < nset; i ++) {
            gen_set_cover(multi_results[i], 1);
            for (int v : multi_results[i]) {
                //int node = cover_pos[v];
                int node = v;
                sum += g_sep[node].size();
            }
            size += multi_results[i].size();
        }
        int sum_indegree = 0;
        for (int i = 0 ; i < n; i ++) {
            sum_indegree += g_sep[i].size();
        }

        int node_covered = 0;
        int edge_covered = 0;
        int pull_size = 0;
        for (int i = 0; i < n; i ++) {
            if (is_node_covered[i]) {
                node_covered ++;
            } else {
                // check pull workload
                int start = nodepointer[1][i];
                int end   = nodepointer[1][i+1];
                pull_size += end - start;
            }
        }
        for (int i = 0; i < m; i ++) {
            if (is_edge_covered[i] == 1) 
                edge_covered ++;
        }
        printf("%d nodes covered\n", node_covered);
        printf("%d edges covered\n", edge_covered);
        printf("Set covers size is %d\n", size);
        printf("Set cover total workload is %d\n", sum);
        printf("Pull total workload is %d\n", pull_size);
        printf("Graph total workload is %d\n", sum_indegree);
    }

    void check_set_covers() {
        for (int i = 0; i < rounds; i ++) {
            // scan the vector
            int * checked = new int[n];
            memset(checked, 0, sizeof(int) * n);
            for (int v = 0; v < multi_results[i].size(); v ++) {
                // get the node in set covers
                int node = multi_results[i][v];

                long long offset = nodepointer[0][node];
                long long  start = offset;
                long long end = nodepointer[0][node+1];
                for (long ii = start; ii < end; ii++) {
                    // update the residues
                    int neighbor = edges[0][ii];
                }
                for (int j : g_sep[node]) {
                    checked[j] += 1;
                }
            }

            for (int j = 0; j < n; j ++) {
                if (checked[j] > 1) {
                    printf("set cover check failed %d %d %d\n", i , j, checked[j]);
                    printf("%ld\n", gr_sep[j].size());
                    exit(0);
                }
            }
            delete[] checked;
        }
    }
#endif  //SETCOVER
    void order_by_pagerank( vector<double>& pg_values) {

        string pagerank_order_file_str = data_folder + "/pagerank_order.txt";
        
        pagerank_order.reserve(n);
       
        for(int i=0; i<n; i++)
            pagerank_order.push_back(0);
        if(!exists_test(pagerank_order_file_str)){
            ofstream pagerank_file(pagerank_order_file_str);
            for(int i=0; i<n; i++){
                pagerank_pair.push_back(make_pair(i, pg_values[i]));
            }
            sort(pagerank_pair.begin(),pagerank_pair.end(), pair_sorter_large_first);
            for(int i=0; i<n; i++){
                int v = pagerank_pair[i].first;
                double pagerank = pagerank_pair[i].second;
                pagerank_order[v] = i;
                pagerank_file<<v<<" "<<pagerank<<endl;
            }
        }else{
            ifstream pagerank_file(pagerank_order_file_str);
            int v;
            double pagerank;
            int line=0;
            printf("read pagerank from file\n");
            while(pagerank_file>>v){
                pagerank_order[v] = line;
                pagerank_file>>pagerank;
                pagerank_pair.push_back(make_pair(v,pagerank));
                line++;
            }
        }
        pr_reverse_map = vector<int>(n, 0);
        for (int i = 0; i < n; i ++) {
            pr_reverse_map[pagerank_order[i]] = i;
        }
    }
    static bool pair_sorter_large_first(pair<int, double> p1, pair<int, double> p2) {
        if (p1.second > p2.second)
            return true;
        else if (p1.second < p2.second)
            return false;
        else if (p1.first < p2.first)
            return true;
        else
            return false;
    }

#ifdef SET_COVER_PR
    void gen_pagerank_set_cover(vector<int>& result_set, bool is_multi) {

        // init
        if (is_in_sets == NULL) {
            is_in_sets  = new bool[n];
            memset(is_in_sets, 0, sizeof(bool) * n);
        }
        if (is_node_covered == NULL) {
            is_node_covered = new bool[n];
            memset(is_node_covered, 0, sizeof(bool) * n);
        }
        if (is_edge_covered == NULL) {
            is_edge_covered = new bool[m];
            memset(is_edge_covered, 0, sizeof(bool) * m);
        }
        // mark the vertex that has been checked (either in set or removed)
        int * checked = new int[n];
        memset(checked, 0, sizeof(int) * n);
        result_set = vector<int>();

        for (int i = 0; i < n; i ++) {
            // get node 
            int node = sep_mapping[pagerank_pair[i].first];
            
            //int node = sep_mapping[degree_pair[i].first];
            if (checked[node] == 0 && is_in_sets[node] == 0) {
                // this vertex is not checked yet
                // mark its outneighbors' in neighbor
                int start = nodepointer[0][node];
                int end = nodepointer[0][node+1];
                for (int j = start; j < end; j ++) {
                    int neighbor = edges[0][j];
                    for (int k = 0; k < gr_sep[neighbor].size(); k ++) {
                        int conflict = gr_sep[neighbor][k];
                        checked[conflict] = 1;
                    }
                }
                for (int k = 0; k < gr_sep[node].size(); k ++) {
                    int conflict = gr_sep[node][k];
                    checked[conflict] = 1;
                }
                // avoid data race in propagation
                //checked[node] = 1;
                result_set.push_back(node);
                //checked[node] = 1;
                if (is_multi)
                    is_in_sets[node] = 1;
                // mark its out neighbors and edges
                for (int j = start; j < end; j ++) {
                    is_node_covered[edges[0][j]] = 1;
                    is_edge_covered[j] = 1;
                }
            }
        }

        // check workload
        int sum = 0;
        for (int v : result_set) {
            sum += g[v].size();
        }
        
        printf("Set cover size is %ld\n", result_set.size());
        printf("Set cover workload is %d\n", sum);
        delete[] checked;
    }
    vector< vector<int> > pr_multi_results;
    void gen_pagerank_multi_cover(int nset) {
        rounds = nset;
        bool* is_in_sets = new bool[n];
        pr_multi_results = vector<vector<int>> (nset, vector<int>());
        int sum = 0;
        int size = 0;
        for (int i = 0; i < nset; i ++) {
            gen_pagerank_set_cover(pr_multi_results[i], 1);
            for (int v : pr_multi_results[i]) {
                //int node = cover_pos[v];
                int node = v;
                sum += g[node].size();
            }
            size += pr_multi_results[i].size();
        }
        int sum_indegree = 0;
        for (int i = 0 ; i < n; i ++) {
            sum_indegree += g[i].size();
        }
        printf("Set covers size is %d\n", size);
        printf("Set cover total workload is %d\n", sum);
        printf("Graph total workload is %d\n", sum_indegree);

        int node_covered = 0;
        int edge_covered = 0;
        int pull_size = 0;
        for (int i = 0; i < n; i ++) {
            if (is_node_covered[i]) {
                node_covered ++;
            } else {
                // check pull workload
                int start = nodepointer[1][i];
                int end   = nodepointer[1][i+1];
                pull_size += end - start;
            }
        }
        for (int i = 0; i < m; i ++) {
            if (is_edge_covered[i] == 1) 
                edge_covered ++;
        }
        printf("%d nodes covered\n", node_covered);
        printf("%d edges covered\n", edge_covered);
        printf("Set covers size is %d\n", size);
        printf("Set cover total workload is %d\n", sum);
        printf("Pull total workload is %d\n", pull_size);
        printf("Graph total workload is %d\n", sum_indegree);

        delete[] is_in_sets;
    }
    void check_pr_set_covers() {
        for (int i = 0; i < rounds; i ++) {
            // scan the vector
            int * checked = new int[n];
            memset(checked, 0, sizeof(int) * n);
            for (int v = 0; v < pr_multi_results[i].size(); v ++) {
                // get the node in set covers
                int node = pr_multi_results[i][v];

                long long offset = nodepointer[0][node];
                long long  start = offset;
                long long end = nodepointer[0][node+1];
                for (long ii = start; ii < end; ii++) {
                    // update the residues
                    int neighbor = edges[0][ii];
                }
                for (int j : g_sep[node]) {
                    checked[j] += 1;
                }
            }

            for (int j = 0; j < n; j ++) {
                if (checked[j] > 1) {
                    printf("set cover check failed %d %d %d\n", i , j, checked[j]);
                    printf("%ld\n", gr_sep[j].size());
                    exit(0);
                }
            }

            delete[] checked;
        }
    }
#endif
    //========================= set cover with duplicates =====================
#ifdef SET_COVER_DUP 
    void gen_one_set_cover_dup(vector<int>& result_set, bool is_multi) {
        //printf("set cover with duplicates begin\n");
        if (is_in_sets == NULL) {
            is_in_sets  = new bool[n];
            memset(is_in_sets, 0, sizeof(bool) * n);
        }
        if (is_node_covered == NULL) {
            is_node_covered = new bool[n];
            memset(is_node_covered, 0, sizeof(bool) * n);
        }
        if (is_edge_covered == NULL) {
            is_edge_covered = new bool[m];
            memset(is_edge_covered, 0, sizeof(bool) * m);
        }
        // TODO ??? what is gen_cover_map() used for???
        gen_cover_map();
        //reset cover map
        // mark the vertex that has been checked (either in set or removed)
        int * checked = new int[n];
        memset(checked, 0, sizeof(int) * n);
        result_set = vector<int>();
        for (int i = 0; i < dividingline; i ++) {
            // get node 
            int node = cover_pos[i];
            if (checked[i] == 0 && is_in_sets[node] == 0) {
                // this vertex is not checked yet
                
                // mark its outneighbors' in neighbor
                int start = nodepointer[0][node];
                int end = nodepointer[0][node+1];
                //for (int j = 0; j < g_sep[node].size(); j ++) {
                for (int j = start; j < end; j ++) {
                    //int neighbor = g_sep[node][j];
                    int neighbor = edges[0][j];
                    for (int k = 0; k < gr_sep[neighbor].size(); k ++) {
                        int conflict = gr_sep[neighbor][k];
                        conflict = cover_map[conflict];
                        checked[conflict] = 1;
                    }
                }
                result_set.push_back(node);
                //checked[node] = 1;
                if (is_multi)
                    is_in_sets[node] = 1;
                // mark its out neighbors and edges
                for (int j = start; j < end; j ++) {
                    is_node_covered[edges[0][j]] = 1;
                    is_edge_covered[j] = 1;
                }
            }
        }

        // check workload
        int sum = 0;
        for (int v : result_set) {
            sum += g_sep[v].size();
        }
        
        printf("Set cover size is %ld\n", result_set.size());
//        printf("Set cover workload is %d\n", sum);
        // chech coverage
        int node_covered = 0;
        int edge_covered = 0;
        for (int i = 0; i < n; i ++) {
            if (is_node_covered[i])
                node_covered ++;
        }
        for (int i = 0; i < m; i ++) {
            if (is_edge_covered[i] == 1) 
                edge_covered ++;
        }
//        printf("%d nodes covered\n", node_covered);
//        printf("%d edges covered\n", edge_covered);
        delete[] checked;
    }
    
    vector< vector<int> > multi_results_dup;
    void gen_multi_cover_dup(int nset) {
        printf("generate %d set covers with duplicates\n", nset);
        multi_results_dup = vector<vector<int>> (nset, vector<int>());
        int sum = 0;
        int size = 0;
        for (int i = 0; i < nset; i ++) {
            gen_one_set_cover_dup(multi_results_dup[i], 1);
            for (int v : multi_results_dup[i]) {
                //int node = cover_pos[v];
                int node = v;
                sum += g_sep[node].size();
            }
            size += multi_results_dup[i].size();
        }
        int sum_indegree = 0;
        for (int i = 0 ; i < n; i ++) {
            sum_indegree += g_sep[i].size();
        }

        int node_covered = 0;
        int edge_covered = 0;
        int pull_size = 0;
        for (int i = 0; i < n; i ++) {
            if (is_node_covered[i]) {
                node_covered ++;
            } else {
                // check pull workload
                int start = nodepointer[1][i];
                int end   = nodepointer[1][i+1];
                pull_size += end - start;
            }
        }
        for (int i = 0; i < m; i ++) {
            if (is_edge_covered[i] == 1) 
                edge_covered ++;
        }

        // check the rest of nodes
        int workload_remain = 0;
        int nodes_in_sets = 0;
        for (int i = 0; i < dividingline; i ++) {
            if (!is_in_sets[i]) {
                workload_remain += get_degree(i); 
            } else {
                nodes_in_sets ++;
            }
        }
        printf("%d nodes covered\n", node_covered);
        printf("%d edges covered\n", edge_covered);
        printf("Set covers size is %d\n", size);
        printf("Set cover total workload is %d\n", sum);
        printf("Pull total workload is %d\n", pull_size);
        printf("Workload remain: %d nodes %d edges avg=%d\n", 
                dividingline - nodes_in_sets, workload_remain,
                workload_remain/(dividingline-nodes_in_sets));
        printf("Graph total workload is %d\n", sum_indegree);

        // check remain nodes coverage

        int* remain_coverage = new int[n];
        memset(remain_coverage, 0, sizeof(int) * n );
        for (int i = 0; i < dividingline; i ++) {
            if (!is_in_sets[i]) {
                int start = nodepointer[0][i];
                int end = nodepointer[0][i+1];
                for (int j = start; j < end; j ++) {
                    int neighbor = edges[0][j];
                    remain_coverage[neighbor] = 1;
                }
            }
        }
        int sum_coverage = 0;
        for (int i = 0; i < n; i ++) {
            if (remain_coverage[i] == 1)
                sum_coverage ++;
        }
        pull_size = 0;
        for (int i = 0; i < n; i ++) {
            if (remain_coverage[i]) {
                // check pull workload
                int start = nodepointer[1][i];
                int end   = nodepointer[1][i+1];
                pull_size += end - start;
            }
        }
        printf("%d nodes covered by remaining nodes, pull workload is %d\n", 
                    sum_coverage, pull_size);
        delete[] remain_coverage;
    }

    void check_set_covers_dup() {
        for (int i = 0; i < rounds; i ++) {
            // scan the vector
            int * checked = new int[n];
            memset(checked, 0, sizeof(int) * n);
            for (int v = 0; v < multi_results_dup[i].size(); v ++) {
                // get the node in set covers
                int node = multi_results_dup[i][v];

                long long offset = nodepointer[0][node];
                long long  start = offset;
                long long end = nodepointer[0][node+1];
                for (long ii = start; ii < end; ii++) {
                    // update the residues
                    int neighbor = edges[0][ii];
                }
                for (int j : g_sep[node]) {
                    checked[j] += 1;
                }
            }

            for (int j = 0; j < n; j ++) {
                if (checked[j] > 1) {
                    printf("set cover check failed %d %d %d\n", i , j, checked[j]);
                    printf("%ld\n", gr_sep[j].size());
                    exit(0);
                }
            }
        }
    }
#endif
    //==================== END OF SETCOVER WITH DUPLICATES ====================
    // workload scheduling: move nodes with most out degrees to critical positions
    

    // random walk indexing
    // This process happens after all other initializations
    
    // random numbers
    double*   random_probs=NULL;
    int*      random_neighbors=NULL;
    long int* random_pointers=NULL;
    
    // results
    long long int* dest_pointers=NULL;  // calculated in generate_random_walks()
    int*      destinations=NULL;

    long long int get_dest_start(int i) {return dest_pointers[i];}
    long long int get_dest_size(int i) {return dest_pointers[i+1]-dest_pointers[i];}

    inline static unsigned long lrand() {
        return rand();
    }
    
    inline static double drand(){
    	return rand()*1.0f/RAND_MAX;
    }

    void generate_permuted_random_walk_index() {
        
        printf("generating permuted random walk destinations thres type is %d\n"
                , thres_type);
        
        dest_pointers = new long long int [n+1];
        memset(dest_pointers, 0, sizeof(long long int) * (n+1));

        // random walk result indices
        long long int result_idx = 0;
        for (int i = 0; i < n; i ++) {
            dest_pointers[i] = result_idx;
            int node = permute_mapping_r[i];
            //int node = i;
            long long int degree_i = get_degree(node);
            long long int num_random_walk = (thres_type==0) ? ceil(rmax*degree_i*Omega) : degree_i;
            num_random_walk*= threshold_scale;
            result_idx += num_random_walk;
        }
        dest_pointers[n] = result_idx;
        printf("destinations length is %lld\n",result_idx);
        destinations = new int[result_idx];
        // random walk
        for (int i = 0 ; i < n; i ++) {
            long long int dest_idx_start = dest_pointers[i];
            long long int num_random_walk = dest_pointers[i+1] - dest_pointers[i];
            //int start_node = permute_mapping_r[i];
            int start_node = i;
            for (long long int j = 0; j < num_random_walk; j ++) {
                int dest = random_walk(start_node);
                destinations[dest_idx_start+j] = dest;        
            }
        }
        printf("total destinations number is %lld\n", result_idx);
        return;
    }
    // create random walk index (load from file if pre-stored)
    // this should be called after separation
    // The index is based on separation
    void generate_random_walks() {
        
        printf("generating random walk destinations thres type is %d\n"
                , thres_type);

        dest_pointers = new long long int [n+1];
        memset(dest_pointers, 0, sizeof(long long int) * (n+1));

		// THIS IS A CRITICAL SECTION : used for test of indexed random walk with thres of rmax
		//thres_type = 1;
		//threshold_scale = 6;
		// END OF CRITICAL SECTION
        // random walk result indices
        //printf("threshold scale in generate_random_walks() is %f\n", threshold_scale);
        cout << "threshold scale in generate_random_walks() is " << threshold_scale << endl;
        long long int result_idx = 0;
        for (int i = 0; i < n; i ++) {
            dest_pointers[i] = result_idx;
            //int node = permute_mapping_r[i];
            #ifdef HOT_SET
            int node = hot_set_map[i];
            #else
            int node = i;
            #endif
            long int degree_i = get_degree(node);
            long long int num_random_walk = (thres_type==0) ? ceil(rmax*degree_i*Omega) : degree_i;
            num_random_walk*= threshold_scale;
            if (num_random_walk != 0) {
                num_random_walk += 2;
            }
            result_idx += num_random_walk;
        }
		// THIS IS A CRITICAL SECTION : used for test of indexed random walk with thres of rmax
		//thres_type = 1;
		//threshold_scale = 1;
		// END OF CRITICAL SECTION
        dest_pointers[n] = result_idx;

        printf("result_idx = %lld\n", result_idx);
        // random walk
        load_random_walk_index();
        // check targets
        //sort_random_walk_index();
        printf("total destinations number is %lld\n", result_idx);
        //generate_graph();

        // consider popuparity
        #ifdef INTEGER_RWUPDATE
        //map_rwidx_pop();
        //map_rwidx_degree();
        //generate_random_map();
        //map_rwidx_random();
        
        #ifndef UINT8_T_RWUPDATE
        load_target_map();
        //#if !defined(FRIENDSTER)
        if (target_divideline != 0)
            map_rwidx_target();
        #endif // UINT8_T_RWUPDATE
        //#endif
        //separate_random_walk_sources();
        #endif

        #ifdef GORDER
        //map_destinations_by_gorder();
        #endif
        //exit(0);
        return;
    }
    static void setWorkers(int n) {
      __cilkrts_end_cilk();
      //__cilkrts_init();
      std::stringstream ss; ss << n;
      if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
        std::cerr << "failed to set worker count!" << std::endl;
        cout << ss.str().c_str() << endl;
        std::abort();
      }
    }
    void check_index_targets() {
        
        atomic<int>* counts = new atomic<int>[n];
        memset(counts, 0, sizeof(atomic<int>)*n);

//        setWorkers(40);
        cilk_for(int i = 0; i < n; i ++) {
            int node = i;
            long long int start_idx = dest_pointers[node];
            long long int end_idx = dest_pointers[node+1];
            for (long long int j = start_idx; j < end_idx; j ++) {
                int t = destinations[j];
                counts[t] ++;
            }
        }
        
        vector<pair<int, double> > count_pair; // residue count pairs
        for (int i = 0; i< n; i ++) {
            count_pair.push_back(make_pair(i, counts[i].load()));
        }
        // sort
        sort( count_pair.begin(), count_pair.end(), pair_sorter_large_first);
        cout << "sorted"<< endl;
        FILE *stats = fopen((data_folder+"/rwidx3_counts.txt").c_str(), "w");
        fprintf(stats, "#%s\n", data_folder.c_str());
        cout << "writing to file..."<< endl;
        for (int i = 0; i < n; i++) {
            int node = count_pair[i].first;
            fprintf(stats, "%d %d %f %ld %ld\n", i,  node, count_pair[i].second, 
                                //pagerank_map[rc_pair[i].first],
                                get_degree(node), get_in_degree(node));
        }
        cout << "written to file"<< endl;
        fprintf(stats, "e\n");
        fflush(stats);
        printf("random index counts checked\n");
        delete[] counts;
        //exit(0);
    }
	void set_worker_number(int number) {

        std::stringstream ss; ss << number;
        if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
          std::cerr << "failed to set worker count!" << std::endl;
            cout << ss.str().c_str() << endl;
          std::abort();
        }
	}
    inline int random_walk(int start ){
        int cur = start;
        unsigned long k;

        while (true) {
            #ifdef USE_DOTMIX 
			double random_num = (double)global_rng.get()/UINTMAX_MAX;
			#else
			double random_num = drand();	
			#endif
            if ( random_num < alpha) {
                return cur;
            }
            int degree = get_degree(cur);
            if (degree){
				#ifdef USE_DOTMIX
				k = global_rng.get() % degree;
				#else
                k = rand() % degree;
				#endif
                cur = edges[0][nodepointer[0][cur] + k];
            } else {
                cur = start;
            }
        }
    }
    
    // =======================PERMUTATION FOR RANDOMWALK================
    // BEGIN OF random walk load balance
    // return largest denominator
    int generate_rw_locations(int* locs, int size) {
        int iter = 0;
        int index = 0;
        locs[index++] = 0;
        for (int i = 0; i < size; i ++) {
            int num_in_layer = pow(2,i);
            int denominator = pow(2,i+1);
            for (int j = 0; j < num_in_layer; j ++) {
                int numerator = 2 * j + 1;
                locs[index] = (int)((long long int)dividingline * numerator / denominator) + 2;
                //printf("%d %d\n", numerator, denominator);
                index ++;
                if (index == size) {
                    for (int k = 0; k < size; k ++) {
                        //printf("locations[%d] = %d\n", k, locs[k]);
                    }
                    return denominator;
                }
            }
        }
        return -1;
    }
    int *is_occupied=NULL;
    int *permute_mapping=NULL;
    int *permute_mapping_r=NULL;
    // assigh each node a new location
    void permute_nodes_for_rw() {
        // init        
        permute_mapping_r = new int[dividingline];
        memset(permute_mapping_r, -1, sizeof(int)*dividingline);

        is_occupied = new int[dividingline]; // mark the location that is occupied by a node
        permute_mapping = new int[n];
        memset(is_occupied, -1, sizeof(int)*dividingline);
        memset(permute_mapping, -1, sizeof(int)*n);
        int* locations = new int[dividingline];
        int* node_is_permuted = new int[n];
        memset(node_is_permuted, -1, sizeof(int) * n);
        // end of init
        // merge_sort(nodes_array, 0, list_size);
        int k = 40;// number of cores
        int num_locs = k*4;
        int denominator= generate_rw_locations(locations, num_locs)*2;
        int segment_size = n/denominator - denominator;
        int counter = 0;
        for (int i = 0; i < dividingline; i ++) {
            
            int node = degree_pair[i].first; 
            node = sep_mapping[node];
            int niter = i / num_locs;
            if (niter > segment_size -1) break;
            int loc_i = locations[i % num_locs] + niter;
            if (loc_i >= dividingline) {
                printf("location %d is illegal i=%d node=%d /=%d...%d\n",loc_i,i,node,niter,i%num_locs);
                exit(0);
            }
            if (is_occupied[loc_i] > 0) {
                printf("location %d is occupied i=%d node=%d /=%d...%d\n",loc_i,i,node,niter,i%num_locs);
                exit(0);
            }
            permute_mapping[node] = loc_i;
            permute_mapping_r[loc_i] = node;
            is_occupied[loc_i] = 1;
            node_is_permuted[node] = 1;
            counter++;
        }
        #ifdef DEBUGINFO
        printf("counter = %d, dividingline - counter = %d\n", counter, dividingline - counter);
        #endif
        // permute nodes whose location is not taken (for memory locality)
        for (int i = 0; i < dividingline; i ++) {
            int node = i;
            // check is this node is mapped
            if (node_is_permuted[node] == -1) {
                // put it in orginial location if possible
                if (is_occupied[i] == -1) {
                    permute_mapping[node] = i;
                    permute_mapping_r[i] = node;
                    is_occupied[i] = 1;
                    node_is_permuted[node] = 1;
                    counter ++;
                }
            }
        }
        #ifdef DEBUGINFO
        printf("counter = %d, dividingline - counter = %d\n", counter, dividingline - counter);
        #endif
        // permute nodes whose location is taken
        int vaccancy_index = 0;
        for (int i = 0; i < dividingline; i ++) {
            int node = i;
            if (node_is_permuted[node] == -1) {
                if (is_occupied[i] == -1) {
                    printf("wrong!!!!!!!!!\n");
                    exit(0);
                }
                while(is_occupied[vaccancy_index] == 1) 
                    vaccancy_index++;

                permute_mapping[node] = vaccancy_index;
                permute_mapping_r[vaccancy_index] = node;
                is_occupied[vaccancy_index] = 1;
                node_is_permuted[node] = 1;
                vaccancy_index++;
                counter++;
            }
        }
        #ifdef DEBUGINFO
        printf("counter = %d, dividingline - counter = %d\n", counter, dividingline - counter);
        #endif
        // checking
        int count_zero = 0;
        for (int i = 0 ; i < dividingline; i ++) {
            if (is_occupied[i] > 1) {
                printf("permutation failed is_occupied[%d]=%d\n", i, is_occupied[i]);
                exit(0);
            }
            if (is_occupied[i] != 1) count_zero++;
        }
        #ifdef DEBUGINFO
        for (int i = 0; i < 10; i ++) {
            
            int node = degree_pair[i].first; 
            node = sep_mapping[node];
            printf("%f\n", (double)permute_mapping[node]/dividingline);
        }
        #endif
        int* counters = new int[n];
        memset(counters, 0, sizeof(int)*n);
        for (int i = 0 ; i < dividingline; i ++) {
            
            counters[permute_mapping_r[i]] ++;
        }
        for (int i = 0 ; i < n; i ++) {
            if (counters[i] > 1) {
                printf("counters[%d] = %d\n", i, counters[i]);
                exit(0);
            }
        }
        delete[] locations;
        delete[] node_is_permuted;
        delete[] counters;
    }
	
    void random_walking() {
        long long dest_count=0;
        struct timespec start, finish;
        double elapsed=0.0;
        
        setWorkers(sysconf(_SC_NPROCESSORS_ONLN));
        
        destinations = new int[dest_pointers[n]];
        //ofstream index_file(rw_index_file_str);
        clock_gettime(CLOCK_REALTIME, &start);
		#ifdef USE_DOTMIX
        #pragma cilk grainsize=1
        cilk_for (int i = 0 ; i < dividingline; i ++) {
		#else
        for (int i = 0 ; i < dividingline; i ++) {
		#endif
        
            long long int dest_start = dest_pointers[i];
            long long int num_random_walk = dest_pointers[i+1] - dest_pointers[i];
            #ifdef RWINDEX_PERMUTE
            int start_node = permute_mapping_r[i];
            #else
            int start_node = i;
            #endif
            for (long long int j = 0; j < num_random_walk; j ++) {
                int dest = random_walk(start_node);
                destinations[dest_start+j] = dest;        
            }
            dest_count += num_random_walk;
            //if (i % 1000000 == 0) printf("%d nodes indexed...%lld\n", i, dest_count);
        }
        //printf("number of destinations is %lld\n", dest_count);
        clock_gettime(CLOCK_REALTIME, &finish);
        elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
        elapsed += (finish.tv_nsec - start.tv_nsec);
        printf("index generated in %.3f milliseconds using %d workers.\n",
                   elapsed/1000000, __cilkrts_get_nworkers());
    }
	// generate random walk index if not generated
    // only relocate index, do not map nodes to peromuted id
    void load_random_walk_index() {
        #ifdef FAST_INDEX
        //deserialize_rw_idx();
        random_walking();
        return;
        #endif
        string rw_index_file_str = data_folder + "/rwIndexPafora3.txt";
        long long dest_count=0;
        struct timespec start, finish;
        double elapsed=0.0;
        if(!exists_test(rw_index_file_str)){
            printf("...index not found, generating index...\n");
            destinations = new int[dest_pointers[n]];
//			set_worker_number(40);
            //ofstream index_file(rw_index_file_str);
            clock_gettime(CLOCK_REALTIME, &start);
			#ifdef USE_DOTMIX
            cilk_for (int i = 0 ; i < dividingline; i ++) {
			#else
            for (int i = 0 ; i < dividingline; i ++) {
			#endif
            
                long long int dest_start = dest_pointers[i];
                long long int num_random_walk = dest_pointers[i+1] - dest_pointers[i];
                #ifdef RWINDEX_PERMUTE
                int start_node = permute_mapping_r[i];
                #else
                int start_node = i;
                #endif
                for (long long int j = 0; j < num_random_walk; j ++) {
                    int dest = random_walk(start_node);
                    destinations[dest_start+j] = dest;        
                }
                dest_count += num_random_walk;
                if (i % 1000000 == 0) printf("%d nodes indexed...%lld\n", i, dest_count);
            }
            //printf("number of destinations is %lld\n", dest_count);
            clock_gettime(CLOCK_REALTIME, &finish);
            elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
            elapsed += (finish.tv_nsec - start.tv_nsec);
            printf("index generated in %.3f milliseconds using %d workers.\n",
                   elapsed/1000000, __cilkrts_get_nworkers());

            clock_gettime(CLOCK_REALTIME, &start);
            FILE *index_file = fopen(rw_index_file_str.c_str(), "w");
            dest_count=0; 
            for(int i=0; i<dividingline; i++){
                long long int dest_start = dest_pointers[i];
                long long int num_random_walk = dest_pointers[i+1] - dest_pointers[i];
                for (long long int j = 0; j < num_random_walk; j ++) {
                    int dest = destinations[dest_start+j];        
                    //index_file << dest <<endl;
                    fprintf(index_file, "%d\n", dest);
                }
                dest_count += num_random_walk;
                if (i % 1000000 == 0) printf("%d nodes index written...\n", i);
            }
            fflush(index_file);
            fclose(index_file);
            clock_gettime(CLOCK_REALTIME, &finish);
            elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
            elapsed += (finish.tv_nsec - start.tv_nsec);
            printf("index written in %.3f milliseconds using %d workers.\n",
                   elapsed/1000000, __cilkrts_get_nworkers());
            exit(0);
        }else{
            #ifdef PARALLEL_LOAD_IDX
            rw_index_file_str = data_folder + "/rwidx3_adj.txt";
            if(!exists_test(rw_index_file_str)){
            #endif
                rw_index_file_str = data_folder + "/rwIndexPafora3.txt";
			    printf("reading index file...\n");
                ifstream index_file(rw_index_file_str);
                int v;
                long long line=0;
                while(index_file>>v){
                    destinations[line++] = v;
                }
                printf("%lld lines read\n", line);
                destinations_length = line;
            #ifdef PARALLEL_LOAD_IDX
                // output to adj list with pointer info
                FILE *adjIdx= fopen((data_folder+"/rwidx3_adj.txt").c_str(), "w");
                fprintf(adjIdx, "AdjacencyGraph\n");
                fprintf(adjIdx, "%d\n", n);
                fprintf(adjIdx, "%lld\n",dest_pointers[n]); // scale is 4
                for (int i = 0; i < n; i++) {
                    long long pointer_i = dest_pointers[i];
                    fprintf(adjIdx, "%lld\n", pointer_i);
                }
                for (long long i = 0; i < line; i ++) {
                    fprintf(adjIdx, "%d\n", destinations[i]);
                    if (i % 100000000 == 0) {
                        printf("%lld lines wirtten\n", i);
                    }
                }
                fflush(adjIdx);
                fclose(adjIdx);
                printf("%lld index edges written\n", line);
                exit(0);
            } else {
                // load
                //map_file_to_memory(rw_index_file_str);
                paidx = readGraphFromFile<int>(rw_index_file_str.c_str());
                printf("ajd idx m = %lld\n", paidx.m);
                printf("adj idx read! \n");
                destinations = paidx.edges;
                destinations_length = paidx.m;
                //exit(0);
                //check
            
                /*
                rw_index_file_str = data_folder + "/rwIndexPafora4.txt";
			    printf("reading index file...\n");
                ifstream index_file(rw_index_file_str);
                int v;
                long long line=0;
                while(index_file>>v){
                    destinations[line++] = v;
                }
                printf("%lld lines read\n", line);
                */
            }
            #endif
        }
    }

    void sort_random_walk_index() {
        std::stringstream ss; ss << 80;
        if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
          std::cerr << "failed to set worker count!" << std::endl;
          std::abort();
        }
        cilk_for (int i = 0; i < dividingline; i ++) {
            merge_sort(destinations, dest_pointers[i], dest_pointers[i+1]-1);
        }
    }

    //######################### Gorder for random walk destinations
    vector<vector<int> > g_dest;
    void generate_graph() {
        g_dest = vector<vector<int>> (n, vector<int>());
        // generate graphs
        for (int i = 0; i < dividingline; i ++) {
            long long int dest_start = dest_pointers[i];
            long long int num_random_walk = dest_pointers[i+1] - dest_pointers[i];
            for (long long int j = 0; j < num_random_walk; j ++) {
                int dest = destinations[dest_start+j];        
                g_dest[i].push_back(dest);
            }
        }
        // sort
        for (int i = 0; i < dividingline; i ++) {
            sort(g_dest[i].begin(), g_dest[i].end());
        }
        //output
        string dest_graph_file_str = data_folder + "/dest_graph.txt";
        ofstream dest_graph_file(dest_graph_file_str);
        for (int i = 0; i < dividingline; i ++) {
            int v = i;
            int last_t;
            int length = g_dest[i].size();
            if (length > 0) {
                last_t = g_dest[i][0];
                for (int j = 0; j < length; j ++) {
                    int t = g_dest[i][j];
                    if (t != last_t || j == 0) {
                        // output
                        dest_graph_file << v << " " << t << endl;
                    }
                    last_t = t;
                }
            }
            if (i % 500000 == 0) {
                printf("%d nodes processed\n",i);
            }
        }
    }

    // use gmap to map the random walk index
    // Map the nodepointers first
    void map_destinations_by_gorder() {
        // load destination gorder map seems unnecessary!
        //load_destination_gorder_map();

//        gmap_destinations = new int[destinations_length];
//        gmap_dest_pointers = new long long[n];
//        // calculate original dest pointers
//        long long* o_dest_pointers = new long long[n+1];
//        o_dest_pointers[0] = 0;
//        for (int i = 0; i < n; i ++) {
//            int g_node = gorder[i];
//            long long degree = dest_pointers[g_node+1] - dest_pointers[g_node];
//            if (degree < 0) {
////                printf("degree < 0, is %lld\n", degree);
//            }
//            o_dest_pointers[i+1] = o_dest_pointers[i] + degree;
//        }
//        if (o_dest_pointers[n] != destinations_length) {
////            printf("length wrong\n");
//        }
////        printf("destinations_length is %lld\n", destinations_length);
//        o_dest_pointers[n] = destinations_length;
/////        printf("last degree is %lld n = %d\n", o_dest_pointers[n] - o_dest_pointers[n-1],n);
//        // map index
////        setWorkers(40);
////        printf("%lld %lld \n", o_dest_pointers[1632803], o_dest_pointers[1632802]);
//        //for (int i = 0; i < n; i ++) {
//        cilk_for (int i = 0; i < n; i ++) {
//            int g_node = gorder[i];
//            long long g_start = dest_pointers[g_node];
//            long long g_end   = dest_pointers[g_node+1];
//            long long start   = o_dest_pointers[i];
//            long long end     = o_dest_pointers[i+1];
//            if (i == n-1) {
////                printf("%lld %lld i = %d\n", o_dest_pointers[1632803], o_dest_pointers[1632802],i);
//                
//            }
////            if (g_end-g_start != end -start) {
////                printf("o_dest_pointers wrong!\n");
////                printf("g_start %lld g_end %lld start %lld end %lld\n",
////                        g_start, g_end, start, end);
////                printf("%lld %lld i = %d, g_node = %d\n", g_end-g_start, end- start,i,g_node);
////                exit(0);
////            }
//            long long g_j;
//            long long j;
//            for (g_j = g_start, j = start;  j < end; g_j ++, j ++) {
//                gmap_destinations[g_j] = gorder[destinations[j]];
//            }
//        }
//        #pragma simd
//        cilk_for(long long i = 0; i < destinations_length; i ++) {
////            destinations[i] = gmap_destinations[i];
//        }
//        return;
        // output sorted gordered map
        //string output_file_str = data_folder + "/dest_sorted_gorder.txt";
        //ofstream output_file(output_file_str);
        //for (int i = 0; i < dividingline; i ++) {
        //    output_file << gmap_destinations[i] << endl;
        //}

    }
    vector<int> dest_gmap;
    vector<int> reverse_d_gmap; // reverse destination gmap
    long long destinations_length=0;
    int *gmap_destinations=NULL;
    long long *gmap_dest_pointers=NULL;

    // this function seems unnecessary
    void load_destination_gorder_map() {
        //string degree_order_file_str = data_folder + FILESEP + "order.txt";
        string dest_gmap_file_str = data_folder + "/dest_graph_Gmap.txt";
        //config.top_heavy = int(log(n)/log(2));
        dest_gmap.reserve(n);
        
        for(int i=0; i<n; i++)
            dest_gmap.push_back(0);
        if(!exists_test(dest_gmap_file_str)){
            printf("dest_graph_Gmap.txt not found\n");
            exit(0);
        }else{
            ifstream dest_gmap_file(dest_gmap_file_str);
            int v;
            int line = 0;
            while(dest_gmap_file>>v){
                //degree_order[v] = line;
                dest_gmap[line] = v;
                line++;
            }
        }
        reverse_d_gmap = vector<int>(n, 0);
        for (int i = 0; i < n; i ++) {
            reverse_d_gmap[dest_gmap[i]] = i;
        }

        printf("dest_graph_Gmap read\n");
        
    }
    //---------------------GORDER FOR RANDOM WALK DESTINATIONS--------------END
    
    //===================== HOT SET SCHEDULEING ===================
    int* hot_set;
    int hot_set_length;
    int hs_block_size;
    vector<int> hot_set_map;
    vector<int> hot_set_map_r;

    void load_hot_set() {
        
        string hot_set_file_str = data_folder + "/hot_nodes_stats.txt";
        
        hot_set_map = vector<int>(n, 0);
        hot_set_map_r = vector<int>(n, 0);

        if (!exists_test(hot_set_file_str)) {
            printf("%s not found\n", hot_set_file_str.c_str());
            exit(0);
        } else {
            FILE* hot_set_file = fopen(hot_set_file_str.c_str(), "r");
            char line[200];
            char* firstline;
            fgets(line, sizeof line, hot_set_file);

            while (fgets(line, sizeof line, hot_set_file) != NULL) {
                int node, idx;
                double hit_i;
                long int din_i, dout_i;
                sscanf(line, "%d %d %f %ld %ld\n", &idx, &node, &hit_i, &dout_i, &din_i);
                hot_set_map[node] = idx;
                hot_set_map_r[idx] = node;
            }
            fclose(hot_set_file);
        }
        printf("hot set file %s loaded\n", hot_set_file_str.c_str());
        // checking loaded file
        //for (int i = 0; i < 10; i ++) {printf("%d \n", hot_set_map_r[n-1-i]);}
        //exit(0);
    }

    void init_hot_set_block_cacheline(int ncore) {
        block_hit = vector<vector<int>>(NBLOCKS, vector<int>());

        if (indices == NULL) indices = new atomic<double>[ncore];
        memset(indices, 0, sizeof(atomic<double>)*ncore);

        if (block_selected == NULL) block_selected = new int[NBLOCKS];
        memset(block_selected, -1, sizeof(int)*NBLOCKS);

        if (nrounds == NULL) {
            nrounds = new int[ncore];
            for (int i = 0; i < ncore; i ++) {
                if (ncore*(NROUNDS-1)+i < NBLOCKS) nrounds[i] = NROUNDS;
                else nrounds[i] = NROUNDS - 1;
            }
        }
    }

    void init_hot_set_block_matrix(int ncore) {
        if (grain_blocks == NULL) {
            grain_blocks = new int*[NROUNDS];
            for (int i = 0; i < NROUNDS; i ++) {
                grain_blocks[i] = new int[ncore];
                memset(grain_blocks[i], -1, sizeof(int)*ncore);
            }
        }
    }
/*
    void compute_hot_set_block_hits() {
        // traverse the graph
        printf("grain size in compute_hot_set_block_hits() is %d\n", grain_size);
        for (int i = 0; i < n; i ++) {
            long start = nodepointer[0][i];
            long end   = nodepointer[0][i+1];
            int block_id = i / grain_size;
            // check each out neighbor
            vector<int> candidates = vector<int>(); // duplications exist in this vector
            candidates.push_back(-1); // to avoid boundary check
            
            for (long j = start; j < end; j ++) {
                int node = edges[0][j];
                int cacheline_id = node/(CACHELINE_SIZE/sizeof(double));
                candidates.push_back(cacheline_id);
            }

            for (int j = 1; j < candidates.size(); j ++) {
                int cacheline_id = candidates[j];
                // check duplicates in each node's out neighbors
                if (cacheline_id != candidates[j-1]) {
                    block_hit[block_id].push_back(cacheline_id);
                }
           }
        }
        // remove duplicates in each block
        for (int i = 0; i < NBLOCKS; i ++) {
            sort(block_hit[i].begin(), block_hit[i].end());
            block_hit[i].erase(unique(block_hit[i].begin(), block_hit[i].end()), block_hit[i].end());
        }
        long sum_cl = 0;
        for (int i = 0; i < NBLOCKS; i ++) {
            sum_cl += block_hit[i].size();
        }
        printf("total cache line size is %ld\n", sum_cl);
        printf("block 0 cache line hit number is %d\n", block_hit[0].size());
        
    }
*/
    
    void reorder_graph(vector<int>& map) {
        // first reorder the graph
        int* tmp_degrees = new int[n];
        int* tmp_degrees_r = new int[n];
        long* tmp_nps = new long[n+1];
        long* tmp_nps_r = new long[n+1];
        int* tmp_edges = new int[m];

        // checking
        ////int checked_node = n/400;
        ////long start = nodepointer[0][checked_node];
        ////long end = nodepointer[0][checked_node+1];
        ////for (long i = start; i < end; i ++) {
        ////    printf("%d\n",map[edges[0][i]]);
        ////}
        // end of checking
        
        // map the degree arrays first
        for (int i = 0; i < n; i ++) {
            tmp_degrees[map[i]] = degree[0][i];
            tmp_degrees_r[map[i]] = in_degree[i];
        }
        // calculate the nodepointers
        tmp_nps[0] = tmp_nps_r[0] = 0;
        for (int i = 0; i < n; i ++) {
            tmp_nps[i+1]  = tmp_nps[i]  + tmp_degrees[i];
            tmp_nps_r[i+1]= tmp_nps_r[i]+ tmp_degrees_r[i];
        }
         
        nodepointer[0][n] = nodepointer[1][n] = m;
        // the edge map does have to relocate
        cilk_for (int i = 0; i < n; i ++) {
            // move old edges to new locations and map each node
            long start_o = nodepointer[0][i];   // old starting position
            long end_o = nodepointer[0][i+1];   // old starting position
            long start = tmp_nps[map[i]];       // new starting position
            long end   = tmp_nps[map[i]+1];       // new ending position
            if (start-end != start_o-end_o) {
                printf("wrong\n");
                abort();
            }
            for (long j = 0; j < end - start; j ++) {
                int node = edges[0][start_o+j]; // get the target of this edge
                tmp_edges[start+j] = map[node]; // map the node
            }
            merge_sort(tmp_edges, start, end-1);
        }

        // checking
        ////int mapped_node = map[checked_node];
        ////start = tmp_nps[mapped_node];
        ////end = tmp_nps[mapped_node+1];
        ////for (long i = start; i < end; i ++) {
        ////    printf("%d\n", tmp_edges[i]);
        ////}
        ////exit(0);
        // end of checking

        cilk_for(long i = 0; i < m; i ++) {
            edges[0][i] = tmp_edges[i];
        }
        // update degrees and nodepointers 
        for (int i = 0; i < n; i ++) {
            degree[0][i] = tmp_degrees[i];
            in_degree[i] = degree[1][i] = tmp_degrees_r[i];
            nodepointer[0][i] = tmp_nps[i];
            nodepointer[1][i] = tmp_nps_r[i];
        }
        // second reorder the random walk index
    }

    // For hot set,This function should be called after reorder the graph
    // This function should be called after index is loaded
    void reorder_random_walk_index(vector<int>& map) {
        long* tmp_dest_pointers = new long[n+1];

        cout << "------threshold scale inreorder_random_walk_index() is " << threshold_scale << endl;
        long long int result_idx = 0;
        for (int i = 0; i < n; i ++) {
            tmp_dest_pointers[i] = result_idx;
            //int node = permute_mapping_r[i];
            int node = i;
            long int degree_i = get_degree(node);
            long long int num_random_walk = (thres_type==0) ? ceil(rmax*degree_i*Omega) : degree_i;
            num_random_walk*= threshold_scale;
            if (num_random_walk != 0) {
                num_random_walk += 2;
            }
            result_idx += num_random_walk;
        }
        printf("result_idx = %lld\n", result_idx);
        tmp_dest_pointers[n] = result_idx;
        int* tmp_destinations = new int[result_idx];

        cilk_for(int i = 0; i < n; i ++) {
            long start_o = dest_pointers[i];
            long end_o = dest_pointers[i+1];
            long start = tmp_dest_pointers[map[i]];
            long end = tmp_dest_pointers[map[i]+1];
            if (start-end != start_o-end_o) { 
                printf("wrong\n");
                printf("%ld %ld\n", start-end, start_o-end_o);
                abort(); 
            }
            for (long j = 0; j < end-start; j ++) {
                int target = destinations[start_o+j];
                tmp_destinations[start+j] = map[target];
            }
        }

        cilk_for(long i = 0; i < result_idx; i ++) {
            destinations[i] = tmp_destinations[i];
        }
        cilk_for(int i = 0; i < n; i ++) {
            dest_pointers[i] = tmp_dest_pointers[i];
        }
        cout << "======random walk index reordered!" << endl;
    }
    // This function should be called after separation
    void hot_set_schedule(int ncore) {
        load_hot_set();
        // reorder the graph (and random walk index) based on hot set statistics
        reorder_graph(hot_set_map);
        //hot_set_length = n/hybridsize;
        hot_set_length = n;

        hs_block_size = hot_set_length/(ncore*8*GRAINSIZE_SCALE);
        grain_size = hs_block_size;
        NTHREADS = ncore;
        NROUNDS  = hot_set_length / (hs_block_size * ncore) + 1;
        NBLOCKS  = hot_set_length / hs_block_size + 1;
        NCACHELINE = n/8+1;

        printf("NBLOCKS = %d\n", NBLOCKS);

        init_hot_set_block_cacheline(ncore);
        init_hot_set_block_matrix(ncore);

        #ifdef GRAINSIZE_SCALE
        string file_str = data_folder + "/grain_hot_set_block_mixed_"+to_string(8*GRAINSIZE_SCALE)+"p.txt";
        #else
        printf("grain size scale not set, please set first\n");
        exit(0);
        #endif
        // timer
        struct timespec start, finish;  clock_gettime(CLOCK_REALTIME, &start);
        double elapsed=0.0;

        if(exists_test(file_str)){
            printf("load grain hot set block...%s\n",file_str.c_str());
            ifstream file(file_str);
            int id;
            int line = 0;
            int* info = new int[NBLOCKS];
            while (file >> id) {
                info[line++] = id;
            }
            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    grain_blocks[i][j] = info[i*NTHREADS+j];
                }
            }
            delete[] info;

        } else {
            printf("file not found...%s\n",file_str.c_str());
            //exit(0);
            #ifdef NOTHING
            for (int round = 0; round < NROUNDS; round ++) {
                int start = round * ncore;
                int end = start + ncore -1;
                for (int i = 0; i < ncore/2; i ++) {
                    grain_blocks[round][i] = start+i*2;
                    grain_blocks[round][ncore-1-i] = start+i*2+1;
                }
            }
            #else // NOTHING
            compute_block_hits();
            KMinSketch sketch = KMinSketch(128*2, NCACHELINE, NBLOCKS, block_hit);
            sketch.calculate_all_block_size();
            sketch.compute_all_conflict_scores();

            // firstly, find blocks with largest size
            vector<pair<int, double> > block_pair; // <block_id, block_size>;
            for (int i = 0; i < NBLOCKS; i ++) {
                block_pair.push_back(make_pair(i, NBLOCKS-i));
            }
            sort(block_pair.begin(),block_pair.end(), pair_sorter_large_first);
//            setWorkers(40);
            // minimize contention in first few rounds
            int thres1 = NTHREADS*(NROUNDS/10);
            int thres2 = NBLOCKS - thres1;
            int div1 = thres1/NTHREADS;
            #ifdef PHASE_TIMER
            phase1 = div1;
            #endif
            int div2 = thres2/NTHREADS;
            for (int round = 0; round < div1; round ++) {
                for (int i = 0; i < NTHREADS; i ++) {
                    if (i == 0) {
                        int first_blk = 0;
                        for (int idx = 0; idx < thres1; idx ++) {
                            int blk = block_pair[idx].first;
                            if (block_selected[blk] != 1) {
                                first_blk = blk;
                                break;
                            }
                        }
                        block_selected[first_blk] = 1;
                        grain_blocks[round][i] = first_blk;
                    } else {
                        int idx = 0;
                        double min_conflict = m*1.0;

                        // traverse all blocks
                        for (int j = 0; j < thres1; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double conflict = 0;
                                for (int prev = 0; prev < i; prev ++) {
                                    int blk_a = grain_blocks[round][prev];
                                    int blk_b = blk;
                                    double score = sketch.calculate_intersection(blk_a, blk_b);     
                                    //printf("score = %f\n", score); exit(0);
                                    conflict += score;
                                }
                                if (conflict < min_conflict) {
                                    idx = blk;
                                    min_conflict = conflict;
                                }
                                //printf("block %d checked\n", blk);
                            }
                        }

                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    }
                }
            }

            printf("phase 1 done\n");
            for ( int i = 0; i < NTHREADS; i ++){
                for (int round = div1; round < div2 ; round ++ ) {
                    if (round == 0) {
                        int idx = 0;
                        // find a block with smallest coverage
                        int min_cvg = m;
                        for (int j = thres1; j < thres2; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double block_size = sketch.calculate_block_size(blk);
                                if (block_size < min_cvg) {
                                    min_cvg = block_size;
                                    idx = blk;
                                }       
                            }
                        }
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        //double min_coverage= NCACHELINE;
                        cilk::reducer<MinMonoid<MinBlk>> min_coverage(MinBlk(NCACHELINE, 0));
                        //MinBlk mb(NCACHELINE,0);
                        //min_coverage = mb;
                        vector<int> existing_blks;
                        for (int prev = 0; prev < round; prev ++) {
                            int blk = grain_blocks[prev][i];
                            existing_blks.push_back(blk);
                        }
                        // traverse all blocks

                        cilk_for (int j = thres1; j < thres2; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double coverage = 0;
                                // compute the union of this blk and existing blks
                                coverage = sketch.compute_union(existing_blks, blk);
                                MinBlk mb(coverage, blk);
                                if ((*(min_coverage.get_value())).number > mb) {
                                    (min_coverage.get_value())->number = mb;
                                }
                            }
                        }

                        idx = min_coverage.get_value()->number.idx;
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    }
                    //printf("phase 2 thread %d round %d done\n", i, round);
                }
                printf("phase 2 thread %d done\n", i);
            }
            printf("phase 2 done\n");
            
            for ( int i = 0; i < NTHREADS; i ++){
                for (int round = div2; round < NROUNDS && round*NTHREADS+i < NBLOCKS ; round ++ ) {
                    if (round == 0) {
                        int idx = 0;
                        // find a block with smallest coverage
                        int min_cvg = m;
                        for (int j = thres2; j < NBLOCKS; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double block_size = sketch.calculate_block_size(blk);
                                if (block_size < min_cvg) {
                                    min_cvg = block_size;
                                    idx = blk;
                                }       
                            }
                        }
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        //double min_coverage= NCACHELINE;
                        cilk::reducer<MinMonoid<MinBlk>> min_coverage(MinBlk(NCACHELINE, 0));
                        //MinBlk mb(NCACHELINE,0);
                        //min_coverage = mb;
                        vector<int> existing_blks;
                        for (int prev = 0; prev < round; prev ++) {
                            int blk = grain_blocks[prev][i];
                            existing_blks.push_back(blk);
                        }
                        // traverse all blocks

                        cilk_for (int j = 0; j < NBLOCKS; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double coverage = 0;
                                // compute the union of this blk and existing blks
                                coverage = sketch.compute_union(existing_blks, blk);
                                MinBlk mb(coverage, blk);
                                if ((*(min_coverage.get_value())).number > mb) {
                                    (min_coverage.get_value())->number = mb;
                                }
                            }
                        }

                        idx = min_coverage.get_value()->number.idx;
                        block_selected[idx] = 1;
//                        printf("round = %d i=%d selected idx is %d\n", round,i, idx);
                        grain_blocks[round][i] = idx;
                    }
                }
            }
            #endif // ifdef NOTHING
            FILE *file = fopen(file_str.c_str(), "w");

            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    printf("%d ", grain_blocks[i][j]);
                    fprintf(file, "%d\n", grain_blocks[i][j]);
                }
                printf("\n");
            }
            fflush(file);
            fclose(file);
            clock_gettime(CLOCK_REALTIME, &finish);
            elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
            elapsed += (finish.tv_nsec - start.tv_nsec);
            printf("block matrixgenerated in %.3f seconds using %d workers.\n",
                   elapsed/1000000000, __cilkrts_get_nworkers());

            printf("%s\n", file_str.c_str());
//            setWorkers(1);
            exit(0);
        }
    }
    
    //--------------------- HOT SET SCHEDULING -----------------------------END
    
    //=====================NUMA AWARE REPLICATES=====================
    //Replicate partial residue array for nodes that is frequently updated
    vector<int> pop_map;    // measures the popularity of each node
    vector<int> pop_map_r;  // reverse of pop_map
    void load_popularity() {
        string pop_file_str = data_folder + "/stats.txt";
        pop_map.reserve(n);
        pop_map_r.reserve(n);

        for (int i = 0; i < n; i ++) { 
            pop_map.push_back(0);
            pop_map_r.push_back(0);
        }
        if (!exists_test(pop_file_str)) {
            printf("stats.txt not found\n");
            exit(0);
        } else {
            FILE *stats = fopen(pop_file_str.c_str(), "r");
            char line[200];
            char* firstline;
            fgets(line, sizeof line, stats);

            while (fgets(line, sizeof line, stats) != NULL) {
                int node, idx;
                double pagerank_i, popularity_i;
                long int din_i, dout_i;
                sscanf(line, "%d %d %lf %lf %ld %ld\n", &idx,  &node, &popularity_i, 
                                &pagerank_i, &dout_i, &din_i);
                pop_map[node] = idx;
                pop_map_r[idx] = node;
            }
            fclose(stats);
        }

        printf("stats loaded!\n");
        // check loaded stats
        //for (int i = 0; i < 20; i ++) {
            //printf("%d \n", pop_map[i]);
        //    printf("%d \n", pop_map_r[i]);
        //}
        //exit(0);
    }
    void map_rwidx_pop() {
        // map each target to its new position        
        for (long long int i = 0; i < dest_pointers[dividingline]; i ++) {
            int t = destinations[i];
            destinations[t] = pop_map[t];
        }
    }
    void map_rwidx_degree() {
        // map each target to its new position        
        for (long long int i = 0; i < dest_pointers[dividingline]; i ++) {
            int t = destinations[i];
            destinations[t] = degree_order[t];
        }
    }
    int * random_map=NULL;
    int * random_map_r=NULL;
    void generate_random_map() {
        random_map = new int[n];
        memset(random_map, 0, sizeof(int) * n);

        int* init_map = new int[n];
        #pragma simd
        for (int i = 0; i < n; i ++) {
            init_map[i] = i;
        }
        // shuffle the positions
        for (int i = 0; i < n; i ++) {
            int pos = rand() % n;
            int tmp = init_map[i];
            init_map[i] = pos;
            init_map[pos] = tmp;
        }
        // get the map
        for (int i = 0; i < n; i ++) {
            random_map[init_map[i]] = i;
        }
        random_map_r = init_map;
    }
    void map_rwidx_random() {
        // map each target to its new position        
        for (long long int i = 0; i < dest_pointers[dividingline]; i ++) {
            int t = destinations[i];
            destinations[t] = random_map[t];
        }
    }

    // target map
    int* target_map=NULL;
    int* target_map_r=NULL;
    int target_divideline = 0;
    void load_target_map() {
        
        target_map = new int[n];
        target_map_r = new int[n];
        memset(target_map, 0, sizeof(int) * n);
        memset(target_map_r,0,sizeof(int) * n);

        string rwcount_file_str = data_folder + "/rwidx3_counts.txt";
        bool isLineSet = false;
        int count_line = 0;

        //if (true) { // always generate new target maps
                    // because each time rw index is different
        if (!exists_test(rwcount_file_str)) {
            printf("rwidx3_counts.txt not found\n");
            check_index_targets();
            //exit(0);
        } 
            FILE *stats = fopen(rwcount_file_str.c_str(), "r");
            char line[200];
            char* firstline;
            fgets(line, sizeof line, stats);
        int maxcount = 0;
            while (fgets(line, sizeof line, stats) != NULL) {
                int node, idx;
                double pagerank_i, count_i;
                long int din_i, dout_i;
                sscanf(line, "%d %d %lf %ld %ld\n", &idx,  &node, &count_i, 
                                &dout_i, &din_i);
                #ifdef HOT_SET
                node = hot_set_map[node];
                #endif
                target_map[node] = idx;
                target_map_r[idx] = node;
                if ( count_i > 255) {
                    target_divideline ++;
                }
                count_line ++;

                if (maxcount < count_i) maxcount = count_i;
                //if ( count_line < 20) {
                //    printf("%d %d %f %ld %ld\n", idx,  node, count_i, 
                //                dout_i, din_i);
                //}
            }
            printf("target counts loaded from %s\n", rwcount_file_str.c_str());
            fclose(stats);

        printf("rwidx stats loaded!\n");
        printf("target dividing line = %d\n", target_divideline);
        printf("count line  = %d\n", count_line);
        // check loaded stats
        //for (int i = 0; i < 20; i ++) {
          //printf("%d \n", pop_map[i]);
        //    printf("%d \n", target_map_r[i]);
        //}
        //exit(0);
        if (target_divideline == 0 ) {
                
            return;
        }
            cout << "target counts max is " << maxcount << endl;
            //exit(0);
        // randomly permute target map (to reduce contention on popular targetrs)
        for (int i = 0; i < target_divideline; i ++) {
            int pos = rand() % target_divideline;
            // swap
            int tmp = target_map_r[pos];
            target_map_r[pos] = target_map_r[i];
            target_map_r[i] = tmp;
        }
        for (int i = 0; i < target_divideline; i ++) {
            target_map[target_map_r[i]] = i;
        }

        int length = n - target_divideline;
        for (int i = target_divideline; i < n; i ++) {
            int pos = rand() % length + target_divideline;
            
            int tmp = target_map_r[pos];
            target_map_r[pos] = target_map_r[i];
            target_map_r[i] = tmp;
        }
        for (int i = target_divideline; i < n; i ++) {
            target_map[target_map_r[i]] = i;
        }
    }
    void map_rwidx_target() {
        // map each target to its new position        
//        setWorkers(40);
        cout << "dest pointers end = " << dest_pointers[dividingline];
        cilk_for (long long int i = 0; i < dest_pointers[dividingline]; i ++) {
            int t = destinations[i];
            destinations[i] = target_map[t];
        }
        cout << "index mapped" << endl;
    }
        
    int rw_dvdngline=0; // separates nodes by whether they hit popular targets
                        // a popular target means hit more than 255 times
    int* pop_target_map=NULL;
    int* pop_target_map_r=NULL;
    // This function should be called after map_idx_target()
    void separate_random_walk_sources() {
        pop_target_map = new int[n];
        pop_target_map_r = new int[n];
        fill_n(pop_target_map, n, -1);
        fill_n(pop_target_map_r, n, -1);
        bool* hit = new bool[n];
        fill_n(hit, n, false);
//        setWorkers(40);
        cilk_for (int node = 0; node < n; node ++) {
            long long start = dest_pointers[node];
            long long end = dest_pointers[node+1];
            for (long long i = start; i < end; i ++) {
                int t = destinations[i];
                if (t < target_divideline) hit[node] = true;
            }
        }
        int idx = 0;
        for (int i = 0; i < n; i ++) {
            if (hit[i]) {
                pop_target_map[i] = idx++;
            } 
        }
        rw_dvdngline = idx;
        for (int i = 0; i < n; i ++) {
            if (!hit[i]) {
                pop_target_map[i] = idx++;
            }
        }
        cilk_for(int i = 0; i < n; i ++) {
            pop_target_map_r[pop_target_map[i]] = i;
        }
        printf("rw_dvdngline is %d\n", rw_dvdngline);
        exit(0);
    }

    // =======================================================================
    // optimizing for grain blocks
    int grain_size=0; //cilk for grainsize
    int NTHREADS = 40;
    // select cache contention less grain blocks
    int** grain_blocks=NULL;
    atomic<double>* indices=NULL;
    int* block_selected=NULL;
    int* nrounds=NULL;
    int NROUNDS;
    int NBLOCKS;
    int NCACHELINE=0;
    vector<vector<int>> block_hit; 
    

    void init_grain_size() {
        //grain_size = n / (8*NTHREADS);
        if (grain_size != 0) return;
        grain_size = dividingline / (8*NTHREADS);
        // round grain size to 2^x
        int idx = -1;
        printf("grain size = %d\n", grain_size);
        while (grain_size > 0) {
            idx++;
            grain_size /= 2;
        }
        grain_size = 1 << (idx);
        #ifdef GRAINSIZE_SCALE
        grain_size /= GRAINSIZE_SCALE;
        #endif
        grain_size = 1024;
        if (dividingline > 40000000) {
            grain_size *= 2;
        }
        if (dividingline > 100000000) {
            grain_size *= 2;
        }
        printf("rounded grain size = %d\n", grain_size);
        NROUNDS = dividingline / (grain_size*NTHREADS) + 1;
        NBLOCKS = dividingline / grain_size + 1;
        NCACHELINE = n/8+1;
    }

    // grain blocks is the block matrix
    void init_grain_blocks() {
        if (grain_blocks == NULL) {
            grain_blocks = new int*[NROUNDS];
            for (int i = 0; i < NROUNDS; i ++) {
                grain_blocks[i] = new int[NTHREADS];
                memset(grain_blocks[i], -1, sizeof(int)*NTHREADS);
            }
        }
    }
    void reset_block_masks() {
        for (int i = 0; i < NTHREADS; i ++) indices[i] = 0;
    }
    void insert_workloads(int blk_id, vector<int>& extg_nodes) {
        int start_offset = blk_id * grain_size;
        for (int i = 0; i < grain_size; i ++) {
            int node = i + start_offset;
            if (node >= n) break;
            extg_nodes.push_back(node);
            // insert out neighbors
            long start = nodepointer[0][node];
            long   end = nodepointer[0][node+1];
            for (long offset = start; offset < end; offset ++) {
                int neighbor = edges[0][offset];
                extg_nodes.push_back(neighbor);
            }
        }
    }

    inline int get_cache_line(int node) {
        return node/8;
    }

    long get_contention_score(int blk_id, int min_ctn_score, vector<int>& extg_nodes) {
        // iterate all nodes in this block
        vector<int> blk_nodes = vector<int>();
        insert_workloads(blk_id, blk_nodes);
        long contention_score = 0;
        for (int node : blk_nodes) {
            for (int extg : extg_nodes) {
                if (get_cache_line(node) == get_cache_line(extg)) {
                    contention_score ++;
                    if (contention_score >= min_ctn_score)
                        return INT_MAX;
                }
            }
        }
        return contention_score;
    }
    void generate_block_orders() {
        // set grain size
        init_grain_size();

        init_grain_blocks();
    
        for (int round = 0; round < NROUNDS; round ++) {
            
            vector<int> extg_nodes = vector<int>(); // all nodes in existing workload
            for (int i = 0; i < NTHREADS && round*NTHREADS+i < NBLOCKS; i ++) {
                // set first grain block
                if (i == 0) {
                    int idx = 0;
                    while (block_selected[idx++] != -1)
                       ;
                    grain_blocks[round][i] = idx;
                    block_selected[idx] = 1;
                    insert_workloads(idx, extg_nodes);
                } else {
                // set the remaining grain blocks
                    long min_contention_score = INT_MAX;
                    int result = -1;
                    int debug_count = 0;
                    // all existing workload are stored in extg_nodes
                    // traverse each block
                    for ( int blkidx = 0; blkidx < NBLOCKS; blkidx ++) {
                        if (block_selected[blkidx] == -1) {
                            // calculate contention score of this block
                            long contention_score = 0;
                            int block_id = blkidx;
                            contention_score = get_contention_score(block_id,
                                                min_contention_score,
                                                extg_nodes);
                            if (contention_score < min_contention_score) {
                                min_contention_score = contention_score;
                                result = block_id;
                            }
                            printf("%d block checked\n", ++debug_count);
                        }   
                    }
                    // set this block
                    grain_blocks[round][i] = result;
                    block_selected[result] = 1;
                    insert_workloads(result, extg_nodes);
                }
                printf("%d blocks set\n", round*NTHREADS+i);
                
                
            }
        }
        printf("block orders set\n");
        exit(0);
    }

    // init several important data structures
    //      block_hit:  the cache lines that a block's out neighbors reside
    //      indices:    the indices for each column in block matrix
    //                  used in forward push
    //      block_selected: mark blocks that has been arranged in the matrix
    //      nrounds:    the length of each column in block matrix
    void init_for_block_cacheline() {
        //cache_line= vector<vector<int>>(NCACHELINE, vector<int>());

        if (indices == NULL) {
            indices = new atomic<double>[NTHREADS];
        }
        for (int i = 0; i < NTHREADS; i ++) indices[i] = 0;
        
        if (block_selected == NULL) {
            block_selected = new int[NBLOCKS];
        }
        memset(block_selected, -1, sizeof(int)*(NBLOCKS));

        printf("NBLOCKS=%d NROUNDS=%d\n",  NBLOCKS, NROUNDS);
        if (nrounds == NULL) {
            nrounds = new int[NTHREADS];
            for (int i = 0; i < NTHREADS; i ++) {
                if (NTHREADS*(NROUNDS-1) + i < NBLOCKS) nrounds[i] = NROUNDS;
                else nrounds[i] = NROUNDS - 1;
                printf("%d ", nrounds[i]);
            }
        }
        printf("\n");
    }
    
    void compute_block_hits() {
        block_hit = vector<vector<int>>(NBLOCKS, vector<int>());
        // traverse the graph
        printf("grain size in compute_block_hits() is %d\n", grain_size);
//        setWorkers(80);
        #pragma cilk grainsize=1
        cilk_for (int blk = 0; blk < NBLOCKS; blk ++) {
            int offset = blk*grain_size;
            for (int i = offset; i<n && i < offset+grain_size; i ++) {
                long start = nodepointer[0][i];
                long end   = nodepointer[0][i+1];
                int block_id = i / grain_size;
                // check each out neighbor
                vector<int> candidates = vector<int>(); // duplications exist in this vector
                candidates.push_back(-1); // to avoid boundary check
                
                for (long j = start; j < end; j ++) {
                    int node = edges[0][j];
                    int cacheline_id = node/(CACHELINE_SIZE/sizeof(double));
                    candidates.push_back(cacheline_id);
                }

                for (int j = 1; j < candidates.size(); j ++) {
                    int cacheline_id = candidates[j];
                    // check duplicates in each node's out neighbors
                    if (cacheline_id != candidates[j-1]) {
                        block_hit[block_id].push_back(cacheline_id);
                    }
                }
            }
            
        }
        // remove duplicates in each block
        #pragma cilk grainsize=1
        cilk_for (int i = 0; i < NBLOCKS; i ++) {
            sort(block_hit[i].begin(), block_hit[i].end());
            block_hit[i].erase(unique(block_hit[i].begin(), block_hit[i].end()), block_hit[i].end());
        }
        long sum_cl = 0;
        for (int i = 0; i < NBLOCKS; i ++) {
            sum_cl += block_hit[i].size();
        }
        printf("total cache line size is %ld\n", sum_cl);
        printf("block 0 cache line hit number is %d\n", block_hit[0].size());
    }

    
    struct KMinSketch {
        int nrows;
        long long ncolumns;
        double* sketches=NULL;
        double* conflict_score=NULL;
        long nblocks;
        double* block_size=NULL;
        double* kmins=NULL;
        double* unions=NULL;
        vector<vector<int>>& block_hits;

        KMinSketch(int r, int c, int nblks, 
                    vector<vector<int>>& blk_hits):block_hits(blk_hits) {
            nrows = r;
            ncolumns = c;
            nblocks = nblks;
            //block_hits = blk_hits;
            cout <<"generating random numbers" << endl;
//            setWorkers(80);
            sketches = new double[nrows*ncolumns];
            cilk_for (long long i = 0; i < nrows*ncolumns; i ++) {
                double rdm = (lrandp()%INT_MAX) * 1.0 / ((double)INT_MAX);
                sketches[i] = rdm;
            }
            conflict_score = new double[nblocks*nblocks];
            fill_n(conflict_score, nblocks*nblocks, -1);
            block_size = new double[nblocks];
            fill_n(block_size, nblocks, -1);
            
            // init k mins for each block
            kmins = new double[nrows*nblocks];
            fill_n(kmins, nrows*nblocks, 1);
            
            cout <<"generating kmins" << endl;
            for (int i = 0; i < nrows; i ++) {
                cilk_for(int blk = 0; blk < nblocks; blk ++) {
                    // get the min value for this block in this row
                    double min = 1;
                    for (int cl_id : block_hits[blk]) {
                        if (min > sketches[i*ncolumns+cl_id])
                            min = sketches[i*ncolumns+cl_id];
                    }
                    kmins[blk*nrows+i] = min;
                }
            }
            unions = new double[nblocks*nblocks];
            fill_n(unions, nblocks*nblocks, -1);
//            setWorkers(1);
            printf("sketch generated\n");
        }
        ~KMinSketch() {
            if (sketches != NULL)        delete[] sketches;
            if (conflict_score  != NULL) delete[] conflict_score;
            if (block_size      != NULL) delete[] block_size;
            if (kmins    != NULL)        delete[] kmins;
        }

        double calculate_block_size(int idx) {
            if (block_size[idx] >= 0) return block_size[idx];

            double sum = 0;
            for (int i = 0; i < nrows; i ++) {
                sum += kmins[idx*nrows+i];
            }
            double result = (1/(sum/nrows) -1);
            block_size[idx] = result;

            //if (result == 0) {
            //   printf("");
            //}
            return result;
        }

        double calculate_union_size(int a, int b) {
            if (unions[a*nblocks+b] >= 0 ) return unions[a*nblocks+b];
            
            double sum = 0;
            for (int i = 0; i < nrows; i ++) {
                double min_a = kmins[a*nrows+i];
                double min_b = kmins[b*nrows+i];
                sum += min_a < min_b ? min_a : min_b;
            }
            double result = (1/(sum/nrows) -1);
            unions[a*nblocks+b] = result;
            unions[b*nblocks+a] = result;
            return result;
        }
        double calculate_intersection(int a, int b) {
            if (conflict_score[a*nblocks+b] >= 0) 
                return get_intersection(a, b);

            double size_a = calculate_block_size(a);
            double size_b = calculate_block_size(b);
            // calculate the union size
            double size_union = calculate_union_size(a, b);
            
            double result = size_a + size_b - size_union;
            if (result < 0) result=0;;
            conflict_score[a*nblocks+b] = result;
            return result;
        }

        double get_intersection(int blk_a, int blk_b) {
            return conflict_score[blk_a*nblocks+blk_b];
        }

        void calculate_all_block_size() {
//            setWorkers(40);
            cilk_for(int i = 0; i < nblocks; i ++) {
                calculate_block_size(i);
            }
           
//            setWorkers(1);
        }
        double get_block_size(int i) {
            return calculate_block_size(i);
        }
        void compute_all_conflict_scores() {
//            setWorkers(80);
            cilk_for (int i = 0; i < nblocks; i ++) {
                calculate_block_size(i);
            }

            for (int i = 0; i < nblocks; i ++) {
                cilk_for (int j = 0; j < nblocks; j ++) {
                    calculate_intersection(i,j);
                }
            }
            //for (int j = 0; j < NBLOCKS/8; j ++) {
            //    printf("%f ", conflict_score[0][j]);
            //}
            //exit(0);
//            setWorkers(1);
        }

        double compute_union(vector<int>& blks, int idx) {
            double sum = 0;
            int limit = 15;
            int blk_size = blks.size();
            for (int i = 0; i < nrows; i ++) {
                double min = 1;
                if (blks.size() <= limit) {
                    for (int blk: blks) {
                        double min_blk = kmins[blk*nrows+i];
                        if (min_blk < min) min = min_blk;
                    }
                } else {
                    for (int id = blk_size-limit; id < blk_size; id ++) {
                        int blk = blks[id];
                        double min_blk = kmins[blk*nrows+i];
                        if (min_blk < min) min = min_blk;
                    }
                }
                double min_idx = kmins[idx*nrows+i];
                if (min_idx < min) min = min_idx;
                sum += min;
            }
            double result = (1/(sum/nrows) -1);
            return result;
        }
    };
    
    void generate_block_order_by_cache_list() {
        init_grain_size();
        init_for_block_cacheline();
        init_grain_blocks();

        //#ifdef GRAINSIZE_SCALE
        //string file_str = data_folder + "/grain_block_matrix40_"+to_string(8*GRAINSIZE_SCALE)+"p.txt";
        //#else
        string file_str = data_folder + "/grain_block_matrix40.txt";
        //#endif

        struct timespec start, finish;
        double elapsed=0.0;
        if(exists_test(file_str)){
            printf("load grain block matrix\n");
            cout << file_str << endl;
            ifstream file(file_str);
            int id;
            int line = 0;
            int* info = new int[NBLOCKS];
            while (file >> id) {
                info[line++] = id;
            }
            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    grain_blocks[i][j] = info[i*NTHREADS+j];
                }
            }
            delete[] info;

        } else {
            compute_block_hits();
            KMinSketch clsketch = KMinSketch(64, NCACHELINE, NBLOCKS, block_hit);
            clsketch.compute_all_conflict_scores();
            for (int round = 0; round < NROUNDS; round ++ ) {
                for ( int i = 0; i < NTHREADS && round*NTHREADS+i < NBLOCKS; i ++){
                    if (i == 0) {
                        int idx = 0;
                        while (block_selected[idx] != -1) {idx++;}
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        double min_conflict = m*1.0;
                        // traverse all blocks
                        for (int blk = 0; blk < NBLOCKS; blk ++) {
                            if (block_selected[blk] != 1) {
                                double conflict = 0;
                                for (int prev = 0; prev < i; prev ++) {
                                    int blk_a = grain_blocks[round][prev];
                                    int blk_b = blk;
                                    double score = clsketch.calculate_intersection(blk_a, blk_b);     
                                    //printf("score = %f\n", score); exit(0);
                                    conflict += score;
                                }
                                if (conflict < min_conflict) {
                                    idx = blk;
                                    min_conflict = conflict;
                                }
                                //printf("block %d checked\n", blk);
                            }
                        }

                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    }
                }
                //printf("%d rounds finished\n", round);
            }
            printf("block order finished\n");

            permute_blocks_for_locality(clsketch);

            FILE *file = fopen(file_str.c_str(), "w");

            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    //printf("%d ", grain_blocks[i][j]);
                    fprintf(file, "%d\n", grain_blocks[i][j]);
                }
                //printf("\n");
            }
            fflush(file);
            fclose(file);
            exit(0);
        }
    }

    void permute_blocks_for_locality(KMinSketch& sketch) {
        // cache conflict score stores in the sketch
        // process each round seperately
        for (int rnd = 1; rnd < NROUNDS; rnd ++) {
            bool* is_taken = new bool[NTHREADS];
            for (int i = 0; i < NTHREADS; i ++) is_taken[i] = false;
            int* new_blocks = new int[NTHREADS];
            for (int i = 0; i < NTHREADS; i ++) new_blocks[i] = -1;
            // find best option for each thread
            for (int thd = 0; thd < NTHREADS; thd ++) {
                int blk_a = grain_blocks[rnd-1][thd];
                // get its block in last round
                double max_conflict = 0;
                int max_idx = -1;
                // check each block in this row
                for (int i = 0; i < NTHREADS; i ++) {
                    if (!is_taken[i]) {
                        int blk_b = grain_blocks[rnd][i];
                        double conflict = sketch.get_intersection(blk_a, blk_b);
                        if (conflict > max_conflict) {
                            max_conflict = conflict;
                            max_idx = i;
                        }
                    }
                }
                // we now have found a block that has the most conlicts
                // which also means best memory locality
                new_blocks[thd] = grain_blocks[rnd][max_idx];
                is_taken[max_idx] = true;
            }
            for (int thd = 0; thd < NTHREADS; thd ++)
                grain_blocks[rnd][thd] = new_blocks[thd];
            //delete[] is_taken;
            //delete[] new_blocks;
        }
    }
    void generate_block_order_locality_first() {
        init_grain_size();
        init_for_block_cacheline();
        init_grain_blocks();
        printf("NCACHELINE is %d\n", NCACHELINE);

        #ifdef GRAINSIZE_SCALE
        string file_str = data_folder + "/grain_block_matrix40_locality_first_"+to_string(8*GRAINSIZE_SCALE)+"p.txt";
        #else
        printf("grain size scale not set, please set first\n");
        exit(0);
        #endif
        struct timespec start, finish;
        double elapsed=0.0;
        clock_gettime(CLOCK_REALTIME, &start);
        if(exists_test(file_str)){
            printf("load grain block matrix...\n");
            ifstream file(file_str);
            int id;
            int line = 0;
            int* info = new int[NBLOCKS];
            while (file >> id) {
                info[line++] = id;
            }
            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    grain_blocks[i][j] = info[i*NTHREADS+j];
                }
            }
            delete[] info;

        } else {
            compute_block_hits();
            KMinSketch clsketch = KMinSketch(64, NCACHELINE, NBLOCKS, block_hit);
            clsketch.calculate_all_block_size();

            for (int i = 0; i < NBLOCKS; i ++) {
                if (i % (NBLOCKS/16) == 0) {
                    int size = block_hit[i].size();
                    double est = clsketch.calculate_block_size(i);
                    printf("block %d coverage=%ld, estimation=%.0f, error= %.0f%\n",
                            i, size, est, (est-size)*100/size);
                }
            }

            double min_size = m;
            double max_size = 0;
            for (int i = 0; i < NBLOCKS; i ++) {
                double size = clsketch.calculate_block_size(i);
                if (size > max_size) max_size = size;
                if (size < min_size) min_size = size;
                //printf("%.0f ", clsketch.calculate_block_size(i));
            }
            //printf("\n");
            printf("min size %f, max size %f\n", min_size, max_size);
            //exit(0);
//            setWorkers(40);
            for ( int i = 0; i < NTHREADS; i ++){
                for (int round = 0; round < NROUNDS && round*NTHREADS+i < NBLOCKS; round ++ ) {
                    if (round == 0) {
                        int idx = 0;
                        // find a block with smallest coverage
                        int min_cvg = m;
                        for (int blk = 0; blk < NBLOCKS; blk ++) {
                            if (block_selected[blk] == -1) {
                                double block_size = clsketch.calculate_block_size(blk);
                                if (block_size < min_cvg) {
                                    min_cvg = block_size;
                                    idx = blk;
                                }       
                            }
                        }
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        //double min_coverage= NCACHELINE;
                        cilk::reducer<MinMonoid<MinBlk>> min_coverage(MinBlk(NCACHELINE, 0));
                        //MinBlk mb(NCACHELINE,0);
                        //min_coverage = mb;
                        vector<int> existing_blks;
                        for (int prev = 0; prev < round; prev ++) {
                            int blk = grain_blocks[prev][i];
                            existing_blks.push_back(blk);
                        }
                        // traverse all blocks

                        cilk_for (int blk = 0; blk < NBLOCKS; blk ++) {
                            if (block_selected[blk] != 1) {
                                double coverage = 0;
                                // compute the union of this blk and existing blks
                                coverage = clsketch.compute_union(existing_blks, blk);
                                MinBlk mb(coverage, blk);
                                if ((*(min_coverage.get_value())).number > mb) {
                                    (min_coverage.get_value())->number = mb;
                                }
                            }
                        }

                        idx = min_coverage.get_value()->number.idx;
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                        //printf("selected idx is %d\n", idx);
                        if (round == nrounds[i]-1) {
                            printf("idx = %d total coverage is %f\n", idx, clsketch.compute_union(existing_blks, idx));
                            printf("intersection is %f\n",clsketch.calculate_intersection(2750, 2763));
                            printf("union is %f\n",clsketch.calculate_union_size(2750, 2763));
                            printf("set size is %f\n",clsketch.calculate_block_size(2763));
                            printf("set size is %f\n",clsketch.calculate_block_size(2750));
                        }
                    }
                }
                printf("%d threads finished\n", i);
            }
            printf("block order finished\n");
//            setWorkers(1);
            FILE *file = fopen(file_str.c_str(), "w");

            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    //printf("%d ", grain_blocks[i][j]);
                    fprintf(file, "%d\n", grain_blocks[i][j]);
                }
                //printf("\n");
            }
            fflush(file);
            fclose(file);
            exit(0);
        }
    }

    int* get_block_matrix() {
        //#ifdef GRAINSIZE_SCALE
        //string file_str = data_folder + "/grain_block_matrix40_mixed_"+to_string(8*GRAINSIZE_SCALE)+"p.txt";
        //#else
        string file_str = data_folder + "/grain_block_matrix40_mixed.txt";
        //#endif
        #ifdef GORDER
        file_str += "gorder";
        #endif
        if(exists_test(file_str)){
            printf("load grain block matrix...%s\n",file_str.c_str());
            ifstream file(file_str);
            int id;
            int line = 0;
            int* info = new int[NBLOCKS];
            while (file >> id) {
                info[line++] = id;
            }
            return info;
        } else {
            printf("block matrix file does not exists!\n %s\n", file_str.c_str());
            permute_blocks_by_contention_locality_pop();
            exit(0);
        }
        return NULL;
    }
    void permute_blocks_by_contention_locality_pop() {
        init_grain_size();
        init_for_block_cacheline();
        init_grain_blocks();
        //{// this chunck is used to check the conflicts of large blocks
        //    compute_block_hits();
        //    KMinSketch sketch = KMinSketch(128*2, NCACHELINE, NBLOCKS, block_hit);
        //    sketch.calculate_all_block_size();
        //    sketch.compute_all_conflict_scores();
        //    //for (int i = 0; i < NTHREADS; i ++) {
        //    //    for (int j = 0; j < NTHREADS; j ++) {
        //    //        printf("%0.f ",sketch.get_block_size(i*NTHREADS+j));
        //    //    }
        //    //    printf("\n");
        //    //}


        //    // firstly, find blocks with largest size
        //    vector<pair<int, double> > block_pair; // <block_id, block_size>;
        //    for (int i = 0; i < NBLOCKS; i ++) {
        //        block_pair.push_back(make_pair(i, sketch.get_block_size(i)));
        //    }
        //    sort(block_pair.begin(),block_pair.end(), pair_sorter_large_first);

        //    // check conflicts
        //    for (int i = 0; i < NTHREADS; i ++) {
        //        int blk = block_pair[i].first;
        //        printf("blk[%d].size=%f \n",blk, sketch.get_block_size(blk));
        //        for (int j = 0; j < NTHREADS; j ++) {
        //            int blk_b = block_pair[j].first;
        //            double size = sketch.block_size[blk_b];
        //            double conflict = sketch.get_intersection(blk, blk_b);
        //            printf("(%d,%0.f,%0.f) ",blk_b,size, conflict);
        //        }
        //        printf("\n");
        //    }
        //    exit(0);
        //}
        printf("NCACHELINE is %d\n", NCACHELINE);
        //#ifdef GRAINSIZE_SCALE
        //string file_str = data_folder + "/grain_block_matrix40_mixed_"+to_string(8*GRAINSIZE_SCALE)+"p.txt";
        //#else
        string file_str = data_folder + "/grain_block_matrix40_mixed.txt";
        //#endif
        #ifdef GORDER
        file_str += "gorder";
        #endif
        struct timespec start, finish;
        clock_gettime(CLOCK_REALTIME, &start);
        double elapsed=0.0;
        if(exists_test(file_str)){
            printf("load grain block matrix...%s\n",file_str.c_str());
            ifstream file(file_str);
            int id;
            int line = 0;
            int* info = new int[NBLOCKS];
            while (file >> id) {
                info[line++] = id;
            }
            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    grain_blocks[i][j] = info[i*NTHREADS+j];
                }
            }
            delete[] info;

        } else {
            compute_block_hits();
            int k = 32;
            //k = 256;
            KMinSketch sketch = KMinSketch(k, NCACHELINE, NBLOCKS, block_hit);
            sketch.calculate_all_block_size();
            sketch.compute_all_conflict_scores();

            // firstly, find blocks with largest size
            vector<pair<int, double> > block_pair; // <block_id, block_size>;
            for (int i = 0; i < NBLOCKS; i ++) {
                block_pair.push_back(make_pair(i, sketch.get_block_size(i)));
            }
            sort(block_pair.begin(),block_pair.end(), pair_sorter_large_first);

            //for (int i = 0 ; i < NBLOCKS; i ++) {
            //    printf("%d %0.f\n", block_pair[i].first, block_pair[i].second);
            //}
            //exit(0);
//            setWorkers(40);
            // minimize contention in first few rounds
            int thres1 = NTHREADS*(NROUNDS/10);
            int thres2 = NBLOCKS - thres1;
            int div1 = thres1/NTHREADS;
            #ifdef PHASE_TIMER
            phase1 = div1;
            #endif
            int div2 = thres2/NTHREADS;
            for (int round = 0; round < div1; round ++) {
                for (int i = 0; i < NTHREADS; i ++) {
                    if (i == 0) {
                        int first_blk = 0;
                        for (int idx = 0; idx < thres1; idx ++) {
                            int blk = block_pair[idx].first;
                            if (block_selected[blk] != 1) {
                                first_blk = blk;
                                break;
                            }
                        }
                        block_selected[first_blk] = 1;
                        grain_blocks[round][i] = first_blk;
                    } else {
                        int idx = 0;
                        double min_conflict = m*1.0;

                        // traverse all blocks
                        for (int j = 0; j < thres1; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double conflict = 0;
                                for (int prev = 0; prev < i; prev ++) {
                                    int blk_a = grain_blocks[round][prev];
                                    int blk_b = blk;
                                    double score = sketch.calculate_intersection(blk_a, blk_b);     
                                    //printf("score = %f\n", score); exit(0);
                                    conflict += score;
                                }
                                if (conflict < min_conflict) {
                                    idx = blk;
                                    min_conflict = conflict;
                                }
                                //printf("block %d checked\n", blk);
                            }
                        }

                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    }
                }
            }

            printf("phase 1 done\n");
            for ( int i = 0; i < NTHREADS; i ++){
                for (int round = div1; round < div2 ; round ++ ) {
                    if (round == 0) {
                        int idx = 0;
                        // find a block with smallest coverage
                        int min_cvg = m;
                        for (int j = thres1; j < thres2; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double block_size = sketch.calculate_block_size(blk);
                                if (block_size < min_cvg) {
                                    min_cvg = block_size;
                                    idx = blk;
                                }       
                            }
                        }
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        //double min_coverage= NCACHELINE;
                        cilk::reducer<MinMonoid<MinBlk>> min_coverage(MinBlk(NCACHELINE, 0));
                        //MinBlk mb(NCACHELINE,0);
                        //min_coverage = mb;
                        vector<int> existing_blks;
                        for (int prev = 0; prev < round; prev ++) {
                            int blk = grain_blocks[prev][i];
                            existing_blks.push_back(blk);
                        }
                        // traverse all blocks

                        cilk_for (int j = thres1; j < thres2; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double coverage = 0;
                                // compute the union of this blk and existing blks
                                coverage = sketch.compute_union(existing_blks, blk);
                                MinBlk mb(coverage, blk);
                                if ((*(min_coverage.get_value())).number > mb) {
                                    (min_coverage.get_value())->number = mb;
                                }
                            }
                        }

                        idx = min_coverage.get_value()->number.idx;
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    }
                    //printf("phase 2 thread %d round %d done\n", i, round);
                }
                printf("phase 2 thread %d done\n", i);
            }
            printf("phase 2 done\n");
            
            for ( int i = 0; i < NTHREADS; i ++){
                for (int round = div2; round < NROUNDS && round*NTHREADS+i < NBLOCKS ; round ++ ) {
                    if (round == 0) {
                        int idx = 0;
                        // find a block with smallest coverage
                        int min_cvg = m;
                        for (int j = thres2; j < NBLOCKS; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double block_size = sketch.calculate_block_size(blk);
                                if (block_size < min_cvg) {
                                    min_cvg = block_size;
                                    idx = blk;
                                }       
                            }
                        }
                        block_selected[idx] = 1;
                        grain_blocks[round][i] = idx;
                    } else {
                        int idx = 0;
                        //double min_coverage= NCACHELINE;
                        cilk::reducer<MinMonoid<MinBlk>> min_coverage(MinBlk(NCACHELINE, 0));
                        //MinBlk mb(NCACHELINE,0);
                        //min_coverage = mb;
                        vector<int> existing_blks;
                        for (int prev = 0; prev < round; prev ++) {
                            int blk = grain_blocks[prev][i];
                            existing_blks.push_back(blk);
                        }
                        // traverse all blocks

                        cilk_for (int j = 0; j < NBLOCKS; j ++) {
                            int blk = block_pair[j].first;
                            if (block_selected[blk] != 1) {
                                double coverage = 0;
                                // compute the union of this blk and existing blks
                                coverage = sketch.compute_union(existing_blks, blk);
                                MinBlk mb(coverage, blk);
                                if ((*(min_coverage.get_value())).number > mb) {
                                    (min_coverage.get_value())->number = mb;
                                }
                            }
                        }

                        idx = min_coverage.get_value()->number.idx;
                        block_selected[idx] = 1;
//                        printf("round = %d i=%d selected idx is %d\n", round,i, idx);
                        grain_blocks[round][i] = idx;
                    }
                }
            }
            FILE *file = fopen(file_str.c_str(), "w");

            for (int i = 0; i < NROUNDS; i ++) {
                for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                    //printf("%d ", grain_blocks[i][j]);
                    fprintf(file, "%d\n", grain_blocks[i][j]);
                }
                //printf("\n");
            }
            fflush(file);
            fclose(file);
            clock_gettime(CLOCK_REALTIME, &finish);
            elapsed = (finish.tv_sec - start.tv_sec)*1000000000;
            elapsed += (finish.tv_nsec - start.tv_nsec);
            printf("k = %d block matrix generated in %.3f seconds using %d workers.\n", k, 
                   elapsed/1000000000, __cilkrts_get_nworkers());
            printf("%s\n", file_str.c_str());
//            setWorkers(1);
            exit(0);
        }
    }
    // bind VMA to numa nodes by block matrix
    long * sched_offsets=NULL;
    int * sched_edges=NULL;
    atomic<double> * residues=NULL;
    atomic<double> * reserves=NULL;
    double* sums= NULL;


    void numa_aware_load() {
        init_grain_size();
        init_for_block_cacheline();
        init_grain_blocks();

        int* block_handled = new int[NBLOCKS];
        memset(block_handled, 0, sizeof(int) * NBLOCKS);

        #ifndef HOT_SET
        int* block_matrix = get_block_matrix();
        #else
        int* block_matrix = new int[NBLOCKS];
        for (int i = 0; i < NROUNDS; i ++) {
            for (int j = 0; j < NTHREADS && i*NTHREADS+j < NBLOCKS; j ++) {
                block_matrix[i*NTHREADS+j] = grain_blocks[i][j] ;
            }
        }

        #endif
        #ifdef CACHE_SCHEDULE
        //numa_set_localalloc();
        //numa_set_interleave_mask(numa_all_nodes_ptr);
        #endif
        sched_offsets = new long[n+1];
        sched_edges = new int[m];

        //numa_set_interleave_mask(numa_all_nodes_ptr);
        residues = new atomic<double>[n];
        reserves = new atomic<double>[n];
        //cilk_for(int i = 0; i < n; i ++) {
        //    residues[0] = 0;
        //    reserves[0] = 0;
        //}
        //numa_set_localalloc();

        sums = new double[NTHREADS];
        memset(sums, 0, sizeof(double)*NTHREADS);
        
        int ncores = 40;
        setWorkers(ncores);
        atomic<int>* marking = new atomic<int>[ncores];
        memset(marking, 0, sizeof(int)*ncores);
#if defined(TWITTER) || defined(FRIENDSTER) || true //&& defined(NOTHING)
        cout << "thd cilk for begin" << endl;
        #pragma cilk grainsize = 1
        cilk_for (int thd = 0; thd < NTHREADS; thd ++) {
            struct timespec start, finish;
            
            cpu_set_t  mask;
            CPU_ZERO(&mask);
            CPU_SET(thd, &mask);
            sched_setaffinity(0, sizeof(mask), &mask);
            clock_gettime(CLOCK_REALTIME, &start);
            int core_id = sched_getcpu()%ncores;
            int check = Atomic_add(marking[core_id], 1);
            if (check != 0) {
                for (int i = 0; i < ncores; i ++) {
                    check = Atomic_add(marking[i], 1);
                    if (check == 0) {
                        core_id = i;
                        break;
                    }
                }
            }
            printf("%d ",core_id);
            
            for (int round = 0; round < nrounds[core_id]; round ++) {
                int block_id = block_matrix[round*NTHREADS+core_id];
                int start = block_id * grain_size; 
                int end = start + grain_size;
                if (end > n) end = n;
                for (int node = start; node < end; node ++) {
                    residues[node] = 0;
                    reserves[node] = 0;
                    sched_offsets[node] = nodepointer[0][node];
                    long long degree = get_degree(node);
                    long long start = sched_offsets[node];
                    long long end   = start + degree;
                    for (long long i = start; i < end; i ++) {
                        sched_edges[i] = edges[0][i];
                    }
                }
                block_handled[block_id] = 1;
            }
            CPU_ZERO(&mask);
            for (int i = 0; i < sysconf(_SC_NPROCESSORS_ONLN); i ++)
                CPU_SET(i, &mask);
            sched_setaffinity(0, sizeof(mask), &mask);
        } 
        cilk_for (int i = dividingline; i < n; i ++) {
            int node = i;
            residues[node] = 0;
            reserves[node] = 0;
            sched_offsets[node] = nodepointer[0][node];
            long long degree = get_degree(node);
            long long start = sched_offsets[node];
            long long end   = start + degree;
            for (long long j = start; j < end; j ++) {
                sched_edges[j] = edges[0][j];
            }
        }
        sched_offsets[n] = nodepointer[0][n];
        
        //cilk_for(long long i = 0; i < m; i ++) {
        //    sched_edges[i] = edges[0][i];
        //}
        //cilk_for(int i =  0; i <= n; i ++) {
        //    sched_offsets[i] = nodepointer[0][i];
        //}
        //sched_offsets[n] = nodepointer[0][n];
        //
#else   //
        ////int nnuma_node = numa_num_configured_nodes();
        ////printf("+++++++++++++++we have %d numa nodes\n", nnuma_node);
        ////exit(0);
        //cilk_for(int thd = 0; thd < NTHREADS; thd ++) {
        //    int core_id = sched_getcpu();
        //    int numa_node = numa_node_of_cpu(core_id);
        //    for (int i = numa_node; i < NBLOCKS; i += nnuma_node) {
        //        for (int j = 0; j < grain_size; j ++) {
        //            int node = i*NBLOCKS + j;
        //            if (node < n) {
        //                residues[node] = 0;
        //                reserves[node] = 0;
        //            }
        //        }
        //    }
        //}
        
        
        cout << "grain size is " << grain_size << endl;;
        cout << "rmax scaled is " << rmax_scaled << endl;;
        for (int i = 0; i < NBLOCKS; i ++) {
//            printf("%d ", block_matrix[i]);
        }
        
        bool* marked = new bool[n];
        memset(marked, 0, sizeof(bool)*n);
        #pragma cilk grainsize=1
        cilk_for(int blk = 0; blk < NBLOCKS; blk ++) {
            int core_id = sched_getcpu()%40;
            int block_id = -1;
            double idx = Atomic_add(indices[core_id], 1);
            int idx_i = (int)lround(idx);
            if (idx_i < nrounds[core_id]) {
                block_id = block_matrix[idx_i*NTHREADS+core_id];
            } else {
            // then try to steal others block
                for (int tid = 0; tid < NTHREADS; tid ++) {
                    double idx = Atomic_add(indices[tid], 1);
                    int idx_i = (int)lround(idx);
                    if (idx_i < nrounds[tid]) {
                        block_id = block_matrix[idx_i*NTHREADS+tid];
                        break;
                    }
                }
            }

            marked[block_id] = true;
            if (block_id >= 0 && block_id < NBLOCKS) {
                int start = block_id *grain_size; 
                int end = start + grain_size;

                for (int i = start; i < end && i < n; i ++) {
                    sched_offsets[i] = nodepointer[0][i];
                    long long degree = get_degree(i);
                    if (degree==0) {

                        residues[i] = 1;
                        reserves[i] = 1;
                    } else {
                        residues[i] = 1/degree;
                        reserves[i] = 1/degree;
                    }
                    if (i == 0) {
                        cout << "node 0 processed " << degree<<endl;
                    }
                    if (residues[i] >  1/10000){// && !graph.istrouble[i]) {
                        int node = i;
                        double residue = Atomic_load_reset(residues[node]);
                        // update reserve of this node
                        double reserve = residue * alpha;
                        atomic_add(reserves[node], reserve);

                        int num_out= nodepointer[0][node+1] - 
                                     nodepointer[0][node];
                        double avg_push_residue = (residue - reserve)/num_out;
                        
                        long long  start = nodepointer[0][node];
                        long long end = start + num_out;
                        //#pragma simd
                        for (long j = start; j < end; j++) {
                            // update the residues
                            int neighbor = edges[0][j];
                            //sched_edges[j] = edges[0][j];
                            //atomic_add(residues[neighbor], avg_push_residue);
                            if (j == 0) cout << "edges[0][0] set to " <<sched_edges[j] <<" " << edges[0][j]<< "by node" <<node<<endl; 
                        }   
                    }
                }
            } else {
                printf("block unhandled !!!!!!!!!!!!!\n");
            }
        }
        cilk_for(int i = dividingline; i <= n; i ++) {
            sched_offsets[i] = m;
            if (i < n) {
                residues[i] = 1;
                reserves[i] = 1;
            }
        }
        cilk_for(long long i = 0; i < m; i++) {
            //sched_edges[i] = edges[0][i];
        }
        for (int i = 0; i < NBLOCKS; i ++) {
            //printf("%d ", marked[i]);
            if (!marked[i]) cout<<"block "<<i<< " unmarked\n";
        }
        for (int i = 0; i < NTHREADS; i ++) {
            //printf("indices = %0.f nrounds = %d\n", indices[i].load(), nrounds[i]);
        }
       
#endif
        // check is block handled
        for (int i = 0; i < NTHREADS; i ++) {
            int core_id = i;
            for (int round = 0; round < nrounds[core_id]; round ++) {
                int block_id = block_matrix[round*NTHREADS+core_id];
                if (block_handled[block_id] == 0) {
                    printf("coreid = %d, round = %d, blockid = %d\n", core_id, round, block_id);
                }
            }
        }
        int count = 0;
        for (int i = 0; i < NBLOCKS; i ++) {
            if (block_handled[i] == 0) {
                printf("unhandled block id = %d\n", i);
                count ++;
            }
        }
        cout  << count << endl;
        // check scheds
        for (int i = 0; i < n; i ++) {
            if ( (sched_offsets[i] < 0 || sched_offsets[i] > m) || sched_offsets[i] != nodepointer[0][i] ) {
                printf("sched_offsets[%d] = %lld nodepointer[%d] = %lld\n", i, sched_offsets[i], i, nodepointer[0][i]);
                abort();
            }
            long long start = sched_offsets[i];
            long long end = sched_offsets[i+1];
            if (end > m ) {printf("node %d end = %lld\n",i, end); exit(0); }
            for (long long j = start; j < end; j ++) {
                    
                //if (sched_edges[j] >= n || sched_edges[j] < 0 || sched_edges[j] != edges[0][j]) {
                //    printf("sched_edges[%lld] = %d, edges[%lld] = %d", j, sched_edges[j],j, edges[0][j]);
                //    abort();
                //}
            }   
        }
        
        sched_offsets[n] = m;
        for (int i = 0; i <= n; i ++) {
            sched_offsets[i] = nodepointer[0][i];
        }
        #if !defined(TWITTER) && !defined(FRIENDSTER)
        #endif
        #ifdef CACHE_SCHEDULE
        nodepointer[0] = sched_offsets;
        edges[0] = sched_edges;
        #endif
        sched_offsets = NULL;
        //sched_edges = NULL;
        cout << "   numa aware residue array loaded" << endl;
        //delete[] marked;
        //delete[] sums;
        //delete[] block_matrix;
    }
    static long lrandp() {
        if (sizeof(int) < sizeof(long))
            return (((long)global_rng.get()) << (sizeof(int) * 8)) |
               global_rng.get();

        return global_rng.get();
    }
    // Compare and Swap
    void atomic_add(atomic<double> &a, double b) {
        auto current = a.load();
        while (!a.compare_exchange_weak(current, current + b));
            ;
    }
    
    // This version has return value that is used for discarding locks
    inline double Atomic_add(atomic<double> &a, double b) {
        while (1) {
            auto current = a.load();
            if (a.compare_exchange_weak(current, current + b)) {
                return current;
            } 
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
    inline double Atomic_add(atomic<int> &a, int b) {
        while (1) {
            auto current = a.load();
            if (a.compare_exchange_weak(current, current + b)) {
                return current;
            } 
        }
    }
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
};

#endif
