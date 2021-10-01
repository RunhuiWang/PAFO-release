#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>


using namespace std;

class Graph {
    double ALPHA_DEFAULT = 0.2;
public:

    vector< vector<int> > g;
    vector< vector<int> > gr;
    string data_folder;

    //vector<double> global_ppr;

    // node rank[100] = 0, means node 100 has first rank
    vector<int> node_rank;
    // node_score[0]
    vector<double> node_score;

    //node order 0 = [100, 34.5], most important node is 100 with score 34.5
    vector< pair<int, double> > node_order;
    vector<int> loc;
    
    vector<double> thres;

    // the tele ratio for random walk
    double alpha;

    static bool cmp(const pair<int, double> &t1, const pair<int, double> &t2) {
        return t1.second > t2.second;
    }

    int n;
    long long m;

    // constructor
    Graph(string data_folder) {
        this->data_folder = data_folder;
        this->alpha = ALPHA_DEFAULT;
        init_graph();
        thres.assign(n, 0.0);
        cout << "init graph n: " << this->n << " m: " << this->m << endl;
    }

    void compute_threshold(double rmax) {
        
        for (int i = 0; i < n; i ++) {
            thres[i] = rmax * g[i].size();
        }
    
    }

    // default constructor
    Graph() {
        // do nothing
    }
    void init_nm() {
        string attribute_file = data_folder + "/" + "attribute.txt";
        assert_file_exist("attribute file", attribute_file);
        ifstream attr(attribute_file.c_str());
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
    }

    void init_graph() {
        init_nm();

        g = vector<vector<int> >(n, vector<int>());
        printf("graph init\n");
        gr = vector<vector<int> >(n, vector<int>());
        string graph_file = data_folder + "/graph.txt";
        //assert_file_exist("graph file", graph_file);
        FILE *fin = fopen(graph_file.c_str(), "r");
        int t1, t2;

        //DEBUG
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            assert(t1 < n);
            assert(t2 < n);
            if(t1 == t2) continue;
            g[t1].push_back(t2);
            gr[t2].push_back(t1);
        }
    }
    
    void check_graph() {
        for (int i =0 ; i < n; i ++) {
            printf("%d\n", (int)g[i].size());
        }
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
};

#endif // __GRAPH_H__
