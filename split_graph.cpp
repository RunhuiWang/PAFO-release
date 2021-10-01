
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include <memory.h>

#include "covergraph.h"

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

int main(int argc, char *argv[]) {

    cout << Green << "--------------start------------" << get_current_time_str() << Reset << endl;
    string defaultdata = "../../../uqswan22/datasets/twitter";
    string datagraph = defaultdata;;
    const char* num_of_workers;
    int pouch_size_cmd = 0;
    int rwscale = 0;
    for(int i=0; i<argc; i++){
        if(string(argv[i]) == "--dataset"){
            datagraph = string(argv[i+1]);   
        }
    }
    Graph graph = Graph(datagraph);
    graph.split_graph(4);
    //graph.free_memory();
    cout << Red << "--------------stop------------" << get_current_time_str() << Reset << endl << endl << endl;
    return 0;
} 
