#include <fstream>
#include <iostream>
#include <string.h>
#include "paLoad.h"
 
int main() {
    string data_folder = "/media/bigdata/uqrwan14/graphdata/friendster/";
    string rw_index_file_str = data_folder + "/rwidx3_adj.txt";
    paGraph<int> paidx = readGraphFromFile<int>(rw_index_file_str.c_str());
}
