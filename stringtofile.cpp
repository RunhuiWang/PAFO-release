#include <string>
#include <iostream>
#include <fstream>
using namespace std;

int main() {
    string test = "begin";
    int n = 100;
    for (int i = 0; i < n; i ++) {
        test+= itoa(i)+"\n";
    }
    test += "\n";
    cout << test;
    ofstream outfile;
    outfile.open("stringOuputTest.txt");
    outfile << test;
    outfile.close();

    return 0;
}
