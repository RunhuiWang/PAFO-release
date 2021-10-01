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

void atomic_add(atomic<uint8_t> &a, uint8_t b) {
    auto current = a.load();
    while (!a.compare_exchange_weak(current, current + b));
        ;
}
uint8_t Atomic_add(atomic<uint8_t> &a, uint8_t b) {
    while (1) {
        auto current = a.load();
        if (a.compare_exchange_weak(current, current + b)) {
            return current;
        } 
    }
}
int main() {
    //uint8_t a = 0;
    atomic<uint8_t> a(0);
    for (int i = 0; i < 600; i ++) {
        printf("%d ",Atomic_add(a, 1) );
    }
    cout << endl;
}
