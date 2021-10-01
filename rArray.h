#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <iostream>
#include <memory.h>
#include <vector>

using namespace std;
extern Graph graph;
class RArray {
    public:
    double* array;
    RArray() {
        array = new double[graph.n];
        memset(array, 0.0, sizeof(double)*graph.n);
    }
    ~RArray() {
        if (array) delete[] array;
    }
};

class RArrayView {
public:
    RArray* residue;    
    RArrayView() {
        residue = new RArray();
        printf("new view size is %d\n", graph.n);
    }

    void add_value(int idx, double value) {
        residue->array[idx] += value;
    }

    void merge_value(RArrayView* other) { 
        double* array = other->get_value()->array;
        for (int i = 0; i < graph.n; i ++) {
            residue->array[i] += array[i];
        }
    }

    RArray* get_value() const { 
        return residue; 
    }
    RArray* view_get_value() const { 
        return residue; 
    }

};

// Monoid class.
//
struct RArrayMonoid : public cilk::monoid_base<RArray*, RArrayView> {

    static void identity(RArrayView* view) {
        view->residue = new RArray();
    }
    
    static void reduce(RArrayView* left, RArrayView* right) {
        left->merge_value(right);
    }
};

