#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <iostream>
#include <memory.h>
#include <vector>

using namespace std;
//extern Graph graph;
class MinBlk{
    public:
    double coverage;
    int idx;
    void operator=(MinBlk b2) {
        coverage = b2.coverage;
        idx = b2.idx;
    }
    bool operator>(const MinBlk& b2) {
        return coverage > b2.coverage;
    }
    MinBlk(double cvg, int id) {
        coverage = cvg;
        idx = id;
    }
    MinBlk() {}
};
template <class T>
class Min{
    public:
    T number;
    Min() {
        number = MinBlk(INT_MAX,0);
    }
    Min(T value) {
        number = value;
    }
};

template<class T>
class MinView {
public:
    Min<T>* num;    
    MinView() {
        num = new Min<T>();
    }
    MinView(T value) {
        num = new Min<T>(value);
    }

    void increase(T inc) {
        num->number += inc;
    }

    void operator+=(T inc) {
        num->number += inc;
    }

    bool operator>(T b2) {
        return num->number > b2;
    }
    MinView<T>& operator=(T value) {
        num->number = value;
    }

    T get_value() const { 
        return num->number; 
    }
    Min<T>* view_get_value() const { 
        return num; 
    }
    void set_view_value(T b2) {
        num->number = b2;
    }

};

// Monoid class.
template<class T>
struct MinMonoid : public cilk::monoid_base<Min<T>*, MinView<T>> {

    static void identity(MinView<T>* view) {
        view->num = new Min<T>();
    }
    
    static void reduce(MinView<T>* left, MinView<T>* right) {
        T l = left->num->number;
        T r = right->num->number;
        //printf("In max reducer: left is %d, right is %d\n", l, r);
        if (l > r) {
            left->num->number = r;
        }
    }
};

