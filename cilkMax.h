#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <iostream>
#include <memory.h>
#include <vector>

using namespace std;
extern Graph graph;
template <class T>
class Max{
    public:
    T number;
    Max() {
        number = 0;
    }
    Max(T value) {
        number = value;
    }
};

template<class T>
class MaxView {
public:
    Max<T>* num;    
    MaxView() {
        num = new Max<T>();
    }
    MaxView(T value) {
        num = new Max<T>(value);
    }

    void increase(T inc) {
        num->number += inc;
    }

    void operator+=(T inc) {
        num->number += inc;
    }
    MaxView<T>& operator=(T value) {
        num->number = value;
    }

    T get_value() const { 
        return num->number; 
    }
    Max<T>* view_get_value() const { 
        return num; 
    }

};

// Monoid class.
template<class T>
struct MaxMonoid : public cilk::monoid_base<Max<T>*, MaxView<T>> {

    static void identity(MaxView<T>* view) {
        view->num = new Max<T>();
    }
    
    static void reduce(MaxView<T>* left, MaxView<T>* right) {
        T l = left->num->number;
        T r = right->num->number;
        //printf("In max reducer: left is %d, right is %d\n", l, r);
        if (l < r) {
            left->num->number = r;
        }
    }
};

