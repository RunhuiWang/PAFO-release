#ifndef __SACK_H__
#define __SACK_H__

//#include <vector>
//#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <map>
#include <memory.h>

extern Graph graph;
extern double rmax;
const int INITIAL_SIZE = 128;
const int INCREASE_TIMES = 8;
// forward declaration
class SackMonoid;

class Queue{

public:
    
    int size;
    int capacity;
    int* queue;

    Queue() {
        size = 0;
        capacity = INITIAL_SIZE;
        queue = new int[capacity];
    }

    Queue(int s, int* q) {
        size = s;
        capacity = s;
        queue = q;
    }

    void enqueue(int node) {
        
        // expand array when it's full
        if (size == capacity) {
            
            capacity *= INCREASE_TIMES;
            int *temp = new int[capacity];
            memcpy(temp, queue, sizeof(int) * size);
            
            // Deletion will never happen after Queue's split 
            // so deletion here is safe
            delete[] queue;
            queue = temp;

        }
        queue[size++] = node;
    }

    Queue* split() {
        
        // calculate the sizes after splitting
        int size_q = size / 2;
        size -= size_q;

        Queue* q = new Queue(size_q, queue+size);

        return q;

    }
    int operator[] (int index) {
        return queue[index];
    }
    int at(int index) {
        return queue[index];
    }
    int get_size() { return size;}

    void clear() {
        size = 0;
    }

/*    // avoid lock
    void remove(int position) {
        queue[position] = -1;
    }
*/    
    ~Queue() {
        if (queue) {
            delete[] queue;
        }
    }
};

class Sack {

public:

    // pairs of node index and its residue
    //map<int, double> residue;    
    // an array for residue
    vector<double>* residue;
    vector<int>* isbig;
    // the residue vector may be used for next iteration
    bool is_occupied;
    Queue* queue;
    Queue* tiny_queue;

    // depricated
    //map<int, bool> inserted;

    Sack() {
        //residue.assign(graph.n, 0.0);
        residue = new vector<double>(0,0.0);
        isbig = new vector<int>(0, 0);
        queue = new Queue();
        tiny_queue = new Queue();
        is_occupied = false;
        //inserted = map<int, bool>();
    }

    ~Sack() {
        
        printf("in sack destructor is_occupied = %d\n ", is_occupied);
        // only delete residue vector when it is not used for next iteration
        if (!is_occupied && residue) {
        printf("in sack destructor is_occupied = %d\n ", is_occupied);
            delete residue;
            if (isbig) delete isbig;
        }
        
        if (queue) {
            delete queue;
        }
        
        if (tiny_queue) {
            delete tiny_queue;
        }

    }

    int size() {
        return queue->get_size();
    }
    
    vector<double>* get_residue() {
        return residue;
    }

    Queue* get_queue() { 
        return queue; 
    }
    
    void set_tiny_queue(Queue* q) { 
        tiny_queue = q;
    }

    void set_queue(Queue* q) { 
        queue = q;
    }


    // lock the residue vector so that the memory will not be freed
    void occupy_residue_vector() {
        is_occupied = true;
    }

    void insert(int key, double value, bool big) {
       
        //printf("sack insert called %d\n", key);
        //if (key != 96787) exit(1);
        // only assign space when first insertion
        // so as to save space when split
        if (residue->empty()) {
            residue->assign(graph.n, 0.0);
            isbig->assign(graph.n, 0);
        }
        // reduce redundancy, if big equals false then do not enqueue
        if (residue->operator[](key) == 0.0) {
            if ( big ) {
                queue->enqueue(key);
                isbig->operator[](key) = 1;
            } else {
                tiny_queue->enqueue(key);
            }
        }
        residue->operator[](key) += value;
        if (isbig->operator[](key) == 0) {
            if (residue->operator[](key) >= graph.thres[key]) {
                isbig->operator[](key) = 1;
                //tiny_queue->enqueue(-key);
                queue->enqueue(key);
            }
        }
//        printf("sack insert %.15f on node %d\n", value, key);
    }

    Sack* split() {
        Sack* sack = new Sack();
        sack->get_half_queue(this);
        return sack;
    }
    
    /*
    bool is_inserted(int key) {
        if (inserted.find(key) == inserted.end()) {
            
            inserted.insert(pair<int, bool>(key, true));
            return false;
        }
        return true;
    }*/

    void get_half_queue(Sack* o_sack) {

        // split the origin queue
        Queue* q = o_sack->get_queue()->split();
        
        // set half queue to this sack
        set_queue(q);
    }

    void clear() {
        residue->clear();
        queue->clear();
    }

};

//View class

class SackView {
    
    // for the identity and reduce functions
    friend class SackMonoid;
    Sack *sack;

public:
   
    SackView() {
        
        sack = new Sack();
    }
    
    int size() {
        return sack->size();
    }
    
    /*
        if enqueue is false then just increase the residude value and
        do not enqueue the node
    */
    void insert(int key, double value, bool isbig) {
        
        //sack->insert(key, value, enqueue);
        // this is a dequeue operation
//        if (key < 0) {
//                
//            sack->residue->operator[](-key) -= value; 
//            return;
//        }
        // only assign space when first insertion
        // so as to save space when split
        if (sack->residue->empty()) {
            sack->residue->assign(graph.n, 0.0);
            sack->isbig->assign(graph.n, 0);
        }
        // 
        // reduce redundancy, if big equals false then do not enqueue
        if (sack->residue->operator[](key) == 0.0) { 
            if ( isbig ) {
                sack->isbig->operator[](key) = 1;
                sack->queue->enqueue(key);
            } else {
                sack->tiny_queue->enqueue(key);
//                if (key == check) {
//                    printf("%d is tiny-enqueued, out degree=%d isbig=%d %.15f %.15f\n", check,
//                      graph.g[key].size(), sack->isbig->operator[](key), value,
//                      sack->residue->operator[](key));
//                }
            }
        }
        sack->residue->operator[](key) += value;
//        if (key == check) {
//            printf("   %d = %.15f\n", key, sack->residue->operator[](key));
//        }

        //printf("isbig on node %d is %d\n", key, sack->isbig->operator[](key));
        if (sack->isbig->operator[](key) == 0) {
            if ((sack->residue->operator[](key)/graph.g[key].size()) >= rmax) {
                sack->isbig->operator[](key) = 1;
                //sack->tiny_queue->enqueue(-key);
                sack->queue->enqueue(key);
                if (sack->residue->at(key) > 1) {
                    printf("residue is %.10ferror!!!!!!!!!!!!!!!!!!!!!!\n",
                        sack->residue->at(key));
                    exit(1);
                }
//                if (key == check) {
//                    printf("%d is big-enqueued out degree is %d isbig=%d\n",
//                       check,graph.g[key].size(), sack->isbig->operator[](key));
//                    printf("rmax = %.15f, r/nout=%.15f\n",rmax, sack->residue->operator[](key) / graph.g[key].size());
//                }
            }
        }
    }
    
    Sack* get_value() const{
        return sack;
    }

    Sack* view_get_value() const{
        return sack;
    }

    // make the sack empty
    void clear() {
        sack->clear();
    }


};

struct SackMonoid : public cilk::monoid_base<Sack*, SackView> {

    // constructs the identity value into the uninitialized memory
    static void identity(SackView* view) {
        //view->sack->residue.assign(graph.n, 0.0);
        view->sack = new Sack();
    }

    // merge two views
    // move the content in the right view to the left
    // after reduce right view is empty
    static void reduce(SackView* left, SackView* right) {
        
        // to reduce the operations of merge
//        int size_left = left->size();
//        int size_right= right->size();
//
//        if (size_left < size_right) {
//            SackView* temp = right;
//            right = left;
//            left = temp;
//        }

        // merge two sacks
        Sack* sack0 = left->get_value();
        Sack* sack1 = right->get_value();
        //vector<int>::size_type i;

        // depricated
        //sack0->inserted.clear();


        //printf("in monoid reduce() left size is %d, right size is %d\n"
        //        ,sack0->size(), sack1->size());
//        printf("left 165070 is %.15f, right 165070 is%.15f\n"
//                ,sack0->residue->at(165070), sack1->residue->at(165070));
        for ( int i = 0; i < sack1->queue->size; i++) {
            int node = sack1->queue->at(i);
            double residue = sack1->residue->at(node);

            //sack0->residue->at(node) = residue;
            //sack0->queue->enqueue(node);
            //if ( !sack0->is_inserted(node)) {
            sack0->insert(node, residue, true);
            //}
        }

        //printf("left tiny size is %d\n", sack1->tiny_queue->size);
        for (int i = 0; i < sack1->tiny_queue->size; i ++) {
            int node = sack1->tiny_queue->at(i);
            double residue = sack1->residue->at(node);
            int isbig = sack1->isbig->operator[](node);
            if (isbig == 0)
                sack0->insert(node, residue, false);
        }
        right->clear();
//        printf("left 165070 is %.15f, right 165070 is%.15f\n"
//                ,sack0->residue->at(165070), sack0->residue->at(165070));
//        printf("end monoid reduce() left size is %d, right size is %d\n"
//                ,sack0->size(), sack1->size());

    }
    

};
#endif // __SACK_H__
