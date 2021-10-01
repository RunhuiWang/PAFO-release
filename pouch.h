
int INITIAL_P_SIZE = 16;
class Queue{

public:
    
    int size;
    int capacity;
    int* queue;
    int workload;
    Queue() {
        size = 0;
        capacity = INITIAL_P_SIZE;
        queue = NULL;
        workload = 0;
    }

    Queue(int s, int* q) {
        size = s;
        capacity = s;
        queue = q;
    }

    void enqueue(int node) {

        if (queue == NULL) {
            queue = new int[INITIAL_P_SIZE]; 
        }
        // expand array when it's full
        if (size == capacity) {
            
            capacity *= 2;
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

    ~Queue() {
        if (queue) {
            delete[] queue;
        }
    }
};

struct Subnode {
    int node;
    int start;
    int end;
    double residue;
};

class Pouch {
    
    int load;
public:
    
    Queue *queue;
    Subnode subnode1;
    Subnode subnode2;

    Pouch() {
        queue = new Queue();
        subnode1 = {-1, -1, -1};
        subnode2 = {-1, -1, -1};
        load = 0;
    }

    void enqueue(int node, int size) {
        queue->enqueue(node);
        load += size;
    }

    void set_subnode_1(int node, int start, int end, double residue) {
        subnode1.node = node;
        subnode1.start = start;
        subnode1.end = end;
        subnode1.residue = residue;
        load += (end-start);
    }

    void set_subnode_2(int node, int start, int end, double residue) {
        subnode2.node = node;
        subnode2.start = start;
        subnode2.end = end;
        subnode2.residue = residue;
        load += (end-start);
    }

    int get_load() { return load;}
    //~Pouch() {
    //    delete queue;
    //}

};
