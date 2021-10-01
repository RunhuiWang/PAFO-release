#ifndef _BENCH_GRAPH_IO
#define _BENCH_GRAPH_IO

#include <iostream>
#include <stdint.h>
#include <cstring>
#include "utils/parallel.h"
#include "utils/utils.h"
//#include "utils/blockRadixSort.h"
//#include "utils/quickSort.h"
using namespace std;

// **************************************************************
//    EDGE ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct edge {
  intT u;
  intT v;
  edge(intT f, intT s) : u(f), v(s) {}
};

template <class intT>
struct edgeArray {
  edge<intT>* E;
  intT numRows;
  intT numCols;
  intT nonZeros;
  void del() {delete[] E;}
  edgeArray(edge<intT> *EE, intT r, intT c, intT nz) :
    E(EE), numRows(r), numCols(c), nonZeros(nz) {}
  edgeArray() {}
};

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class intT>
struct vertex {
  intT* Neighbors;
  intT degree;
  void del() {delete[] Neighbors;}
  vertex(intT* N, intT d) : Neighbors(N), degree(d) {}
  vertex(){}
};


template <class intT>
struct paGraph {
  vertex<intT> *V=NULL;
  intT n;
  long long m;
  long long * offsets=NULL;
  int* edges=NULL;
  paGraph(){}
  paGraph(vertex<intT>* VV, intT nn, long long mm) 
    : V(VV), n(nn), m(mm), allocatedInplace(NULL) {}
  paGraph(vertex<intT>* VV, intT nn, long long mm, intT* ai, long long* attrr )
    : V(VV), n(nn), m(mm), edges(ai), offsets(attrr){
    
  }
  paGraph copy() {
    prinf("paGraph copy() calledn\n");
    //vertex<intT>* VN = newA(vertex<intT>,n);
    //intT* _allocatedInplace = newA(intT,n+m+2);
    vertex<intT>* VN = new vertex<intT>[n];
    intT* _allocatedInplace = new intT[n+m+2];

    _allocatedInplace[0] = n;
    _allocatedInplace[1] = m;
    intT* Edges = _allocatedInplace+n+2;
    intT k = 0;
    for (intT i=0; i < n; i++) {
      _allocatedInplace[i+2] = allocatedInplace[i+2];
      VN[i] = V[i];
      VN[i].Neighbors = Edges + k;
      for (intT j =0; j < V[i].degree; j++) 
	Edges[k++] = V[i].Neighbors[j];
    }
    return paGraph(VN, n, m, _allocatedInplace);
  } 
  void del() {
    printf("paGraph del() called\n");
    delete[] V;
    delete[] offsets;
    delete[] edges;
    printf("paGraph del() ended\n");
  }
};


// **************************************************************
//    GRAPH UTILITIES
// **************************************************************

template <class intT>
paGraph<intT> paGraphFromEdges(edgeArray<intT> EA, bool makeSym) {
  edgeArray<intT> A;
  if (makeSym) A = makeSymmetric<intT>(EA);
  else {  // should have copy constructor
    edge<intT> *E = newA(edge<intT>,EA.nonZeros);
    parallel_for (intT i=0; i < EA.nonZeros; i++) E[i] = EA.E[i];
    A = edgeArray<intT>(E,EA.numRows,EA.numCols,EA.nonZeros);
  }
  intT m = A.nonZeros;
  intT n = max<intT>(A.numCols,A.numRows);
  intT* offsets = newA(intT,n*2);
  intSort::iSort(A.E,offsets,m,n,getuF<intT>());
  intT *X = newA(intT,m);
  vertex<intT> *v = newA(vertex<intT>,n);
  parallel_for (intT i=0; i < n; i++) {
    intT o = offsets[i];
    intT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].degree = l;
    v[i].Neighbors = X+o;
    for (intT j=0; j < l; j++) {
      v[i].Neighbors[j] = A.E[o+j].v;
    }
  }
  A.del();
  delete[] offsets;
  return paGraph<intT>(v,n,m,X);
}


// **************************************************************
//    BASIC I/O
// **************************************************************
namespace benchIO {
  using namespace std;

  // A structure that keeps a sequence of strings all allocated from
  // the same block of memory
  struct words {
    long n; // total number of characters
    char* Chars;  // array storing all strings
    long m; // number of substrings
    char** Strings; // pointers to strings (all should be null terminated)
    words() {}
    words(char* C, long nn, char** S, long mm)
      : Chars(C), n(nn), Strings(S), m(mm) {}
    void del() {delete[] Chars; delete[] Strings;}
  };
 
  inline bool isSpace(char c) {
    switch (c)  {
    case '\r': 
    case '\t': 
    case '\n': 
    case 0:
    case ' ' : return true;
    default : return false;
    }
  }

  struct toLong { long operator() (bool v) {return (long) v;} };

  // parallel code for converting a string to words
  words stringToWords(char *Str, long n) {

    set_Workers(40);
    parallel_for (long i=0; i < n; i++) 
      if (isSpace(Str[i])) Str[i] = 0; 
    // mark start of words
    bool *FL = newA(bool,n);
    FL[0] = Str[0];
    parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];
    //cout << "-------stringToWords start of words marked" << endl;
    // offset for each start of word
    _seq<long> Off = sequence::packIndex<long>(FL, n);
    long m = Off.n;
    long *offsets = Off.A;
    //cout << "-------stringToWords offset for words set" << endl;

    // pointer to each start of word
    char **SA = newA(char*, m);
    //cout << "-------------new array SA allocated m = " << m << endl;
    int rounds = 1;
    if (m > 5000000000) rounds = 50;
    for (int i = rounds-1; i >= 0; i --) {
        long start = (m/rounds) * i;
        long end = (m/rounds) * (i+1);
        if (i == rounds-1) end = m;
        parallel_for (long j=start; j<end; j++) SA[j] = Str+offsets[j];
        //printf("------------round %d\n", i);
    }
    //cout << "-------stringToWords pointer to each start of work set" << endl;

    set_Workers(1);

    delete[] offsets; delete[] FL;
    return words(Str,n,SA,m);
  }

  int writeStringToFile(char* S, long n, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file.write(S, n);
    file.close();
    return 0;
  }

  inline int xToStringLen(long a) { return 21;}
  inline void xToString(char* s, long a) { sprintf(s,"%ld",a);}

  inline int xToStringLen(unsigned long a) { return 21;}
  inline void xToString(char* s, unsigned long a) { sprintf(s,"%lu",a);}

  inline uint xToStringLen(uint a) { return 12;}
  inline void xToString(char* s, uint a) { sprintf(s,"%u",a);}

  inline int xToStringLen(int a) { return 12;}
  inline void xToString(char* s, int a) { sprintf(s,"%d",a);}

  inline int xToStringLen(double a) { return 18;}
  inline void xToString(char* s, double a) { sprintf(s,"%.11le", a);}

  inline int xToStringLen(char* a) { return strlen(a)+1;}
  inline void xToString(char* s, char* a) { sprintf(s,"%s",a);}

  template <class A, class B>
  inline int xToStringLen(pair<A,B> a) { 
    return xToStringLen(a.first) + xToStringLen(a.second) + 1;
  }
  template <class A, class B>
  inline void xToString(char* s, pair<A,B> a) { 
    int l = xToStringLen(a.first);
    xToString(s,a.first);
    s[l] = ' ';
    xToString(s+l+1,a.second);
  }

  struct notZero { bool operator() (char A) {return A > 0;}};

  template <class T>
  _seq<char> arrayToString(T* A, long n) {
    long* L = newA(long,n);
    {parallel_for(long i=0; i < n; i++) L[i] = xToStringLen(A[i])+1;}
    long m = sequence::scan(L,L,n,addF<long>(),(long) 0);
    char* B = newA(char,m);
    parallel_for(long j=0; j < m; j++) 
      B[j] = 0;
    parallel_for(long i=0; i < n-1; i++) {
      xToString(B + L[i],A[i]);
      B[L[i+1] - 1] = '\n';
    }
    xToString(B + L[n-1],A[n-1]);
    B[m-1] = '\n';
    delete[] L;
    char* C = newA(char,m+1);
    long mm = sequence::filter(B,C,m,notZero());
    C[mm] = 0;
    delete[] B;
    return _seq<char>(C,mm);
  }

  template <class T>
  void writeArrayToStream(ofstream& os, T* A, long n) {
    long BSIZE = 1000000;
    long offset = 0;
    while (offset < n) {
      // Generates a string for a sequence of size at most BSIZE
      // and then wrties it to the output stream
      _seq<char> S = arrayToString(A+offset,min(BSIZE,n-offset));
      os.write(S.A, S.n);
      S.del();
      offset += BSIZE;
    }    
  }

  template <class T>
    int writeArrayToFile(string header, T* A, long n, char* fileName) {
    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file << header << endl;
    writeArrayToStream(file, A, n);
    file.close();
    return 0;
  }

  _seq<char> readStringFromFile(const char *fileName) {
    ifstream file (fileName, ios::in | ios::binary | ios::ate);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      abort();
    }
    long end = file.tellg();
    file.seekg (0, ios::beg);
    long n = end - file.tellg();
    char* bytes = newA(char,n+1);
    file.read (bytes,n);
    file.close();
    return _seq<char>(bytes,n);
  }
};


// **************************************************************
//    GRAPH I/O
// **************************************************************
using namespace benchIO;

template <class intT>
int xToStringLen(edge<intT> a) { 
  return xToStringLen(a.u) + xToStringLen(a.v) + 1;
}

template <class intT>
void xToString(char* s, edge<intT> a) { 
  int l = xToStringLen(a.u);
  xToString(s, a.u);
  s[l] = ' ';
  xToString(s+l+1, a.v);
}

namespace benchIO {
  using namespace std;

  string AdjGraphHeader = "AdjacencyGraph";
  string WghAdjGraphHeader = "WeightedAdjacencyGraph";

  template <class intT>
  int writeGraphToFile(paGraph<intT> G, char* fname) {
    long m = G.m;
    long n = G.n;
    long totalLen = 2 + n + m;
    intT *Out = newA(uintT, totalLen);
    Out[0] = n;
    Out[1] = m;
    parallel_for (long i=0; i < n; i++) {
      Out[i+2] = G.V[i].degree;
    }
    long total = sequence::scan(Out+2,Out+2,n,addF<intT>(),(intT)0);
    for (long i=0; i < n; i++) {
      intT *O = Out + (2 + n + Out[i+2]);
      vertex<intT> v = G.V[i];
      for (long j = 0; j < v.degree; j++) 
	O[j] = v.Neighbors[j];
    }
    int r = writeArrayToFile(AdjGraphHeader, Out, totalLen, fname);
    delete[] Out;
    return r;
  }


  template <class intT>
  edgeArray<intT> readSNAP(char* fname) {
    _seq<char> S = readStringFromFile(fname);
    char* S2 = newA(char,S.n);
    //ignore starting lines with '#' and find where to start in file 
    long k=0;
    while(1) {
      if(S.A[k] == '#') {
	while(S.A[k++] != '\n') continue;
      }
      if(k >= S.n || S.A[k] != '#') break; 
    }
    parallel_for(long i=0;i<S.n-k;i++) S2[i] = S.A[k+i];
    S.del();

    words W = stringToWords(S2, S.n-k);
    long n = W.m/2;
    edge<intT> *E = newA(edge<intT>,n);
    {parallel_for(long i=0; i < n; i++)
      E[i] = edge<intT>(atol(W.Strings[2*i]), 
		  atol(W.Strings[2*i + 1]));}
    W.del();

    long maxR = 0;
    long maxC = 0;
    for (long i=0; i < n; i++) {
      maxR = max<intT>(maxR, E[i].u);
      maxC = max<intT>(maxC, E[i].v);
    }
    long maxrc = max<intT>(maxR,maxC) + 1;
    return edgeArray<intT>(E, maxrc, maxrc, n);
  }


  template <class intT>
  paGraph<intT> readGraphFromFile(const char* fname) {
    cout << "read string from file..."<< endl;
        //numa_set_interleave_mask(numa_all_nodes_ptr);
    _seq<char> S = readStringFromFile(fname);
        //numa_set_localalloc();
    cout << "read string from file done"<< endl;

    cout << "string to words..." << endl;
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != AdjGraphHeader) {
      cout << "Bad input file: missing header: " << AdjGraphHeader << endl;
      abort();
    }
    cout << "string to words done" << endl;

    long len = W.m -1;
    set_Workers(80);
    long n = atol(W.Strings[1]);
    long long m = atol(W.Strings[2]);
    long long * offsets = new long long [n];
    memset(offsets, 0, sizeof(long long)*n);
    {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i+3]);}
    uintT * In = newA(uintT, len-n-2);
    {parallel_for(long i=0; i < len-n-2; i++) In[i] = atol(W.Strings[i+n+3]);}
    W.del();
    printf("atol finished\n");
    //long n = In[0];
    //long long m = In[1];

    printf("n=%ld m=%lld\n", n, m);
    if (len != n + m + 2) {
      cout << "Bad input file: length = "<<len<< " n+m+2 = " << n+m+2 << endl;
      abort();
    }
    vertex<intT> *v = newA(vertex<intT>,n);
    uintT* edges = In;

    parallel_for (uintT i=0; i < n; i++) {
      uintT o = offsets[i+2];
      uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
      v[i].degree = l;
      v[i].Neighbors = (intT*)(edges+o);
    }
    set_Workers(1);
    return paGraph<intT>(v,(intT)n,(long long)m,(intT*)In, offsets);
  }

};

#endif // _BENCH_GRAPH_IO
