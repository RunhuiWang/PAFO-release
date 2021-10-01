#ifndef __PAGERANK_H__
#define __PAGERANK_H__
//#include "mylib.h"
//#include "config.h"
//#include "covergraph.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <memory.h>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>

const int VectorDefaultSize=20;


using namespace std;
template <typename _T>
class iVector
{
public:
    unsigned int m_size;
    _T* m_data;
    unsigned int m_num;

    void free_mem()
    {
        delete[] m_data;
    }

    iVector()
    {
        //printf("%d\n",VectorDefaultSize);
        m_size = VectorDefaultSize;
        m_data = new _T[VectorDefaultSize];
        m_num = 0;
    }
    iVector( unsigned int n )
    {
        if ( n == 0 )
        {
            n = VectorDefaultSize;
        }
//      printf("iVector allocate: %d\n",n);
        m_size = n;
        m_data = new _T[m_size];
        m_num = 0;
    }
    void push_back( _T d )
    {
        if ( m_num == m_size )
        {
            re_allocate( m_size*2 );
        }
        m_data[m_num] = d ;
        m_num++;        
    }
    void push_back( const _T* p, unsigned int len )
    {
        while ( m_num + len > m_size )
        {
            re_allocate( m_size*2 );
        }
        memcpy( m_data+m_num, p, sizeof(_T)*len );
        m_num += len;
    }

    void re_allocate( unsigned int size )
    {
        if ( size < m_num )
        {
            return;
        }
        _T* tmp = new _T[size];
        memcpy( tmp, m_data, sizeof(_T)*m_num );
        m_size = size;
        delete[] m_data;
        m_data = tmp;
    }
    void Sort()
    {
        if ( m_num < 20 )
        {
            int k ;
            _T tmp;
            for ( int i = 0 ; i < m_num-1 ; ++i )
            {
                k = i ;
                for ( int j = i+1 ; j < m_num ; ++j )
                    if ( m_data[j] < m_data[k] ) k = j ;
                if ( k != i )
                {
                    tmp = m_data[i];
                    m_data[i] = m_data[k];
                    m_data[k] = tmp;
                }
            }
        }
        else sort( m_data, m_data+m_num );
    }
    void unique()
    {
        if ( m_num == 0 ) return;
        Sort();
        unsigned int j = 0;
        for ( unsigned int i = 0 ; i < m_num ; ++i )
            if ( !(m_data[i] == m_data[j]) )
            {
                ++j;
                if ( j != i ) m_data[j] = m_data[i];
            }
        m_num = j+1;
    }
    int BinarySearch( _T& data )
    {
        for ( int x = 0 , y = m_num-1 ; x <= y ; )
        {
            int p = (x+y)/2;
            if ( m_data[p] == data ) return p;
            if ( m_data[p] < data ) x = p+1;
            else y = p-1;
        }
        return -1;
    }
    void clean()
    {
        m_num = 0;
    }
    void assign( iVector& t )
    {
        m_num = t.m_num;
        m_size = t.m_size;
        delete[] m_data;
        m_data = t.m_data;
    }

    bool remove( _T& x )
    {
        for ( int l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;

            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
            else if ( m_data[m] < x ) l = m+1;
            else r = m;
        }
        return false;
    }

    void sorted_insert( _T& x )
    {
        if ( m_num == 0 )
        {
            push_back( x );
            return;
        }

        if ( m_num == m_size ) re_allocate( m_size*2 );

        int l,r;

        for ( l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;
            if ( m_data[m] < x ) l = m+1;
            else r = m;
        }

        if ( l < m_num && m_data[l] == x )
        {
            //printf("Insert Duplicate....\n");
            //cout<<x<<endl;
    //      break;
        }
        else
        {
            if ( m_num > l )
            {
                memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
            }
            m_num++;
            m_data[l] = x;
        }
    }

    bool remove_unsorted( _T& x )
    {
        for ( int m = 0 ; m < m_num ; ++m )
        {
            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
        }
        return false;
    }

    _T& operator[]( unsigned int i )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[i];
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //close range check for [] in iVector if release

};

template <typename _T>
struct iMap
{   
    _T* m_data;
    int m_num;
    int cur;
    iVector<int> occur;
    _T nil; 
    iMap()
    {
        m_data = NULL;
        m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
        //nil = 1073741834;
    }
    iMap(int size){
        initialize(size);
    }
    void free_mem()
    {
        delete[] m_data;
        occur.free_mem();
    }

    void initialize( int n )
    {
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = nil;
        cur = 0;
    }
    void clean()
    {
        for ( int i = 0 ; i < occur.m_num ; ++i )
        {
            m_data[occur[i]] = nil;
        }
        occur.clean();
        cur = 0;
    }
    
    //init keys 0-n, value as 0
    void init_keys(int n){
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i ){
            m_data[i] = 0;
            occur.push_back( i );
            cur++;
        }
    }
    //reset all values to be zero
    void reset_zero_values(){
        // for ( int i = 0 ; i < m_num ; ++i )
            // m_data[i] = 0.0;
        memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    void reset_one_values(){
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = 1.0;
        // memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    _T get( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //  return -8;
        //}
        return m_data[p];
    }
    _T& operator[](  int p )
    {
        //if ( i < 0 || i >= m_num ) 
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[p];
    }
    void erase( int p )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap get out of range!!!\n");
        //}
        m_data[p] = nil;
        cur--;
    }
    bool notexist( int p )
    {
        return m_data[p] == nil ;
    }
    bool exist( int p )
    {
        return !(m_data[p] == nil);
    }
    void insert( int p , _T d )
    {
        //if ( p < 0 || p >= m_num ) 
        //{
        //  printf("iMap insert out of range!!!\n");
        //}
        if ( m_data[p] == nil )
        {
            occur.push_back( p );
            cur++;
        }
        m_data[p] = d;
    }
    void inc( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p]++;
    }
    void inc( int p , int x )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p] += x;
    }
    void dec( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("dec some unexisted point\n" );
        //}
        m_data[p]--;
    }
    //close range check when release!!!!!!!!!!!!!!!!!!!!    
};
class Graph;
extern double alpha;
class PageRank{
public:
    vector<double> pg_values;
    map<int, map<int,double > > ppr_results;
    iMap<int> ppr_marker;
    double pr_min;
    int n;
    vector<vector<int>> buckets;
    vector<double> min_in_bucket;
    PageRank(){}

    void load_exact_pagerank(const string path){
        std::string pagerank_file_str = path + ".exact.pagerank";
        ////string pagerank_file_str = path;;
        //INFO(pagerank_file_str);
        if(!exists_test(pagerank_file_str)){
            //INFO("No exact PageRank ");
            printf("No exact PageRank");
            exit(0);
        }
        std::ifstream ifs(pagerank_file_str);
        boost::archive::text_iarchive ia(ifs);
        ia>>pg_values;
        n = pg_values.size();
        int counter1 =0;
        int counter2=0;
        int counter3=0;
        for(int i=0; i<n; i++){
            if(pg_values[i]*0.1<=1) 
                counter1++;
            if(pg_values[i]*0.01<=1)
                counter2++;
            if(pg_values[i]*0.005<=1)
                counter3++;
        }
        printf("The ratio of PageRank scores no larger than 1 %f %f %f\n", counter1*1.0/n, counter2*1.0/n, counter3*1.0/n);
        int num_bucket = ceil(log(n*1.0/alpha)/log(2));
        pr_min = n;
        for(int j=0; j<num_bucket; j++){
            vector<int> tl;
            buckets.push_back(tl);
        }
        for(int i=0; i<pg_values.size(); i++){
            double pr = pg_values[i];
            if(pr< pr_min) pr_min = pr;
            int bucket_pos = floor(log(n*1.0/pr)/log(2));
            //if(pr<=0.22) 
            //    INFO(pr, bucket_pos);
            buckets[bucket_pos].push_back(i);
        }
        double base_2 = n;
        for(int i=0; i<buckets.size(); i++){
            double min = base_2;
            for(int j=0; j<buckets[i].size(); j++){
                if(pg_values[buckets[i][j]]<min) min = pg_values[buckets[i][j]];
            }
            min_in_bucket.push_back(min);
            base_2/=2.0;
        }
        printf("num_bucket = %d\n", num_bucket);
    }
    inline bool exists_test(const std::string &name) {
        ifstream f(name.c_str());
        if (f.good()) {
            f.close();
            return true;
        } else {
            f.close();
            return false;
        }
    }
/*
    void save_exact_pagerank(){
        string pagerank_file_str = config.exact_pagerank_folder+FILESEP+config.graph_alias+".exact.pagerank";
       std::ofstream ofs(pagerank_file_str);
       boost::archive::text_oarchive oa(ofs);
        oa<<pg_values;
    }


    void save_exact_ppr(){
        string ppr_file_str = config.exact_pagerank_folder+FILESEP+config.graph_alias+".exact.ppr";
       std::ofstream ofs(ppr_file_str);
       boost::archive::text_oarchive oa(ofs);
        oa<<ppr_results;
    }

    void load_exact_ppr(){
        string ppr_file_str = config.exact_pagerank_folder+FILESEP+config.graph_alias+".exact.ppr";
        INFO(ppr_file_str);
        if(!exists_test(ppr_file_str)){
            INFO("No exact ppr ");
            exit(0);
        }
        std::ifstream ifs(ppr_file_str);
        boost::archive::text_iarchive ia(ifs);
        ia>>ppr_results;
        INFO("finished loading exact ppr");
        ppr_marker.nil = -1;
        ppr_marker.initialize(n);
        for(auto x: ppr_results){
            ppr_marker.insert(x.first, 1);
        }
    }
*/
    void calculate_exact_pagerank(Graph& graph){
        //INFO("Calculating exact PageRank");
        pg_values = vector<double>(graph.n, 1.0/graph.n);
        int iteration = 0;
        //Timer timer(PAGERANK_TIME);
        //int power_iterations = 500;
        int power_iterations = 100;
        while(iteration<power_iterations){
            iteration++;
            //INFO(iteration, config.graph_alias);
            vector<double> previous_pg_value = pg_values;
            double pr_sum =0;
            double dead_end_scores=0;
            for(int i=0; i<graph.n; i++){
                pg_values[i]=graph.alpha/graph.n;
                int in_degree_i = graph.degree[1][i];
                //INFO(i, in_degree_i);
                for(int j=0; j<in_degree_i; j++){
                    int v = graph.get(1, i, j);
                    int out_deg_v = graph.degree[0][v];
                    pg_values[i]+=previous_pg_value[v]*(1-graph.alpha)/out_deg_v;
                }
                if(graph.degree[0][i]==0){
                    pg_values[i]+=previous_pg_value[i]*(1-graph.alpha);
                }
            }
            for(int i=0; i<graph.n; i++){
                pr_sum+=pg_values[i];
            }
            //INFO(pr_sum);
        }

        for(int i=0; i<graph.n; i++){
            pg_values[i]*=graph.n;
        }
    }
/*
    void calculate_exact_personalized_pagerank(const Graph& graph, Config& config){
        iMap<double> reserve;
        iMap<double> residual;
        vector<pair<int,double>> degree_pair;
        for(int i=0; i<graph.n; i++){
            degree_pair.push_back(make_pair(i, graph.degree[1][i]));
        }
        sort(degree_pair.begin(),degree_pair.end(), pair_sorter_large_first);

        reserve.nil = -1;
        residual.nil= -1;
        reserve.initialize(graph.n);
        residual.initialize(graph.n);
        
        long pre_store_size = graph.m;
        long used_size = 0;
        int pre_store_num = max( graph.get_avg_degree(), log(graph.n));
        int count_10000=0;
        int count_5000=0;
        int count_1000=0;
        for(int i=0; i<graph.n; i++){
            if(degree_pair[i].second>10000){
                //INFO(degree_pair[i].second);
                count_10000++;
            }
            if(degree_pair[i].second > 5000){
                count_5000++;
            }
            if(degree_pair[i].second > 1000){
                count_1000++;
            }
        }
        INFO(count_10000, count_5000, count_1000);
        int stored_num = 0;
        for(int i=0;  ; i++){
            if(used_size > pre_store_size) break;
            int target = degree_pair[i].first;
            INFO(degree_pair[i].second);
            residual.clean();
            reserve.clean();
            residual.insert(target, 1);
            vector<int> q;
            vector<int> idx = vector<int>(graph.n, 0);
            q.reserve(graph.n);
            idx[target]=1;
            q.push_back(target);
            double rmax = 1.0/graph.n;
            for(long i=0; i<q.size(); i++){
                int v = q[i];
                idx[v] = false;
                if(residual[v]>= rmax){
                }else
                    break;
                //INFO(residual[v]);
                if(reserve.exist(v))
                    reserve[v]+=config.alpha*residual[v];
                else{
                    reserve.insert(v, config.alpha*residual[v]);
                }
                double remain_residual = (1-config.alpha)*residual[v];
                residual.erase(v);
                int in_degree_v = graph.degree[1][v];
                for(int j=0; j<in_degree_v; j++){
                    int u = graph.get(1, v, j);
                    int out_deg_u = graph.degree[0][u];
                    if(residual.exist(u))
                        residual[u]+=remain_residual/out_deg_u;
                    else
                        residual.insert(u, remain_residual/out_deg_u);
                    if(idx[u]!=true && residual[u]>rmax){
                        idx[u] = true;
                        q.push_back(u);
                    }
                }
            }

            map<int,double> pprs;
            double threshold = 1.0/graph.n;
            for(int i=0; i<reserve.occur.m_num; i++){
                int v = reserve.occur[i];
                double ppr_score = reserve[v];
                if(ppr_score<threshold){
                    //INFO("too small");
                    continue;
                }
                pprs.insert(make_pair(v, ppr_score));
            }
            ppr_results.insert(make_pair(target, pprs));
            INFO(pprs.size());
            used_size+=pprs.size();
            stored_num++;
        }
        INFO(stored_num);
    }
    */
};

PageRank pr;
#endif
