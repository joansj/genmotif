#pragma once
#ifndef __GENMOTIF_H__
#define __GENMOTIF_H__

#include <chrono>
#include <unordered_map>
using namespace std;

#include "utils.h"
#include "rng.h"
#include "clustering.h"

#undef _WORK_AS_RANDOM_SEARCH_
//#define _WORK_AS_RANDOM_SEARCH_
#undef _FORCE_INITIAL_VALID_SOLUTIONS_
#define _FORCE_INITIAL_VALID_SOLUTIONS_
#undef _ZNORMALIZE_PATTERNS_
#define _ZNORMALIZE_PATTERNS_
#undef _USE_INVALID_MARK_
//#define _USE_INVALID_MARK_                  -10000  // Less than or equal to this number is considered a break mark
#define _PATTERN_DISTANCE_                  _CLUSTERING_SQUAREDEUCLIDEAN_
#define _CLUSTERING_INDEX_                  _CLUSTERING_DAVIESBOULDIN_

class cindividual{
    public:
        vector<int> start,winlen;
        vector<float> cluind;
    public:
        cindividual(){};
        ~cindividual(){};
};

class cgenmotif{
    private:
        bool _verb;
		int _trackfitfreq,_tslen,_tsdim,_wmin,_wmax,_k,_f;
        int _popsize;
        float _prob,_sigma_start,_sigma_len,_sigma_cluster;
        crng _rng;
        cclustering _clu;
    public:
        cgenmotif(){};
        ~cgenmotif(){};
        void init(int n,int dim,int wmin,int wmax,int nmotifs,int fmotifs,bool verbose=false,int track_fitness_freq=-1);
        void load_params(string & fn_config);
        void config(int popsize,float sigma_factor);
        void close();
        void run(MatrixXf & stream,double tmax,string & fn_out);    // time in seconds
    private:
        cindividual new_individual(MatrixXf & stream);
        void crossover(cindividual & x,cindividual & y);
        void mutate(cindividual & x);
        double mutation_random_number();
        void shuffle_patterns(cindividual & x);
        float evaluate_fitness(cindividual & x,MatrixXf & stream);
        bool valid_solution(cindividual & x);
        void write_solution(string & fn_out,vector<pair<double,float> > & fc,cindividual & x,MatrixXf & stream);
};

#endif
