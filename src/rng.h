#pragma once
#ifndef __RNG_H__
#define __RNG_H__

#include <algorithm>
#include <random>
#include <math.h>
using namespace std;

#include "utils.h"

class crng{
     private:
          mt19937 gen;
          uniform_real_distribution<double> U;
          normal_distribution<double> N;
          double levy_betainv,levy_sigma;
     public:
          crng(){};
          ~crng(){};
          void init(int seed=-1,unsigned int warm_ups=100000,double beta_levy=1.5);
          void close();
          double uniform();
          void uniform(int size,VectorXf & v);
          void uniform(int size,VectorXd & v);
          void uniform(int rows,int cols,MatrixXd & m);
          double normal();
          void normal(int size,VectorXf & v);
          void normal(int size,VectorXd & v);
          void normal(int rows,int cols,MatrixXd & m);
          double levy();
          double levy1();
          void shuffle(VectorXi & x);
};


#endif
