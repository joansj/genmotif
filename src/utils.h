#pragma once
#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <math.h>
using namespace std;

#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include<Eigen/Dense>
using namespace Eigen;

#define IINFINITE   (int)2147483646
#define IPINFINITE  (int)2047483646
#define FINFINITE   numeric_limits<float>::infinity()
#define FPINFINITE  (float)1e30
#define DINFINITE   numeric_limits<double>::infinity()
#define LINFINITE   (long long)9223372036854775800

double extract_time(string runtimetext);		// in seconds

void load_txt(string filename,MatrixXf & stream,int n,int d,bool verb=false,char delimiter=',');
void save_txt(string filename,vector<int> & v);
void save_txt(string filename,vector<float> & v);
void save_txt(string filename,VectorXf & v);
void save_txt(string filename,MatrixXf & x,char delimiter=',');

void moments(MatrixXf & x,VectorXf & mu,VectorXf & sigma);
void znormalize(MatrixXf & x);
void brute_force_upsample(MatrixXf & x,int w);
void smooth_moving_average(MatrixXf & x,int w);
void pointwise_derivative(MatrixXf & x);

int modulo(int a,int b);
float modulo(float a,float b);
void argsort(vector<float> & in,vector<int> & out,bool decreasing_order=false);
void argsort(vector<int> & in,vector<int> & out,bool decreasing_order=false);
int argmin(vector<float> & x);
int argmin(VectorXf & x);
int argmax(vector<float> & x);
int argmax(VectorXf & x);

void print_vector(vector<int> & v);
void print_vector(vector<long> & v);
void print_vector(vector<long long> & v);
void print_vector(vector<float> & v);
void print_vector(vector<double> & v);

#endif
