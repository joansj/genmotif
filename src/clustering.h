#pragma once
#ifndef __CLUSTERING_H__
#define __CLUSTERING_H__

#include "utils.h"

#define _CLUSTERING_SQUAREDEUCLIDEAN_	1
#define _CLUSTERING_DTW_				2
#define _CLUSTERING_DTW_BW_				0.025
#define _CLUSTERING_SPECIAL_CALLS_		3
#define _CLUSTERING_LINF_				4
#define _CLUSTERING_TESTDIST_			5
#define _CLUSTERING_SQEUCMED_			6

#define _CLUSTERING_DAVIESBOULDIN_		1
#define _CLUSTERING_KMEANSOBJECTIVE_	2
#define _CLUSTERING_SILHOUETTE_			3
#define _CLUSTERING_MODIFIED_DBI_		4

class cclustering{
	private:
		int _disttype,_cluindex;
		float _dtw_w;
	public:
		cclustering(){};
		~cclustering(){};
		void init(int dist,int index,float dtw_w=_CLUSTERING_DTW_BW_);
		void close();
		void get_assignation(vector<float> & clufit,int minfreq,vector<int> & cluster);
		void get_centers(vector<int> & cluster,vector<MatrixXf> & X,vector<MatrixXf> & centers);
		float get_clustering_index(vector<int> & cluster,vector<MatrixXf> & X);
	private:
		float distance(MatrixXf & x,MatrixXf & y);
		float davis_bouldin_index(vector<int> & cluster,vector<MatrixXf> & X);
		float davis_bouldin_index_modified(vector<int> & cluster,vector<MatrixXf> & X);
		float kmeans_objective(vector<int> & cluster,vector<MatrixXf> & X);
		float silhouette(vector<int> & cluster,vector<MatrixXf> & X);
};

#endif
