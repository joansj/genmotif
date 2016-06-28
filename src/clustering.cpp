#include "clustering.h"

void cclustering::init(int dist,int index,float dtw_w){
    _disttype=dist;
    _cluindex=index;
    _dtw_w=dtw_w;
}

void cclustering::close(){
}

void cclustering::get_assignation(vector<float> & clufit,int freq,vector<int> & cluster){
     int i,j,k,n=(int)clufit.size();
     vector<int> o;

     // Get assignation
     argsort(clufit,o);
     cluster.resize(n);
     for(i=0;i<n;i++) cluster[i]=-1;
     k=0;
     j=0;
     for(i=0;i<n;i++){
          cluster[o[i]]=k;
          j++;
          if(j>=freq){
               k++;
               j=0;
          }
     }
}

void cclustering::get_centers(vector<int> & cluster,vector<MatrixXf> & X,vector<MatrixXf> & centers){
     int i,j,k,numclusters,cnt,n=(int)cluster.size(),w=X[0].rows(),d=X[0].cols();
     vector<int> cand;
     MatrixXf average,dist;
     VectorXf sum;

     // Get number of clusters
     numclusters=-1;
     for(i=0;i<n;i++) if(cluster[i]>numclusters) numclusters=cluster[i];
     numclusters++;

     // Centers
     centers.clear();
     switch(_disttype){
         case _CLUSTERING_TESTDIST_:
         case _CLUSTERING_LINF_:
         case _CLUSTERING_SQUAREDEUCLIDEAN_:
             // Compute averages
             for(i=0;i<numclusters;i++){
                  average.setZero(w,d);
                  cnt=0;
                  for(j=0;j<n;j++){
                       if(cluster[j]!=i) continue;
                       average+=X[j];
                       cnt++;
                  }
                  average/=(float)cnt;
                  centers.push_back(average);
             }
             break;
        default:
            // Get medoids
            for(i=0;i<numclusters;i++){
                cand.clear();
                for(j=0;j<n;j++) if(cluster[j]==i) cand.push_back(j);
                dist=MatrixXf::Zero(cand.size(),cand.size());
                for(j=0;j<dist.rows();j++){
                    for(k=j+1;k<dist.cols();k++){
                        dist(j,k)=this->distance(X[cand[j]],X[cand[k]]);
                        dist(k,j)=dist(j,k);
                    }
                }
                sum=dist.rowwise().sum();
                j=argmin(sum);
                centers.push_back(X[cand[j]]);
            }
            break;
    }
}

float cclustering::get_clustering_index(vector<int> & cluster,vector<MatrixXf> & X){
    switch(_cluindex){
        case _CLUSTERING_DAVIESBOULDIN_:
            return this->davis_bouldin_index(cluster,X);
            break;
        case _CLUSTERING_MODIFIED_DBI_:
            return this->davis_bouldin_index_modified(cluster,X);
            break;
        case _CLUSTERING_KMEANSOBJECTIVE_:
            return this->kmeans_objective(cluster,X);
            break;
        default:
            return 1-this->silhouette(cluster,X);
            break;
    }
}

float cclustering::davis_bouldin_index(vector<int> & cluster,vector<MatrixXf> & X){
    int i,j;
    float dbi,aux;
    vector<MatrixXf> centers;

    // Get centers
    this->get_centers(cluster,X,centers);
    int numclusters=(int)centers.size();
    if(numclusters<2) return FINFINITE;

    // Calculate deviance for each cluster
    VectorXf numelem=VectorXf::Zero(numclusters);
    VectorXf deviance=VectorXf::Zero(numclusters);
    for(i=0;i<(int)cluster.size();i++){
        numelem(cluster[i])+=1.0f;
        deviance(cluster[i])+=this->distance(X[i],centers[cluster[i]]);
    }
    deviance.array()/=numelem.array();

    // DB index
    VectorXf mx=-FINFINITE*VectorXf::Ones(numclusters);
    for(i=0;i<numclusters;i++){
        for(j=i+1;j<numclusters;j++){
            aux=(deviance(i)+deviance(j))/this->distance(centers[i],centers[j]);
            if(aux>mx(i)) mx(i)=aux;
            if(aux>mx(j)) mx(j)=aux;
        }
    }
    return mx.sum()/((float)mx.size());
}

float cclustering::davis_bouldin_index_modified(vector<int> & cluster,vector<MatrixXf> & X){
    int i,j,k,w;
    float dbi,aux,dist;
    MatrixXf tmp;
    vector<MatrixXf> C,Chalf,Xhalf;

    w=floor(X[0].rows()*0.5);
    for(i=0;i<(int)X.size();i++){
        tmp=X[i].block(0,0,w,X[i].cols());
        znormalize(tmp);
        Xhalf.push_back(tmp);
    }

    // Get centers
    this->get_centers(cluster,X,C);
    int numclusters=(int)C.size();
    if(numclusters<2) return FINFINITE;
    this->get_centers(cluster,Xhalf,Chalf);

    // Calculate deviance for each cluster
    VectorXf numelem=VectorXf::Zero(numclusters);
    VectorXf deviance=VectorXf::Zero(numclusters);
    for(i=0;i<(int)cluster.size();i++){
        numelem(cluster[i])+=1.0f;
        deviance(cluster[i])+=this->distance(Xhalf[i],Chalf[cluster[i]]);
    }
    deviance.array()/=numelem.array();

    // DB index --- Modified!
    VectorXf mx=-FINFINITE*VectorXf::Ones(numclusters);
    for(i=0;i<numclusters;i++){
        for(j=i+1;j<numclusters;j++){
            dist=FINFINITE;
            for(k=0;k<w;k++){
                tmp=C[j].block(k,0,w,C[j].cols());
                znormalize(tmp);
                aux=this->distance(Chalf[i],tmp);
                if(aux<dist) dist=aux;
                tmp=C[i].block(k,0,w,C[i].cols());
                znormalize(tmp);
                aux=this->distance(Chalf[j],tmp);
                if(aux<dist) dist=aux;
            }
            aux=(deviance(i)+deviance(j))/pow(dist,0.5f);
            if(aux>mx(i)) mx(i)=aux;
            if(aux>mx(j)) mx(j)=aux;
        }
    }
    return mx.sum()/((float)mx.size());
}


float cclustering::kmeans_objective(vector<int> & cluster,vector<MatrixXf> & X){
    int i;
    float obj;
    vector<MatrixXf> centers;

    this->get_centers(cluster,X,centers);
    int numclusters=(int)centers.size();
    if(numclusters<1) return FINFINITE;

    obj=0.0f;
    for(i=0;i<(int)cluster.size();i++) obj+=this->distance(X[i],centers[cluster[i]]);

    return obj;
}

float cclustering::silhouette(vector<int> & cluster,vector<MatrixXf> & X){
    int i,j,n=(int)cluster.size(),numclusters;
    float silhouette,ai,bi;
    VectorXf avgdist,avgnum;

	// Get number of clusters
    numclusters=-1;
	for(i=0;i<n;i++) if(cluster[i]>numclusters) numclusters=cluster[i];
	numclusters++;

	// Compute silhouette
    silhouette=0.0f;
    for(i=0;i<n;i++){
        avgdist=VectorXf::Zero(numclusters);
        avgnum=VectorXf::Zero(numclusters);
        for(j=0;j<n;j++){
            if(i==j) continue;
            avgdist(cluster[j])+=this->distance(X[i],X[j]);
            avgnum(cluster[j])+=1.0f;
        }
        ai=avgdist(cluster[i])/avgnum(cluster[i]);
        avgdist(cluster[i])=FINFINITE;
        bi=(avgdist.array()/avgnum.array()).minCoeff();
        silhouette+=(bi-ai)/max(ai,bi);
    }
    return silhouette/(float)n;
}

float cclustering::distance(MatrixXf & x,MatrixXf & y){
    int i,j,w,ii,im;
    float dist_xy,dist_yx,mindist,aux;
    vector<float> cx,cy;
    VectorXf tmp;
    MatrixXf D;

    switch(_disttype){
        case _CLUSTERING_LINF_:
            return (x-y).array().abs().maxCoeff();
            break;
        case _CLUSTERING_SQEUCMED_:
        case _CLUSTERING_SQUAREDEUCLIDEAN_:
            return (x-y).array().square().sum();
            break;
        case _CLUSTERING_DTW_:
            w=floor(_dtw_w*x.rows());
            // Correct, linear-mem. working version!
            D=FINFINITE*MatrixXf::Ones(2,y.rows());
            D(0,0)=(x.row(0)-y.row(0)).array().square().sum();
            for(j=1;j<=w;j++) D(0,j)=D(0,j-1)+(x.row(0)-y.row(j)).array().square().sum();
            for(i=1;i<x.rows();i++){
                ii=i%2;
                im=(i-1)%2;
                if(i<=w){
                    D(ii,0)=D(im,0)+(x.row(i)-y.row(0)).array().square().sum();
                }else{
                    D(ii,i-w-1)=FINFINITE;
                }
                for(j=max(1,i-w);j<min((int)y.rows(),i+w+1);j++){
                    D(ii,j)=min(D(im,j-1),min(D(im,j),D(ii,j-1)))+(x.row(i)-y.row(j)).array().square().sum();
                }
            }
            return D((x.rows()-1)%2,y.rows()-1);
            break;
        case _CLUSTERING_SPECIAL_CALLS_:
            cx.push_back(x(0,0));
            for(i=1;i<x.rows();i++) cx.push_back(x(i,0)+cx[i-1]);
            cy.push_back(y(0,0));
            for(j=1;j<y.rows();j++) cy.push_back(y(j,0)+cy[j-1]);
            dist_xy=0.0f;
            for(i=0;i<x.rows();i++){
                mindist=FINFINITE;
                for(j=0;j<y.rows();j++){
                    aux=abs(cx[i]-cy[j]);
                    //aux+=1.0f; aux*=1+sqrt(abs(x(i,1)-y(j,1)));
                    if(aux<mindist) mindist=aux;
                }
                dist_xy+=mindist;
            }
            dist_xy/=(float)x.rows();
            dist_yx=0.0f;
            for(j=0;j<y.rows();j++){
                mindist=FINFINITE;
                for(i=0;i<x.rows();i++){
                    aux=abs(cx[i]-cy[j]);
                    //aux+=1.0f; aux*=1+sqrt(abs(x(i,1)-y(j,1)));
                    if(aux<mindist) mindist=aux;
                }
                dist_yx+=mindist;
            }
            dist_yx/=(float)y.rows();
            return dist_xy+dist_yx;
            break;
        case _CLUSTERING_TESTDIST_:
            tmp=x.col(0);
            i=argmin(tmp);
            tmp=y.col(0);
            j=argmin(tmp);
            return (x-y).array().square().sum()*pow(1+i+j,0.3333);
            break;
        default:
            return FINFINITE;
            break;
    }
}
