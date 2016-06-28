#include "rng.h"

void crng::init(int seed,unsigned int warm_ups,double beta){
     if(seed<0){
          random_device rd;
          gen.seed(rd());
     }else{
          gen.seed(seed);
     }
     uniform_real_distribution<double> aux1(0.0,1.0);
     U.param(aux1.param());
     normal_distribution<double> aux2(0.0,1.0);
     N.param(aux2.param());
     for(unsigned int i=0;i<warm_ups;i++){
          this->uniform();
          this->normal();
     }
     levy_betainv=1.0/beta;
     levy_sigma=pow((tgamma(1+beta)*sin(M_PI*beta/2.0))/(tgamma((1+beta)/2.0)*beta*pow(2.0,(beta-1)/2.0)),levy_betainv);
}

void crng::close(){
}

double crng::uniform(){
     return U(gen);
}

void crng::uniform(int size,VectorXf & v){
     v.setZero(size);
     for(int i=0;i<size;i++) v(i)=(float)U(gen);
}

void crng::uniform(int size,VectorXd & v){
     v.setZero(size);
     for(int i=0;i<size;i++) v(i)=U(gen);
}

void crng::uniform(int rows,int cols,MatrixXd & m){
     m.setZero(rows,cols);
     for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) m(i,j)=U(gen);
}

double crng::normal(){
     return N(gen);
}

void crng::normal(int size,VectorXf & v){
     v.setZero(size);
     for(int i=0;i<size;i++) v(i)=(float)N(gen);
}

void crng::normal(int size,VectorXd & v){
     v.setZero(size);
     for(int i=0;i<size;i++) v(i)=N(gen);
}

void crng::normal(int rows,int cols,MatrixXd & m){
     m.setZero(rows,cols);
     for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) m(i,j)=N(gen);
}

double crng::levy(){
    return levy_sigma*this->normal()/pow(abs(this->normal()),levy_betainv);
}

double crng::levy1(){
    return this->normal()/abs(this->normal());
}

void crng::shuffle(VectorXi & x){
     int i=((int)x.size());
     while(i>1){
          swap(x(i-1),x(floor(i*this->uniform())));
          i--;
     }
}
