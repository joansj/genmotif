#include "utils.h"

double extract_time(string runtimetext){
	double t;
	string num,tim,time_chars("smh");

	size_t found=runtimetext.find_first_of(time_chars);
	if(string::npos==found) return -1.0;
	num=runtimetext.substr(0,found);
	tim=runtimetext.substr(found);
	t=atof(num.c_str());
	if(tim=="m"){
		t*=60;
	}else{
		if(tim=="h"){
			t*=60*60;
		}
	}
	return t;
}

// *************************************************************************************************************************

void load_txt(string filename,MatrixXf & stream,int rows,int cols,bool verb,char sep){
     int i,j;
     int verb_freq_print=1000;
     ifstream myfile;
     string line_in,sub_in;
     string::size_type loc1,loc2;
     vector<float> aux;

     stream.setZero(rows,cols);
     myfile.open(filename.c_str());
     if(myfile.is_open()){
          for(i=0;i<rows;i++){
               getline(myfile,line_in);
               loc1=0;
               for(j=0;j<cols;j++){
                    loc2=line_in.find(sep,loc1);
                    sub_in=line_in.substr(loc1,loc2-loc1);
                    stream(i,j)=strtof(sub_in.c_str(),NULL);
                    loc1=loc2+1;
               }
               if(verb){ if(i%verb_freq_print==0){ cout<<"\rLoading ["<<filename<<"]: "<<i; cout.flush(); } }
          }
          if(verb){ cout<<"\rLoading ["<<filename<<"]: "<<i<<endl; }
     }else{
          cout<<"ERROR: Could not open file ["<<filename<<"]"<<endl;
          exit(1);
     }
     myfile.close();
}

void save_txt(string filename,vector<int> & v){
     int i;
     ofstream myfile;

     myfile.open(filename.c_str());
     if(myfile.is_open()){
          for(i=0;i<(int)v.size();i++) myfile<<v[i]<<endl;
     }else{
          cout<<"ERROR: Could not open file ["<<filename<<"]"<<endl;
          exit(1);
     }
     myfile.close();
}

void save_txt(string filename,vector<float> & v){
     int i;
     ofstream myfile;

     myfile.open(filename.c_str());
     if(myfile.is_open()){
          for(i=0;i<(int)v.size();i++) myfile<<v[i]<<endl;
     }else{
          cout<<"ERROR: Could not open file ["<<filename<<"]"<<endl;
          exit(1);
     }
     myfile.close();
}

void save_txt(string filename,VectorXf & v){
     int i;
     ofstream myfile;

     myfile.open(filename.c_str());
     if(myfile.is_open()){
          for(i=0;i<v.size();i++) myfile<<v(i)<<endl;
     }else{
          cout<<"ERROR: Could not open file ["<<filename<<"]"<<endl;
          exit(1);
     }
     myfile.close();
}

void save_txt(string filename,MatrixXf & x,char delimiter){
     int i,j;
     ofstream myfile;

     myfile.open(filename.c_str());
     if(myfile.is_open()){
          for(i=0;i<x.rows();i++){
               for(j=0;j<x.cols();j++){
                    myfile<<x(i,j);
                    if(j<x.cols()-1) myfile<<delimiter;
               }
               myfile<<endl;
          }
     }else{
          cout<<"ERROR: Could not open file ["<<filename<<"]"<<endl;
          exit(1);
     }
     myfile.close();
}

// *************************************************************************************************************************

void moments(MatrixXf & x,VectorXf & mu,VectorXf & sigma){
     int i,j;
     float delta;

     mu.setZero(x.cols());
     sigma.setZero(x.cols());
     for(j=0;j<x.cols();j++){
          for(i=0;i<x.rows();i++){
               delta=x(i,j)-mu(j);
               mu(j)+=delta/((float)(i+1));
               sigma(j)+=delta*(x(i,j)-mu(j));
          }
          sigma(j)=sqrt(sigma(j)/((float)(x.rows()-1)));
     }
}

void znormalize(MatrixXf & x){
     int j;
     VectorXf mu,sigma;

     moments(x,mu,sigma);
     for(j=0;j<x.cols();j++){
          x.col(j)=x.col(j).array()-mu(j);
          if(sigma(j)>0) x.col(j)=x.col(j).array()/sigma(j);
     }
}

void brute_force_upsample(MatrixXf & x,int w){
     int i;
     float factor=((float)x.rows())/((float)w);
     MatrixXf aux=x;

     x.setZero(w,x.cols());
     for(i=0;i<w;i++){
          x.row(i)=aux.row((int)floor(i*factor));
     }
}

void smooth_moving_average(MatrixXf & x,int winlen){
	int i,j,w,ww;
	MatrixXf y;

	w=floor(0.5*winlen);
	if(w>1){
		y=x;
		for(i=0;i<y.rows();i++){
			ww=min(w,min(i,((int)y.rows())-i-1));
			x.row(i)=VectorXf::Zero(y.cols());
			for(j=i-ww;j<=i+ww;j++){
				x.row(i)+=y.row(j);
			}
			x.row(i)/=(float)(2*ww+1);
		}
	}
}

void pointwise_derivative(MatrixXf & x){
	int i;

	for(i=0;i<x.rows()-1;i++) x.row(i)=x.row(i+1)-x.row(i);
	x.row(x.rows()-1)=VectorXf::Zero(x.cols());
}

// *************************************************************************************************************************

int modulo(int a,int b){
     return (a%b+b)%b;
}

float modulo(float a,float b){
     return fmod(fmod(a,b)+b,b);
}

void argsort(vector<float> & in,vector<int> & out,bool decreasing_order){
     int i;
     vector<pair<float,int> > vp(in.size());
     for(i=0;i<(int)in.size();i++) vp[i]=make_pair(in[i],i);
     sort(vp.begin(),vp.end());
     out.resize(in.size());
     for(i=0;i<(int)in.size();i++) out[i]=vp[i].second;
     if(decreasing_order) reverse(out.begin(),out.end());
}

void argsort(vector<int> & in,vector<int> & out,bool decreasing_order){
     int i;
     vector<pair<int,int> > vp(in.size());
     for(i=0;i<(int)in.size();i++) vp[i]=make_pair(in[i],i);
     sort(vp.begin(),vp.end());
     out.resize(in.size());
     for(i=0;i<(int)in.size();i++) out[i]=vp[i].second;
     if(decreasing_order) reverse(out.begin(),out.end());
}

int argmin(vector<float> & x){
	int i,id;
	float mn=FINFINITE;

	for(i=0;i<x.size();i++){
		if(x[i]<mn){
			mn=x[i];
			id=i;
		}
	}
	return id;
}

int argmin(VectorXf & x){
	int i,id;
	float mn=FINFINITE;

	for(i=0;i<x.size();i++){
		if(x(i)<mn){
			mn=x(i);
			id=i;
		}
	}
	return id;
}

int argmax(vector<float> & x){
	int i,id;
	float mx=-FINFINITE;

	for(i=0;i<(int)x.size();i++){
		if(x[i]>mx){
			mx=x[i];
			id=i;
		}
	}
	return id;
}

int argmax(VectorXf & x){
	int i,id;
	float mx=-FINFINITE;

	for(i=0;i<x.size();i++){
		if(x(i)>mx){
			mx=x(i);
			id=i;
		}
	}
	return id;
}

// *************************************************************************************************************************

void print_vector(vector<int> & v){
     cout<<"[ ";
     for(unsigned int i=0;i<v.size();i++){
          cout<<v[i]<<" ";
     }
     cout<<"]"<<endl;
}

void print_vector(vector<long> & v){
     cout<<"[ ";
     for(unsigned int i=0;i<v.size();i++){
          cout<<v[i]<<" ";
     }
     cout<<"]"<<endl;
}

void print_vector(vector<long long> & v){
     cout<<"[ ";
     for(unsigned int i=0;i<v.size();i++){
          cout<<v[i]<<" ";
     }
     cout<<"]"<<endl;
}

void print_vector(vector<float> & v){
     cout<<"[ ";
     for(unsigned int i=0;i<v.size();i++){
          cout<<v[i]<<" ";
     }
     cout<<"]"<<endl;
}

void print_vector(vector<double> & v){
     cout<<"[ ";
     for(unsigned int i=0;i<v.size();i++){
          cout<<v[i]<<" ";
     }
     cout<<"]"<<endl;
}

// *************************************************************************************************************************
