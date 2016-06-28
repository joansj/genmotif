#include "main.h"

int main(int argc,char *argv[]){
    bool verbose=false;
    int i,n,dim,wmin,wmax,k,freq,aux;
	double tmax;
    string fn_stream,runtimetext,fn_out,fn_conf;

    // Load input parameters
    VectorXi inputs=VectorXi::Zero(6);
    for(i=1;i<argc;i++){
        if(strcmp(argv[i],"-s")==0){
            fn_stream.assign(argv[i+1]);
            inputs(0)=1;
        }
        if(strcmp(argv[i],"-l")==0){
            n=atoi(argv[i+1]);
            dim=atoi(argv[i+2]);
            inputs(1)=1;
        }
        if(strcmp(argv[i],"-w")==0){
            wmin=atoi(argv[i+1]);
            wmax=atoi(argv[i+2]);
            inputs(2)=1;
        }
        if(strcmp(argv[i],"-m")==0){
            k=atoi(argv[i+1]);
        	freq=atoi(argv[i+2]);
            inputs(3)=1;
        }
        if(strcmp(argv[i],"-t")==0){
            runtimetext.assign(argv[i+1]);
            tmax=extract_time(runtimetext);
        	if(tmax<=0) tmax=DINFINITE;
            inputs(4)=1;
        }
        if(strcmp(argv[i],"-o")==0){
            fn_out.assign(argv[i+1]);
            inputs(5)=1;
        }
        if(strcmp(argv[i],"-c")==0){
            fn_conf.assign(argv[i+1]);
        }
        if(strcmp(argv[i],"-v")==0){
            aux=atoi(argv[i+1]);
            if(aux>0) verbose=true;
        }
    }
    if(inputs.sum()!=6){
        cout<<"ERROR: Wrong call. Use [optional]:"<<endl;
        cout<<argv[0]<<" -s <stream_filename.txt> -l <stream_len> <stream_dim> -w <wmin> <wmax> -m <num_motifs> <frequency> -t <running_time(000smh)> -o <output_filename(no_extension)> [-c <config_filename> -v <verbose(0/1)>]"<<endl;
        return 1;
    }

    // Load stream time series
    MatrixXf stream;
    load_txt(fn_stream,stream,n,dim,verbose);

    // Run motif discovery
    cgenmotif motif_finder;
    motif_finder.init(n,dim,wmin,wmax,k,freq,verbose,_TRACKFITS_);
    if(!fn_conf.empty()) motif_finder.load_params(fn_conf);
    motif_finder.run(stream,tmax,fn_out);
    motif_finder.close();

    return 0;
}
