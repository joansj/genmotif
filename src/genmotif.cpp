#include "genmotif.h"

// =========================================================================================================

void cgenmotif::init(int n,int dim,int wmin,int wmax,int nmotifs,int fmotifs,bool verbose,int trackfitfreq){
    // Assign values
    _tslen=n;
    _tsdim=dim;
    _wmin=wmin;
    _wmax=wmax;
    _k=nmotifs;
    _f=fmotifs;
    _verb=verbose;
	_trackfitfreq=trackfitfreq;
    if(_verb) cout<<"Task specifications: n=["<<_tslen<<","<<_tsdim<<"], w=["<<_wmin<<","<<_wmax<<"], m=["<<_k<<","<<_f<<"]"<<endl;
    if(_CLUSTERING_INDEX_==_CLUSTERING_MODIFIED_DBI_){
        _wmin*=2;
        _wmax*=2;
    }

    // Default parameters
    this->config(15,0.01);

    // Inits
    _prob=1.0f/((float)(_k*_f));
    _rng.init();
    _clu.init(_PATTERN_DISTANCE_,_CLUSTERING_INDEX_);
}

void cgenmotif::close(){
    // Closes
    _clu.close();
    _rng.close();
}

// ------------------------------------------------------------------------------------------------------------

void cgenmotif::load_params(string & fn_conf){
    int popsize;
    float pc_elite,prob_f,sigma,extrapar;
    string foo;
    ifstream myfile;

    myfile.open(fn_conf.c_str());
    if(myfile.is_open()){
        myfile>>foo>>popsize;
        myfile>>foo>>sigma;
        myfile.close();
    }else{
        cout<<"ERROR: Could not open configuration file ["<<fn_conf<<"]"<<endl;
		exit(1);
    }
    this->config(popsize,sigma);
}

void cgenmotif::config(int popsize,float sigma_factor){
    _popsize=1+2*floor(popsize/2.0f);
    _sigma_start=sigma_factor*_tslen;
    _sigma_len=sigma_factor*(_wmax+1-_wmin);
    _sigma_cluster=sigma_factor;
}

// =========================================================================================================

void cgenmotif::run(MatrixXf & stream,double t_max,string & fn_out){
    int i;
    long long it;
    float best_fit;
    chrono::time_point<std::chrono::system_clock> t_start;
    chrono::duration<double> t_elapsed;
    vector<cindividual> population,new_population;
    vector<float> fitness;
    vector<pair<double,float> > trackfitcurve;

    // Print configuration
    if(_verb){
        cout<<"Algorithm parameters: popsize="<<_popsize<<", sigma="<<_sigma_cluster<<endl;
        cout<<"Algorithm definitions: initial_valid=";
#ifdef _FORCE_INITIAL_VALID_SOLUTIONS_
        cout<<"true";
#else
        cout<<"false";
#endif
        cout<<", znormalize=";
#ifdef _ZNORMALIZE_PATTERNS_
        cout<<"true";
#else
        cout<<"false";
#endif
        cout<<", use_mark=";
#ifdef _USE_INVALID_MARK_
        cout<<"true ("<<_USE_INVALID_MARK_<<")";
#else
        cout<<"false";
#endif
        cout<<", distance="<<_PATTERN_DISTANCE_<<", clustering_index="<<_CLUSTERING_INDEX_<<endl;
    }

    // Start counting time
    t_start=chrono::system_clock::now();

    // Init population
    for(i=0;i<_popsize;i++){
        population.push_back(this->new_individual(stream));
        fitness.push_back(this->evaluate_fitness(population[i],stream));
    }
    best_fit=FINFINITE;
    new_population.resize(_popsize);

    // Init counters
    it=0;
    t_elapsed=chrono::system_clock::now()-t_start;

    // Loop
    while(t_elapsed.count()<t_max){
        if(_verb) cout<<"\rRun:\tt="<<fixed<<setprecision(1)<<t_elapsed.count()/60.0<<" min\tit="<<it<<flush;

        // Identify best + Fitness tracking
        i=argmin(fitness);
        if(fitness[i]<best_fit){
            best_fit=fitness[i];
            trackfitcurve.push_back(make_pair(t_elapsed.count(),best_fit));
            if(_verb) cout<<"\t[F="<<fixed<<setprecision(6)<<best_fit<<"]"<<endl;
            // Save intermediate solution?
			if((_trackfitfreq>=0)&&(it%_trackfitfreq==0)) this->write_solution(fn_out,trackfitcurve,population[i],stream);
        }

        // Directly transfer best individual
        new_population[0]=population[i];
        fitness[0]=fitness[i];

        //population[argmax(fitness)]=this->new_individual(stream);

#ifdef _WORK_AS_RANDOM_SEARCH_
        // Just perform random search
        for(i=1;i<_popsize;i++){
            new_population[i]=this->new_individual(stream);
            fitness[i]=this->evaluate_fitness(new_population[i],stream);
        }
#else

        // Generate the rest of the new population + evaluate its fitness
        for(i=1;i<_popsize;i+=2){
            new_population[i]=population[floor(_popsize*_rng.uniform())];
            new_population[i+1]=population[floor(_popsize*_rng.uniform())];
            this->crossover(new_population[i],new_population[i+1]);
            this->mutate(new_population[i]);
            this->mutate(new_population[i+1]);
            fitness[i]=this->evaluate_fitness(new_population[i],stream);
            fitness[i+1]=this->evaluate_fitness(new_population[i+1],stream);
        }

        // Shuffle patterns of each individual
        for(i=0;i<_popsize;i++) this->shuffle_patterns(new_population[i]);

#endif

        // Substitute old population by new one
        population=new_population;

        // Update counters
        it++;
        t_elapsed=chrono::system_clock::now()-t_start;
    }

	// Identify best + Fitness tracking
    i=argmin(fitness);
    if(_verb) cout<<"\rRun:\tt="<<fixed<<setprecision(1)<<t_max/60.0<<" min\tit="<<it<<"\t[F="<<fixed<<setprecision(6)<<fitness[i]<<"]"<<endl;
	trackfitcurve.push_back(make_pair(t_max,fitness[i]));

    // Save best solution
	this->write_solution(fn_out,trackfitcurve,population[i],stream);
	if(_verb) cout<<"End: solutions written ["<<fn_out<<".sol.txt, "<<fn_out<<".fit.txt]"<<endl;

}

// =========================================================================================================

cindividual cgenmotif::new_individual(MatrixXf & stream){
    int i;
    cindividual x;

    // Uniform sampling of the search space
    for(i=0;i<_k*_f;i++){
        x.cluind.push_back(_rng.uniform());
        x.winlen.push_back(_wmin+floor((_wmax+1-_wmin)*_rng.uniform()));
        x.start.push_back(floor((_tslen-x.winlen[i])*_rng.uniform()));
#ifdef _FORCE_INITIAL_VALID_SOLUTIONS_
        // Force that the initial solution is a valid one
        bool valid=false;
        while(!valid){
            x.start[i]=floor((_tslen-x.winlen[i])*_rng.uniform());
            valid=true;
            for(int j=0;j<i;j++){
                if((x.start[i]+x.winlen[i]<=x.start[j])||(x.start[j]+x.winlen[j]<=x.start[i])){
                    // No overlap
                }else{
                    // There is overlap
                    valid=false;
                    break;
                }
            }
#ifdef _USE_INVALID_MARK_
            if(!valid) continue;
            for(int j=x.start[i];j<x.start[i]+x.winlen[i];j++){
                if(stream(j,0)<=_USE_INVALID_MARK_){
                    valid=false;
                    break;
                }
            }
#endif
        }
#endif
    }
    return x;
}

void cgenmotif::crossover(cindividual & x,cindividual & y){
    int i;

    // Swap one position with _prob probability
    for(i=0;i<(int)x.start.size();i++){
         if(_rng.uniform()<_prob){
              swap(x.start[i],y.start[i]);
              swap(x.winlen[i],y.winlen[i]);
              swap(x.cluind[i],y.cluind[i]);
         }
    }
}

void cgenmotif::mutate(cindividual & x){
    int i;

    // Mutate one position with _prob probability
    for(i=0;i<(int)x.start.size();i++){
         if(_rng.uniform()<_prob){
             x.cluind[i]+=_sigma_cluster*this->mutation_random_number();
             x.cluind[i]=modulo(x.cluind[i],1.0f);
             x.winlen[i]+=round(_sigma_len*this->mutation_random_number());
             x.winlen[i]=modulo(x.winlen[i]-_wmin,_wmax+1-_wmin)+_wmin;
             x.start[i]+=round(_sigma_start*this->mutation_random_number());
             x.start[i]=modulo(x.start[i],_tslen-x.winlen[i]);
         }
    }
}

double cgenmotif::mutation_random_number(){
    // Get random number from Cauchy distribution (Levy flight with scale parameter beta=1)
    return _rng.levy1();
}

void cgenmotif::shuffle_patterns(cindividual & x){
     int i;
     VectorXi o;
     cindividual aux=x;

     // Shuffle positions
     o.resize(x.start.size());
     for(i=0;i<o.size();i++) o(i)=i;
     _rng.shuffle(o);
     for(i=0;i<o.size();i++){
         x.start[i]=aux.start[o(i)];
         x.winlen[i]=aux.winlen[o(i)];
         x.cluind[i]=aux.cluind[o(i)];
     }
}

// ------------------------------------------------------------------------------------------------------------

float cgenmotif::evaluate_fitness(cindividual & x,MatrixXf & stream){
	int i;
	MatrixXf tmp;
    vector<int> cluster;
	vector<MatrixXf> patterns;

	// Check validity
	if(!this->valid_solution(x)) return FINFINITE;

	// Resample all subsequences to same length
    for(i=0;i<(int)x.start.size();i++){
#ifdef _USE_INVALID_MARK_
        for(int j=x.start[i];j<x.start[i]+x.winlen[i];j++){
            if(stream(j,0)<=_USE_INVALID_MARK_){
                return FINFINITE;
            }
        }
#endif
        tmp=stream.block(x.start[i],0,x.winlen[i],stream.cols());
#ifdef _ZNORMALIZE_PATTERNS_
        znormalize(tmp);
#endif
        brute_force_upsample(tmp,_wmax);
        patterns.push_back(tmp);
    }

    // Get cluster assignation
    _clu.get_assignation(x.cluind,_f,cluster);

	// Calculate clustering index
    return _clu.get_clustering_index(cluster,patterns);
}

bool cgenmotif::valid_solution(cindividual & x){
	int i,j;
    vector<int> o;

    // Check if some overlap exists (the smart way)
    argsort(x.start,o);
    for(i=0;i<(int)o.size();i++){
        for(j=1;j<((int)o.size())-i;j++){
            if(x.start[o[i+j]]-x.start[o[i]]>=_wmax) break;
            if(x.start[o[i]]+x.winlen[o[i]]>=x.start[o[i+j]]) return false;
        }
    }
	return true;
}

// ------------------------------------------------------------------------------------------------------------

void cgenmotif::write_solution(string & fn_out,vector<pair<double,float> > & fc,cindividual & x,MatrixXf & stream){
	int i,w;
	string fn;
	vector<int> cluster,o;
	ofstream myfile;

    // Write fitness curve
	fn.assign(fn_out);
	fn.append(".fit.txt");
	myfile.open(fn.c_str());
	if(myfile.is_open()){
		for(i=0;i<(int)fc.size();i++){
			myfile<<fc[i].first<<"\t"<<fc[i].second<<endl;
		}
		myfile.close();
	}else{
		cout<<"ERROR: Could not open file ["<<fn<<"]"<<endl;
		exit(1);
	}

    // Get cluster assignation (same as before)
	_clu.get_assignation(x.cluind,_f,cluster);
    argsort(cluster,o);

    // Write solution
	fn.assign(fn_out);
	fn.append(".sol.txt");
	myfile.open(fn.c_str());
	if(myfile.is_open()){
		for(i=0;i<(int)x.start.size();i++){
            w=x.winlen[o[i]];
            if(_CLUSTERING_INDEX_==_CLUSTERING_MODIFIED_DBI_) w/=2;
			myfile<<x.start[o[i]]<<"\t"<<w<<"\t"<<cluster[o[i]]<<endl;
        }
		myfile.close();
	}else{
		cout<<"ERROR: Could not open file ["<<fn<<"]"<<endl;
		exit(1);
	}
}

// =========================================================================================================
