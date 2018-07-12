#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chisq.h"

chisquared::chisquared(){
    death_knell("meaningless constructor");
};

chisquared::chisquared(int id){
    _ncenters=1;
    _dim=id;
    _characteristic_width=1.0;
    _dice=NULL;
    _seed=22;
    _chisq_initialized=0;
    _initialize();
    initialize();
};

chisquared::chisquared(int id, int ic){
    _dim=id;
    _ncenters=ic;
    _characteristic_width=1.0;
    _dice=NULL;
    _seed=22;
    _chisq_initialized=0;
    _initialize();
    initialize();
};


chisquared::chisquared(int id, int ic, double ww){
    _dim=id;
    _ncenters=ic;
    _characteristic_width=ww;
    _dice=NULL;
    _seed=22;
    _chisq_initialized=0;
    _initialize();
    initialize();
}

void chisquared::_initialize(){
    _with_logging=0;
    _time_spent=0.0;
    _called=0;
    _maxs.set_name("chisq_maxs");
    _mins.set_name("chisq_mins");
    _bases.set_name("chisq_bases");
    _widths.set_name("chisq_widths");
    _centers.set_name("chisq_centers");

    int i;
    for(i=0;i<_dim;i++){
        _mins.set(i,-2.0*exception_value);
        _maxs.set(i,2.0*exception_value);
    }
}

void chisquared::initialize(){
    printf("calling chisq initialize ww %e\n",_characteristic_width);

    _dice=new Ran(_seed);

    array_1d<double> trial;
    trial.set_name("chisq_initialize_trial");

    int i,j;
    double mu,norm;
    while(_bases.get_rows()<_dim){
        for(i=0;i<_dim;i++){
           trial.set(i,normal_deviate(_dice,0.0,1.0));
        }
        trial.normalize();
        for(i=0;i<_bases.get_rows();i++){
            mu=0.0;
            for(j=0;j<_dim;j++){
                mu+=trial.get_data(j)*_bases.get_data(i,j);
            }
            for(j=0;j<_dim;j++){
                trial.subtract_val(j,mu*_bases.get_data(i,j));
            }
        }
        norm=trial.normalize();
        if(norm>1.0e-20){
            _bases.add_row(trial);
        }
    }

    int k;

    for(i=0;i<_dim;i++){
        for(j=i;j<_dim;j++){
            mu=0.0;
            for(k=0;k<_dim;k++){
                mu+=_bases.get_data(i,k)*_bases.get_data(j,k);
            }

            if(i==j){
                if(fabs(mu-1.0)>0.001){
                    printf("WARNING basis normalization error %e\n",mu);
                    exit(1);
                }
            }
            else{
                if(fabs(mu)>0.001){
                    printf("WARNING basis orthogonalization error %e\n",mu);
                    exit(1);
                }
            }
        }
    }

    while(_centers.get_rows()<_ncenters){
        for(i=0;i<_dim;i++){
            trial.set(i,_mins.get_data(i)+
                        0.5*_dice->doub()*(_maxs.get_data(i)-_mins.get_data(i))+
                        0.25*(_maxs.get_data(i)-_mins.get_data(i)));
        }
        _centers.add_row(trial);
    }

    double mean,stdev;

    while(_widths.get_rows()<_ncenters){
        for(i=0;i<_dim;i++){
            mu=0.0;
            while(mu<1.0e-20){
                if(i%2==0){
                    mean=_characteristic_width;
                    stdev=0.25*_characteristic_width;
                }
                else{
                    mean=10.0*_characteristic_width;
                    stdev=2.5*_characteristic_width;
                }
                mu=normal_deviate(_dice, mean, stdev);
            }
            trial.set(i,mu);
        }
        _widths.add_row(trial);
    }

    _chisq_initialized=1;

}

void chisquared::set_max(int dex, double nn){
    _maxs.set(dex,nn);
}

void chisquared::set_min(int dex, double nn){
    _mins.set(dex,nn);
}

double chisquared::get_min(int dex){
    return _mins.get_data(dex);
}

double chisquared::get_max(int dex){
    return _maxs.get_data(dex);
}

double chisquared::get_time_spent(){
    return _time_spent;
}


void chisquared::print_mins_maxs(){
    int i;
    double nn;
    nn=0.0;
    printf("mins and maxs\n");
    for(i=0;i<_dim;i++){
        printf("p%d -- %e %e -- %e\n",i,_mins.get_data(i),_maxs.get_data(i),_maxs.get_data(i)-_mins.get_data(i));
	nn+=power(_maxs.get_data(i)-_mins.get_data(i),2);
    }
    printf("\nfiducial distance %e\n",sqrt(nn));
}


double chisquared::get_width(int ic, int ix){

    return _widths.get_data(ic,ix);

}

double chisquared::get_center(int ic, int ix){
    return _centers.get_data(ic,ix);
}

double chisquared::get_real_center(int ic, int ix){
    if(ic>=_ncenters || ix>=_dim){
        return exception_value;
    }

    int i;
    double ans=0.0;
    for(i=0;i<_dim;i++){
        ans+=_centers.get_data(ic,i)*_bases.get_data(i,ix);
    }
    return ans;

}


void chisquared::death_knell(char *word)const{
    printf("%s\n",word);
    exit(1);
}

int chisquared::get_dim(){
    return _dim;
}

int chisquared::get_ncenters(){
    return _ncenters;
}

int chisquared::get_called(){
    return _called;
}

void chisquared::reset_timer(){
    _called=0;
    _time_spent=0.0;
}

double chisquared::operator()(const array_1d<double> &v){
    death_knell("meaningless chisq operator");
    return -1.0;
}

double chisquared::get_basis(int ix, int id){
    return _bases.get_data(ix,id);
}

void chisquared::get_basis(int ix, array_1d<double> &v){
    int i;
    if(ix<_dim){
        for(i=0;i<_dim;i++)v.set(i,_bases.get_data(ix,i));
    }
    else{
        printf("WARNING called get_basis with %d %d\n",ix,_dim);
	exit(1);
    }
}

void chisquared::project_to_basis(const array_1d<double> &in, array_1d<double> &out) const{
    int ix,iy;
    for(ix=0;ix<_dim;ix++){
        out.set(ix,0.0);
        for(iy=0;iy<_dim;iy++){
            out.add_val(ix,in.get_data(iy)*_bases.get_data(ix,iy));
        }
    }

}

double chisquared::project_to_basis(int ix, const array_1d<double> &vv) const{
    int i;
    double nn=1.0e30;
    if(ix<_dim){
        nn=0.0;
	for(i=0;i<_dim;i++)nn+=vv.get_data(i)*_bases.get_data(ix,i);
    }
    else{
        printf("WARNING called project_to_basis with %d %d\n",ix,_dim);
    }

    return nn;
}
