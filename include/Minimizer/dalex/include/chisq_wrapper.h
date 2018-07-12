#ifndef CHISQ_WRAPPER_H
#define CHISQ_WRAPPER_H

#include <math.h>
#include <time.h>
#include <stdio.h>

#include "goto_tools.h"
#include "containers.h"
#include "kd.h"
#include "wrappers.h"
#include "chisq.h"

#define _type_init 0
#define _type_refine 1
#define _type_init_tendril 2
#define _type_tendril 3
#define _type_find_bases 4
#define _type_tendril_seed 5
#define _type_tendril_fill 6

class chisq_wrapper : public function_wrapper{

public:
    chisq_wrapper();
    ~chisq_wrapper();
    void copy(chisq_wrapper&);

    void initialize(int);

    void write_pts();
    int get_search_type_log(int ii){
        return _search_type_log.get_data(ii);
    }
    int get_search_type(){
        return _search_type;
    }
    void set_search_type(int ii){
        if(ii!=_type_init && ii!=_type_refine &&
           ii!=_type_init_tendril && ii!=_type_tendril &&
           ii!=_type_find_bases && ii!=_type_tendril_seed &&
           ii!=_type_tendril_fill && ii!=-1){

            if(ii>=0){
                printf("WARNING search type %d is not allowed\n", ii);
                exit(1);
            }
        }
        _search_type=ii;
    }
    void set_write_every(int ii){
        _write_every=ii;
    }
    void set_outname(char *nn){
        int i;
        for(i=0;i<letters-1 && nn[i]!=0;i++){
            _outname[i]=nn[i];
        }
        _outname[i]=0;
    }

    void set_timingname(char *nn){
        int i;
        for(i=0;i<letters-1 && nn[i]!=0;i++){
            _timingname[i]=nn[i];
        }
        _timingname[i]=0;
    }

    void set_chisquared(chisquared*);

    void set_target(double);
    void set_seed(int);
    void set_deltachi(double);
    void set_characteristic_length(int, double);
    double get_characteristic_length(int);
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_ddmin(double);

    int could_it_go_lower(double);

    double target();
    double chimin();
    int mindex();
    double get_deltachi();
    int get_pts();
    virtual int get_dim();
    virtual int get_called();
    virtual double get_time_spent();

    double random_double();
    int random_int();

    double raw_evaluate(const array_1d<double>&);
    void evaluate(const array_1d<double>&, double*, int*);
    virtual double operator()(const array_1d<double>&);
    double get_fn(int);
    double get_pt(int,int);
    array_1d<double> get_pt(int);

    double distance(array_1d<double>&,int);
    double distance(array_1d<double>&,array_1d<double>&);
    double distance(int,int);

    void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);

    Ran* get_dice();

    void get_min(array_1d<double>&);
    void get_max(array_1d<double>&);
    virtual double get_min(int);
    virtual double get_max(int);

    void find_gradient(array_1d<double>&,array_1d<double>&);

    int in_bounds(const array_1d<double>&);
    int in_bounds(int, double);

    void set_confidence_limit(double);
    void set_dof(int);
    int get_seed();

    kd_tree* get_tree();
    array_1d<double>* get_fn_arr();

private:
    double _chimin,_deltachi,_target,_ddmin;
    double _expected_min,_expected_delta,_confidence_limit;
    int _adaptive_target,_seed,_called,_mindex,_iWhere;
    int _dof;
    int _search_type;
    array_1d<int> _search_type_log;

    array_1d<double> _characteristic_length,_range_min,_range_max,_fn;

    chisquared *_chifn;
    kd_tree *_kptr;
    Ran *_dice;

    int is_valid(const array_1d<double>&, int*);
    void is_it_safe(char*);

    array_1d<double> _valid_dd;
    array_1d<int> _valid_neigh;

    chisquared_distribution _distribution;

    char _outname[letters],_timingname[letters];
    int _last_written,_write_every;
    double _time_started,_last_time_spent,_time_batch;

};

#endif
