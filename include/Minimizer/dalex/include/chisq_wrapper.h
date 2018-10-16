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

    void initialize(index_t);

    void write_pts();
    index_t get_search_type_log(index_t ii){
        return _search_type_log.get_data(ii);
    }
    index_t get_search_type(){
        return _search_type;
    }
    void set_search_type(index_t ii){
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
    void set_write_every(index_t ii){
        _write_every=ii;
    }
    void set_outname(char *nn){
        index_t i;
        for(i=0;i<letters-1 && nn[i]!=0;i++){
            _outname[i]=nn[i];
        }
        _outname[i]=0;
    }

    void set_timingname(char *nn){
        index_t i;
        for(i=0;i<letters-1 && nn[i]!=0;i++){
            _timingname[i]=nn[i];
        }
        _timingname[i]=0;
    }

    void set_chisquared(chisquared*);

    void set_target(value_t);
    void set_seed(index_t);
    void set_deltachi(value_t);
    void set_characteristic_length(index_t, value_t);
    value_t get_characteristic_length(index_t);
    void set_min(array_1d<value_t>&);
    void set_max(array_1d<value_t>&);
    void set_ddmin(value_t);

    index_t could_it_go_lower(value_t);

    value_t target();
    value_t chimin();
    index_t mindex();
    value_t get_deltachi();
    index_t get_pts();
    virtual index_t get_dim();
    virtual index_t get_called();
    virtual value_t get_time_spent();

    value_t random_double();
    index_t random_int();

    value_t raw_evaluate(const array_1d<value_t>&);
    void evaluate(const array_1d<value_t>&, value_t*, index_t*);
    virtual value_t operator()(const array_1d<value_t>&);
    value_t get_fn(index_t);
    value_t get_pt(index_t,index_t);
    array_1d<value_t> get_pt(index_t);

    value_t distance(array_1d<value_t>&,index_t);
    value_t distance(array_1d<value_t>&,array_1d<value_t>&);
    value_t distance(index_t,index_t);

    void nn_srch(array_1d<value_t>&,index_t,array_1d<index_t>&,array_1d<value_t>&);

    Ran* get_dice();

    void get_min(array_1d<value_t>&);
    void get_max(array_1d<value_t>&);
    virtual value_t get_min(index_t);
    virtual value_t get_max(index_t);

    void find_gradient(array_1d<value_t>&,array_1d<value_t>&);

    index_t in_bounds(const array_1d<value_t>&);
    index_t in_bounds(index_t, value_t);

    void set_confidence_limit(value_t);
    void set_dof(index_t);
    index_t get_seed();

    kd_tree* get_tree();
    array_1d<value_t>* get_fn_arr();

private:
    value_t _chimin,_deltachi,_target,_ddmin;
    value_t _expected_min,_expected_delta,_confidence_limit;
    index_t _adaptive_target,_seed,_called,_mindex,_iWhere;
    index_t _dof;
    index_t _search_type;
    array_1d<index_t> _search_type_log;

    array_1d<value_t> _characteristic_length,_range_min,_range_max,_fn;

    chisquared *_chifn;
    kd_tree *_kptr;
    Ran *_dice;

    index_t is_valid(const array_1d<value_t>&, index_t*);
    void is_it_safe(char*);

    array_1d<value_t> _valid_dd;
    array_1d<index_t> _valid_neigh;

    chisquared_distribution _distribution;

    char _outname[letters],_timingname[letters];
    index_t _last_written,_write_every;
    value_t _time_started,_last_time_spent,_time_batch;

};

#endif
