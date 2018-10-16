#ifndef CHISQ_H
#define CHISQ_H

#include "goto_tools.h"
#include "containers.h"
#include "kd.h"
#include "wrappers.h"
#include <math.h>
#include <stdio.h>

class chisquared : public function_wrapper{

public:
    chisquared();
    chisquared(index_t);
    chisquared(index_t,index_t);
    chisquared(index_t,index_t,value_t);
    ~chisquared(){
        if(_dice!=NULL){
            delete _dice;
        }
    }

    virtual value_t operator()(const array_1d<value_t>&);

    virtual index_t get_called();
    void reset_timer();

    virtual value_t get_time_spent();

    void set_max(index_t,value_t);
    void set_min(index_t,value_t);

    virtual value_t get_min(index_t);
    virtual value_t get_max(index_t);

    void get_basis(index_t,array_1d<value_t>&);
    value_t get_basis(index_t,index_t);

    value_t project_to_basis(index_t,const array_1d<value_t>&) const;
    void project_to_basis(const array_1d<value_t>&, array_1d<value_t>&) const;

    value_t get_width(index_t,index_t);

    value_t get_center(index_t,index_t);

    value_t get_real_center(index_t,index_t);

    void print_mins_maxs();

    index_t get_ncenters();

    virtual index_t get_dim();

    void enable_logging(){
        _with_logging=1;
    }

    array_2d<value_t> pt_log;
    array_1d<value_t> fn_log;

protected:

    index_t _dim,_ncenters,_seed;
    value_t _characteristic_width;
    index_t _chisq_initialized;


    mutable index_t _called;

    mutable value_t _time_spent;

    array_1d<value_t> _maxs,_mins;

    array_2d<value_t> _bases,_widths,_centers;

    Ran *_dice;

    void death_knell(char*) const;

    void _initialize();
    void initialize();

    index_t _with_logging;

    void _log_point(const array_1d<value_t> &pt, value_t mu){
        if(_with_logging==0){
             return;
        }

        fn_log.add(mu);
        pt_log.add_row(pt);
    }

};

#endif
