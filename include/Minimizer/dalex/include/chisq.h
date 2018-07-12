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
    chisquared(int);
    chisquared(int,int);
    chisquared(int,int,double);
    ~chisquared(){
        if(_dice!=NULL){
            delete _dice;
        }
    }

    virtual double operator()(const array_1d<double>&);

    virtual int get_called();
    void reset_timer();

    virtual double get_time_spent();

    void set_max(int,double);
    void set_min(int,double);

    virtual double get_min(int);
    virtual double get_max(int);

    void get_basis(int,array_1d<double>&);
    double get_basis(int,int);

    double project_to_basis(int,const array_1d<double>&) const;
    void project_to_basis(const array_1d<double>&, array_1d<double>&) const;

    double get_width(int,int);

    double get_center(int,int);

    double get_real_center(int,int);

    void print_mins_maxs();

    int get_ncenters();

    virtual int get_dim();

    void enable_logging(){
        _with_logging=1;
    }

    array_2d<double> pt_log;
    array_1d<double> fn_log;

protected:

    int _dim,_ncenters,_seed;
    double _characteristic_width;
    int _chisq_initialized;


    mutable int _called;

    mutable double _time_spent;

    array_1d<double> _maxs,_mins;

    array_2d<double> _bases,_widths,_centers;

    Ran *_dice;

    void death_knell(char*) const;

    void _initialize();
    void initialize();

    int _with_logging;

    void _log_point(const array_1d<double> &pt, double mu){
        if(_with_logging==0){
             return;
        }

        fn_log.add(mu);
        pt_log.add_row(pt);
    }

};

#endif
