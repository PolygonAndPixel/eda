#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <math.h>
#include <time.h>
#include <stdio.h>

#include "containers.h"

class function_wrapper{
public:
    function_wrapper();
    ~function_wrapper();
    virtual double operator()(const array_1d<double>&);
    virtual int get_called();
    virtual double get_min(int i){
        printf("unimplemented get_min in function_wrapper");
        exit(1);
    }

    virtual double get_max(int i){
        printf("unimplemented get_max in function_wrapper");
        exit(1);
    }

    virtual int get_dim(){
        printf("unimplemented get_dim in function_wrapper");
        exit(1);
    }

    virtual double get_time_spent(){
        printf("unimplemented get_time_spent in function_wrapper");
        exit(1);
    }
};

#endif
