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
    virtual value_t operator()(const array_1d<value_t>&);
    virtual index_t get_called();
    virtual value_t get_min(index_t i){
        printf("unimplemented get_min in function_wrapper");
        exit(1);
    }

    virtual value_t get_max(index_t i){
        printf("unimplemented get_max in function_wrapper");
        exit(1);
    }

    virtual index_t get_dim(){
        printf("unimplemented get_dim in function_wrapper");
        exit(1);
    }

    virtual value_t get_time_spent(){
        printf("unimplemented get_time_spent in function_wrapper");
        exit(1);
    }
};

#endif
