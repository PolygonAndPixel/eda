#include "wrappers.h"

function_wrapper::function_wrapper(){}

function_wrapper::~function_wrapper(){}

value_t function_wrapper::operator()(const array_1d<value_t> &vv){
    printf("WARNING calling un-implemented function_wrapper operator\n");
    exit(1);
}

index_t function_wrapper::get_called(){
    printf("WARNING calling un-implemented function_wrapper get_called\n");
    exit(1);
}
