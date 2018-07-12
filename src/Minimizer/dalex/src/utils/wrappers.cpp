#include "wrappers.h"

function_wrapper::function_wrapper(){}

function_wrapper::~function_wrapper(){}

double function_wrapper::operator()(const array_1d<double> &vv){
    printf("WARNING calling un-implemented function_wrapper operator\n");
    exit(1);
}

int function_wrapper::get_called(){
    printf("WARNING calling un-implemented function_wrapper get_called\n");
    exit(1);
}
