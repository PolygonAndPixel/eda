//#include <iostream.h>

#include "goto_tools.h"
#include "containers.h"

#ifndef EIGEN_H
#define EIGEN_H

extern "C" void dsytrd_(char*,index_t*,value_t*,index_t*,value_t*,value_t*,value_t*,value_t*,index_t*,index_t*);

// extern "C" void dstebz_(char*,char*,index_t*,value_t*,value_t*,index_t*,index_t*,value_t*,value_t*,value_t*,index_t*,index_t*,value_t*,index_t*,index_t*,value_t*,index_t*,index_t*);

// extern "C" void dstein_(index_t*,value_t*,value_t*,index_t*,value_t*,index_t*,index_t*,value_t*,index_t*,value_t*,index_t*,index_t*,index_t*);

// extern "C" void dormtr_(char*,char*,char*,index_t*,index_t*,value_t*,index_t*,value_t*,value_t*,index_t*,value_t*,index_t*,index_t*);

extern "C" void dgetrf_(index_t*,index_t*,value_t*,index_t*,index_t*,index_t*);

// extern "C" void dgetrs_(char*,index_t*,index_t*,value_t*,index_t*,index_t*,value_t*,index_t*,index_t*);

// extern "C" void dgetri_(index_t*,value_t*,index_t*,index_t*,value_t*,index_t*,index_t*);

extern "C" void dsaupd_(index_t*,char*,index_t*,char*,index_t*,value_t*,value_t*,index_t*,value_t*,index_t*,index_t*,index_t*,value_t*,value_t*,index_t*,index_t*);

extern "C" void dseupd_(index_t*,char*,index_t*,value_t*,value_t*,index_t*,value_t*,char*,index_t*,char*,index_t*,value_t*,value_t*,index_t*,value_t*,index_t*,index_t*,index_t*,value_t*,value_t*,index_t*,index_t*);

void matrix_multiply(value_t**, index_t, index_t, value_t**, index_t, index_t, value_t**);

value_t check_inversion(array_2d<value_t>&,array_2d<value_t>&);

value_t determinant(value_t*,index_t);

value_t trace(value_t*,index_t);

void invert_lapack(array_2d<value_t>&,array_2d<value_t>&,index_t);

void solve_lapack_nbyn(array_2d<value_t>&,array_1d<value_t>&,array_1d<value_t>&);

void eval_symm(array_2d<value_t>&, array_2d<value_t>&, array_1d<value_t>&);

void eval_symm(array_2d<value_t>&, array_2d<value_t>&, array_1d<value_t>&, value_t);

void eval_symm_guts(array_2d<value_t>&,array_2d<value_t>&,array_1d<value_t>&,index_t,index_t,index_t);

void eval_symm_guts(array_2d<value_t>&,array_2d<value_t>&,array_1d<value_t>&,index_t,index_t,index_t,value_t);

value_t eigen_check(array_2d<value_t>&,array_1d<value_t>&,value_t,index_t);

value_t eigen_check_open(value_t**,value_t*,index_t);

void eigen_solve(value_t**,index_t,index_t,value_t*, value_t**);

#endif
