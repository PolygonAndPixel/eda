/*
This file provides the chisquared function wmap_likelihood which allows APS to
interface with the WMAP 7 year likelihood function.

In order to compile with this option, the user will have to download CAMB and
the WMAP likelihood code and adjust the library and include paths in the
Makefile accordintly

*/

#ifndef WMAP_LIKELIHOOD_H
#define WMAP_LIKELIHOOD_H

#include "chisq.h"
#include <math.h>
#include <stdio.h>


/*
the two external function definitions below are only needed to
interface with CAMB and the WMAP 7 likelihood function
as provided in aps_cmb_module.cpp
*/

extern "C" void \
camb_wrap_(value_t*,value_t*,value_t*,value_t*,value_t*,value_t*);

extern "C" void wmaplikeness_(value_t*,value_t*,value_t*,value_t*,value_t*);

class wmap_likelihood : public chisquared{

public:
    wmap_likelihood();
    ~wmap_likelihood();
    virtual value_t operator()(array_1d<value_t>&);
};


class wmap_2d_likelihood : public chisquared{

public:
    wmap_2d_likelihood();
    ~wmap_2d_likelihood();
    virtual value_t operator()(array_1d<value_t>&);

private:
    wmap_likelihood _wmap;

};

#endif
