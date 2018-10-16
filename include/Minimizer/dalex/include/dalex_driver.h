#ifndef DALEX_DRIVER_H
#define DALEX_DRIVER_H

#include <stdio.h>
#include <time.h>
#include "eigen_wrapper.h"
#include "chisq_wrapper.h"
#include "simplex.h"
#include "cost_fn.h"
#include "dalex.h"
#include "dalex_initializer.h"

class dalex_driver{

public:

    dalex_driver();
    ~dalex_driver();

    void initialize(index_t);

    void set_seed(index_t);
    void set_min(array_1d<value_t>&);
    void set_max(array_1d<value_t>&);
    void set_characteristic_length(index_t,value_t);
    void set_deltachi(value_t);
    void set_target(value_t);
    void set_write_every(index_t);
    void set_outname(char*);
    void set_timingname(char*);

    void set_chisquared(chisquared*);

    void search(index_t);

    void mcmc_init();

    index_t get_dim();
    index_t get_called();
    value_t get_chimin();

    void set_confidence_limit(value_t);
    void set_dof(index_t);

    value_t evaluate(array_1d<value_t>&, index_t*);
    void assess_good_points();
    void assess_good_points(index_t);
    void assess_good_points(index_t,index_t);

private:

    chisq_wrapper _chifn;
    index_t _ct_dalex;

    dalex _cloud;

    char _outname[letters],_timingname[letters];

    array_1d<index_t> _good_points;

    index_t _last_did_min;

};

#endif
