#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "containers.h"
#include "goto_tools.h"
#include "wrappers.h"

class simplex_minimizer{

public:

    simplex_minimizer();
    simplex_minimizer(value_t,value_t,value_t);
    ~simplex_minimizer();

    void set_chisquared(function_wrapper*);
    void set_cost(function_wrapper*);
    void set_dice(Ran*);
    void set_minmax(array_1d<value_t>&, array_1d<value_t>&);
    void set_abort_max_factor(index_t);
    void use_gradient();
    void do_not_use_gradient();
    void freeze_temp();
    void unfreeze_temp();
    void get_minpt(array_1d<value_t>&);
    void get_pt(index_t,array_1d<value_t>&);
    void is_a_model();

    /*
    the array_2d will be the input array of points;
    the array_1d will be the output minimum point
    */
    void find_minimum(array_2d<value_t>&, array_1d<value_t>&);
    value_t get_minimum();

    void set_limit(index_t ii){
        _limit=ii;
    }

private:

    value_t evaluate(const array_1d<value_t>&);
    value_t evaluate_cost(const array_1d<value_t>&);
    void cool_off();

    void gradient_minimizer();
    void gradient_cloud();
    void calculate_gradient(const array_1d<value_t>&, array_1d<value_t>&);
    value_t get_dx(index_t);
    void expand();

    void find_il();
    void paranoia();

    void initialize();

    void is_it_safe(char*);

    value_t _temp,_min_ff,_true_min_ff,_fstar,_fstarstar,_min_temp;
    value_t _alpha,_beta,_gamma;
    index_t _il,_ih,_called_cost,_freeze_temp,_use_gradient;
    index_t _freeze_called,_last_called_gradient;
    index_t _last_found,_called_evaluate,_abort_max_factor;
    index_t _last_cooled_off;
    index_t _is_a_model;
    index_t _limit;
    index_t _n_gradients;
    array_1d<value_t> _transform, _origin,_ff,_pstar,_pstarstar,_min_pt;
    array_1d<value_t> _last_improved_ff;
    array_1d<value_t> _min,_max;
    array_2d<value_t> _pts,_last_improved_pts;

    Ran *_dice;
    function_wrapper *_chisquared, *_cost;

    /*
    cost will need to be a function_wrapper sub-class
    that has pointers to all of the aps nodes.
    */

};

#endif
