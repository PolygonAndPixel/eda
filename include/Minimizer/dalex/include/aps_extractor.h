/*
This file defines the class that converts the outputs of APS into Frequentist
and Bayesian limits on parameter space.

For an example of how to use this class and its methods, se

aps_extraction_runner.cpp

*/

#include "goto_tools.h"
#include "kd.h"

#ifndef APS_EXTRACTOR_H
#define APS_EXTRACTOR_H

class aps_extractor{
    /*
    This class will store the name of the file in which the APS output being
    processed is stored in the global variable filenam.
    
    The file will automatically learn the number of parameters in parameter
    space (provided that the user has not changed the columns output by
    aps.write_pts() in aps.cpp)
    */
    
    public:
        aps_extractor();
        ~aps_extractor();
        
        /*set the name of the file where the APS output is stored*/
        void set_filename(char*);
        
        /*if chisquared_lim = chisquared_min + delta chisquared, set 
        delta chisquared*/
        void set_delta_chi(value_t);
        
        /*write all points for which chisquared<chisquared_lim to the specified
        file*/
        void write_good_points(char*);
        
        /*same as above, but for a 2-dimensional subspace of parameter space*/
        void write_good_points(char*,index_t,index_t);
        
        /*write a file that charts how the minimum discovered value of 
        chisquared changes as a function of the number of samples*/
        void plot_chimin(char*);
        
        /*print the minimum value of chisquared and the corresponding
        point in parameter space*/
        void show_minpt();
        
        /*
        Assuming that the hyperbox method of drawing the Bayesian constraints
        described in the paper works, draw random samples from that posterior
        */
        void sample_posterior(char*,index_t);
        void sample_posterior(array_2d<value_t>&,index_t);
        void sample_posterior(char*,array_2d<value_t>&,index_t,index_t);
        
        /*draw bayesian constraints on a 2-dimensional subspace as described
        in the paper*/
        void draw_bayesian_bounds(char*,index_t,index_t,value_t);
        
        /*set a maximum number of APS samples to use when drawing constraints*/
        void set_cutoff(index_t);
        
        /*set the target value of chisquared by hand (i.e. do not use the
        chisquared_lim = chisquared_min + delta chisquared formalism)*/
        void set_target(value_t);

    private:

        char filename[letters];
        value_t chi_target;
        value_t chi_min,delta_chi,global_tol;
        index_t nparams,extra_words;
        index_t cutoff,asserted;
        
        array_1d<index_t> l_prob_dexes;
        array_1d<value_t> min_pt,l_probability,chisq;
        array_2d<value_t> box_max,box_min;
        
        /*find the minimum value of chisquared in the set of samples*/
        void learn_chimin();
        
        /*find the number of parameters in parameter space*/
        void learn_nparams();
        
        /*make sure that filename has been set before you try to do any
        real work*/
        void validate();
        
        /*thin out the points plotted by write_good_points and
        draw_bayesian_bounds so that plots of the data do not
        create obscenely large image files*/
        void plot_thinned_data(array_2d<value_t>&,char*);
        
        /*backend for draw_bayesian_bounds and sample_posterior*/
        void make_boxes();
};

#endif
