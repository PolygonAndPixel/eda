#ifndef JELLY_BEAN_H
#define JELLY_BEAN_H

#include "chisq.h"

class chiSquaredData : public chisquared{

    public:
        chiSquaredData(index_t, index_t, value_t, index_t, value_t);
        ~chiSquaredData(){};

        void set_width(index_t, index_t, value_t);

        virtual value_t operator()(const array_1d<value_t>&);
        void write_data();
        void print_mins();
        
        virtual void convert_params(const array_1d<value_t>&, array_1d<value_t>&, index_t);

    protected:
        index_t _ndata;
        value_t _sig;
        value_t _xmax,_dx;
        array_1d<value_t> _x_values,_y_values,_sigma;
        array_1d<value_t> _wave_phase,_wave_lambda,_env_d_center,_env_width;
        array_1d<value_t> _wave_amp;
        array_1d<value_t> _param_buffer;

        value_t data_function(array_1d<value_t>&, value_t);
        
        void initialize_data();

};

class jellyBeanData : public chiSquaredData{

    public:
        ~jellyBeanData(){}
        jellyBeanData(index_t,index_t,value_t,index_t,value_t);
        
        virtual void convert_params(const array_1d<value_t>&, array_1d<value_t>&, index_t);


    protected:
        value_t _parabola_curvature;
        array_2d<value_t> _parabola_centers,_parabola_x,_parabola_y;

};

class ellipseData : public chiSquaredData{

    public:
        ~ellipseData(){}
        ellipseData(index_t, index_t, index_t, value_t);

        virtual void convert_params(const array_1d<value_t>&, array_1d<value_t>&, index_t);
    
    protected:
        array_1d<value_t> _dir,_projected;
        

};

class nonGaussianEllipseData : public ellipseData{
    public:
        ~nonGaussianEllipseData(){}
        nonGaussianEllipseData(index_t, index_t, index_t, value_t);
        virtual void convert_params(const array_1d<value_t>&, array_1d<value_t>&, index_t);
};

#endif
