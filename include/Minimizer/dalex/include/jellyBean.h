#ifndef JELLY_BEAN_H
#define JELLY_BEAN_H

#include "chisq.h"

class chiSquaredData : public chisquared{

    public:
        chiSquaredData(int, int, double, int, double);
        ~chiSquaredData(){};

        void set_width(int, int, double);

        virtual double operator()(const array_1d<double>&);
        void write_data();
        void print_mins();
        
        virtual void convert_params(const array_1d<double>&, array_1d<double>&, int);

    protected:
        int _ndata;
        double _sig;
        double _xmax,_dx;
        array_1d<double> _x_values,_y_values,_sigma;
        array_1d<double> _wave_phase,_wave_lambda,_env_d_center,_env_width;
        array_1d<double> _wave_amp;
        array_1d<double> _param_buffer;

        double data_function(array_1d<double>&, double);
        
        void initialize_data();

};

class jellyBeanData : public chiSquaredData{

    public:
        ~jellyBeanData(){}
        jellyBeanData(int,int,double,int,double);
        
        virtual void convert_params(const array_1d<double>&, array_1d<double>&, int);


    protected:
        double _parabola_curvature;
        array_2d<double> _parabola_centers,_parabola_x,_parabola_y;

};

class ellipseData : public chiSquaredData{

    public:
        ~ellipseData(){}
        ellipseData(int, int, int, double);

        virtual void convert_params(const array_1d<double>&, array_1d<double>&, int);
    
    protected:
        array_1d<double> _dir,_projected;
        

};

class nonGaussianEllipseData : public ellipseData{
    public:
        ~nonGaussianEllipseData(){}
        nonGaussianEllipseData(int, int, int, double);
        virtual void convert_params(const array_1d<double>&, array_1d<double>&, int);
};

#endif
