#ifndef CONTROL_INTEGRATOR_H
#define CONTROL_INTEGRATOR_H

#include "goto_tools.h"
#include "containers.h"
#include "wrappers.h"
#include "kd.h"
#include "simplex.h"

class iteration_parameters{

    public:
        iteration_parameters(){}
        ~iteration_parameters(){}

        virtual void initialize(array_1d<double> &min,
                                array_1d<double> &max,
                                array_1d<double> &dx){

            printf("calling initialize on iteration_parameters base class\n");
            exit(1);

        }

        virtual int get_pt(array_1d<double> &pt, array_1d<int> &idx){
            printf("calling get_pt on iteration_parameters base class\n");
            exit(1);
            return 0;
        }

        virtual long int get_current_ct(){
            return _current_ct;
        }

        virtual long int get_total_ct(){
            printf("cannot call default get_total_ct\n");
            exit(1);
        }


        protected:
            long int _current_ct;

};


class default_iteration_parameters : public iteration_parameters{

    public:

        default_iteration_parameters(){
            _grid_ct.set_name("default_grid_ct");
        }
        ~default_iteration_parameters(){}

        virtual void initialize(array_1d<double> &min,
                               array_1d<double> &max,
                               array_1d<double> &dx){

            _grid_ct.reset();

            int ix, iy;
            _total_ct=1;
            _current_ct=0;
            for(ix=0;ix<min.get_dim();ix++){
                _min.set(ix,min.get_data(ix));
                _max.set(ix,max.get_data(ix));
                _dx.set(ix,dx.get_data(ix));
                iy=int((max.get_data(ix)-min.get_data(ix))/(dx.get_data(ix)));
                _grid_ct.set(ix, iy);
                _total_ct*=iy;
                printf("factor %d %e %e %e\n",
                       iy,min.get_data(ix),max.get_data(ix),dx.get_data(ix));
            }
            printf("total_ct %ld\n",_total_ct);

        }

        virtual int get_pt(array_1d<double> &pt, array_1d<int> &idx){
            expand_grid(_current_ct, _grid_ct, idx);
            int ix;
            for(ix=0;ix<_min.get_dim();ix++){
                pt.set(ix,_min.get_data(ix)+idx.get_data(ix)*_dx.get_data(ix));
            }

            _current_ct++;
            if(_current_ct<_total_ct){
                return 1;
            }
            else{
                return 0;
            }
        }

        virtual long int get_total_ct(){
            return _total_ct;
        }

    private:
        array_1d<int> _grid_ct;
        array_1d<double> _min,_max,_dx;
        long int _total_ct;

};


class control_integrator{

    public:
        ~control_integrator(){
            if(_iter_p!=NULL){
                delete _iter_p;
            }
        }

        control_integrator(function_wrapper&,array_1d<double>&,array_1d<double>&,array_1d<double>&,char*);

        void set_d_threshold(double dd){
            _d_threshold=dd;
        }

        void set_max_chi_lim_freq(double mm){
            _max_chi_lim_freq=mm;
        }

        void set_chi_lim_freq(double);

        void run_analysis(array_1d<double>&);
        void run_analysis(double);
        void run_analysis();

        int get_dex(double, double, double, int);
        void write_output(int, int,
                        array_1d<double>&,
                        array_1d<double>&,
                        array_1d<int>&, array_2d<int>&,
                        double,
                        double,
                        char*);

        void find_chi_min(array_1d<double>&, array_1d<double>&);

        double get_min(int);
        double get_max(int);

    protected:
        array_1d<double> _min,_max,_dx;
        double _chi_min;
        double _d_threshold;
        double _max_chi_lim_freq;
        function_wrapper *_chisq;
        char _name_root[letters];

        virtual void _initialize_iterate();

        iteration_parameters *_iter_p;
};

#endif
