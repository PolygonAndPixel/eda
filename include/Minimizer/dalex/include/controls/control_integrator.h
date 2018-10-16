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

        virtual void initialize(array_1d<value_t> &min,
                                array_1d<value_t> &max,
                                array_1d<value_t> &dx){

            printf("calling initialize on iteration_parameters base class\n");
            exit(1);

        }

        virtual index_t get_pt(array_1d<value_t> &pt, array_1d<index_t> &idx){
            printf("calling get_pt on iteration_parameters base class\n");
            exit(1);
            return 0;
        }

        virtual long index_t get_current_ct(){
            return _current_ct;
        }

        virtual long index_t get_total_ct(){
            printf("cannot call default get_total_ct\n");
            exit(1);
        }


        protected:
            long index_t _current_ct;

};


class default_iteration_parameters : public iteration_parameters{

    public:

        default_iteration_parameters(){
            _grid_ct.set_name("default_grid_ct");
        }
        ~default_iteration_parameters(){}

        virtual void initialize(array_1d<value_t> &min,
                               array_1d<value_t> &max,
                               array_1d<value_t> &dx){

            _grid_ct.reset();

            index_t ix, iy;
            _total_ct=1;
            _current_ct=0;
            for(ix=0;ix<min.get_dim();ix++){
                _min.set(ix,min.get_data(ix));
                _max.set(ix,max.get_data(ix));
                _dx.set(ix,dx.get_data(ix));
                iy=index_t((max.get_data(ix)-min.get_data(ix))/(dx.get_data(ix)));
                _grid_ct.set(ix, iy);
                _total_ct*=iy;
                printf("factor %d %e %e %e\n",
                       iy,min.get_data(ix),max.get_data(ix),dx.get_data(ix));
            }
            printf("total_ct %ld\n",_total_ct);

        }

        virtual index_t get_pt(array_1d<value_t> &pt, array_1d<index_t> &idx){
            expand_grid(_current_ct, _grid_ct, idx);
            index_t ix;
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

        virtual long index_t get_total_ct(){
            return _total_ct;
        }

    private:
        array_1d<index_t> _grid_ct;
        array_1d<value_t> _min,_max,_dx;
        long index_t _total_ct;

};


class control_integrator{

    public:
        ~control_integrator(){
            if(_iter_p!=NULL){
                delete _iter_p;
            }
        }

        control_integrator(function_wrapper&,array_1d<value_t>&,array_1d<value_t>&,array_1d<value_t>&,char*);

        void set_d_threshold(value_t dd){
            _d_threshold=dd;
        }

        void set_max_chi_lim_freq(value_t mm){
            _max_chi_lim_freq=mm;
        }

        void set_chi_lim_freq(value_t);

        void run_analysis(array_1d<value_t>&);
        void run_analysis(value_t);
        void run_analysis();

        index_t get_dex(value_t, value_t, value_t, index_t);
        void write_output(index_t, index_t,
                        array_1d<value_t>&,
                        array_1d<value_t>&,
                        array_1d<index_t>&, array_2d<index_t>&,
                        value_t,
                        value_t,
                        char*);

        void find_chi_min(array_1d<value_t>&, array_1d<value_t>&);

        value_t get_min(index_t);
        value_t get_max(index_t);

    protected:
        array_1d<value_t> _min,_max,_dx;
        value_t _chi_min;
        value_t _d_threshold;
        value_t _max_chi_lim_freq;
        function_wrapper *_chisq;
        char _name_root[letters];

        virtual void _initialize_iterate();

        iteration_parameters *_iter_p;
};

#endif
