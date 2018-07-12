#ifndef COST_FN_H
#define COST_FN_H

#include "chisq_wrapper.h"

class cost_fn : public function_wrapper{
    public:
        cost_fn(chisq_wrapper*, array_1d<int>&);
        cost_fn(chisq_wrapper*, array_1d<int>&, int);
        cost_fn(){_chifn=NULL;};
        ~cost_fn(){};
        void build(chisq_wrapper*, array_1d<int>&, int);
        virtual double operator()(const array_1d<double>&);
        virtual int get_called();
        void multiply_norm(double dd){
            _scalar_norm*=dd;
        }
        double nn_distance(const array_1d<double>&);
        void set_envelope(double dd){
            _envelope=dd;
        }

        int get_cached_values(const int dex, double *fn){

            int i;
            for(i=0;i<_pt_cache.get_dim();i++){
                if(_pt_cache.get_data(i)==dex){
                    fn[0]=_fn_cache.get_data(i);
                    return 1;
                }
            }
            return 0;

        }

    private:
        array_1d<int> _associates;
        array_1d<double> _median_associate;
        double _scalar_norm;
        chisq_wrapper *_chifn;
        int _called;
        double _envelope;

        array_1d<int> _pt_cache;
        array_1d<double> _fn_cache;

};

#endif
