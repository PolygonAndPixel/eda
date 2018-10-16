#ifndef COST_FN_H
#define COST_FN_H

#include "chisq_wrapper.h"

class cost_fn : public function_wrapper{
    public:
        cost_fn(chisq_wrapper*, array_1d<index_t>&);
        cost_fn(chisq_wrapper*, array_1d<index_t>&, index_t);
        cost_fn(){_chifn=NULL;};
        ~cost_fn(){};
        void build(chisq_wrapper*, array_1d<index_t>&, index_t);
        virtual value_t operator()(const array_1d<value_t>&);
        virtual index_t get_called();
        void multiply_norm(value_t dd){
            _scalar_norm*=dd;
        }
        value_t nn_distance(const array_1d<value_t>&);
        void set_envelope(value_t dd){
            _envelope=dd;
        }

        index_t get_cached_values(const index_t dex, value_t *fn){

            index_t i;
            for(i=0;i<_pt_cache.get_dim();i++){
                if(_pt_cache.get_data(i)==dex){
                    fn[0]=_fn_cache.get_data(i);
                    return 1;
                }
            }
            return 0;

        }

    private:
        array_1d<index_t> _associates;
        array_1d<value_t> _median_associate;
        value_t _scalar_norm;
        chisq_wrapper *_chifn;
        index_t _called;
        value_t _envelope;

        array_1d<index_t> _pt_cache;
        array_1d<value_t> _fn_cache;

};

#endif
