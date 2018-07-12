#ifndef DALEX_INIT_H
#define DALEX_INIT_H

#include "chisq_wrapper.h"
#include "simplex.h"

class dalex_initializer{

    public:
    
        ~dalex_initializer(){}
    
        dalex_initializer(){
            _chifn=NULL;
            _particles.set_name("dalex_init_particles");
            _local_min.set_name("dalex_init_local_min");
            _abs_min.set_name("dalex_init_abs_min");
            _since_min.set_name("dalex_init_since_min");
        }

        void set_chifn(chisq_wrapper *cc){
            _chifn=cc;
        }

        void safety_check(){
            if(_chifn==NULL){
                printf("WARNING _chifn is null\n");
                exit(1);
            }
        }

        double evaluate(array_1d<double> &pt, int *i_found, int i_point, double *mu_true_out){
            double mu;
            _chifn->evaluate(pt, &mu, i_found);
            mu_true_out[0]=mu;
            if(i_point>=0 && i_found[0]>=0){
                if(i_point>=_abs_min.get_dim() || mu_true_out[0]<_chifn->get_fn(_abs_min.get_data(i_point))){
                    _abs_min.set(i_point,i_found[0]);
                }
                if(i_point>=_local_min.get_dim() || mu_true_out[0]<_chifn->get_fn(_local_min.get_data(i_point))){
                    _local_min.set(i_point,i_found[0]);
                    _since_min.set(i_point,0);
                }
                else{
                    _since_min.add_val(i_point,1);
                }
            }
            return mu;
        }


        void search();

    private:
        chisq_wrapper *_chifn;
        array_1d<int> _particles,_local_min,_abs_min,_since_min;
};


#endif
