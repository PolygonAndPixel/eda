#ifndef EXPLORERS_H
#define EXPLORERS_H

#include "goto_tools.h"
#include "chisq_wrapper.h"
#include "cost_fn.h"
#include "ellipse.h"


class explorers{

    public:
        ~explorers(){}

        explorers(){
            _chifn=NULL;
            _mindex=-1;
            _temp=1.0;
            _scalar_acceptance=0;
            _attempted=0;
            _scalar_steps=0;
            _associates.set_name("explorers_associates");
            _median_associate.set_name("explorers_mean_associate");
            _particles.set_name("explorers_particles");
            _bases.set_name("explorers_bases");
            _norm.set_name("explorers_norm");
            _min.set_name("explorers_min");
            _max.set_name("explorers_max");
            _req_temp.set_name("explorers_req_temp");
            _mu_arr.set_name("explorers_mu_arr");
            _envelope=1.0;
            _target_rate=0.5;
        }

        void set_target_rate(value_t dd){
            _target_rate=dd;
        }

        void set_envelope(value_t dd){
            _envelope=dd;
        }

        void set_particle(index_t dex, array_1d<value_t> &pp){
            if(dex>=_particles.get_rows()){
                printf("WARNING no %d particle\n",dex);
                exit(1);
            }
            index_t i;
            for(i=0;i<pp.get_dim();i++){
                _particles.set(dex,i,pp.get_data(i));
            }
        }

        void set_chifn(chisq_wrapper *cc){
            _chifn=cc;
        }

        void set_associates(array_1d<index_t> &aa){
            index_t i;
            _associates.reset_preserving_room();
            for(i=0;i<aa.get_dim();i++){
                _associates.set(i,aa.get_data(i));
            }
        }

        void get_min_pt(array_1d<value_t> &pp){
            index_t i;
            for(i=0;i<_chifn->get_dim();i++){
                pp.set(i,_particles.get_data(_mindex,i));
            }
        }

        index_t get_n_particles(){
            return _particles.get_rows();
        }

        value_t get_mu(index_t dex){
            return _mu_arr.get_data(dex);
        }

        const array_1d<value_t> get_pt(index_t dex){
            return _particles(dex);
        }

        value_t get_pt(index_t dex, index_t i){
            return _particles.get_data(dex,i);
        }

        void get_pt(index_t dex, array_1d<value_t> &pp){
            index_t i;
            for(i=0;i<_chifn->get_dim();i++){
                pp.set(i,_particles.get_data(dex,i));
            }
        }

        void add_particle(const array_1d<value_t> &pt){
            _particles.add_row(pt);
        }

        void get_seed(array_2d<value_t>&);

        void set_norm();
        void reset();
        void initialize();
        void bump_particles();
        void kick(index_t);
        void sample(index_t,index_t);
        void sample(index_t);

    private:
        chisq_wrapper *_chifn;
        array_1d<index_t> _associates;
        array_1d<value_t> _median_associate;
        array_2d<value_t> _bases;
        array_1d<value_t> _norm,_min,_max;
        array_1d<value_t> _mu_arr;
        array_1d<value_t> _req_temp;
        index_t _mindex;
        value_t _mu_min;
        value_t _temp;
        index_t _attempted;
        array_2d<value_t> _particles;
        value_t _scalar_acceptance;
        value_t _scalar_steps;
        value_t _envelope;
        value_t _target_rate;
};

#endif
