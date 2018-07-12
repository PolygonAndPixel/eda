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

        void set_target_rate(double dd){
            _target_rate=dd;
        }

        void set_envelope(double dd){
            _envelope=dd;
        }

        void set_particle(int dex, array_1d<double> &pp){
            if(dex>=_particles.get_rows()){
                printf("WARNING no %d particle\n",dex);
                exit(1);
            }
            int i;
            for(i=0;i<pp.get_dim();i++){
                _particles.set(dex,i,pp.get_data(i));
            }
        }

        void set_chifn(chisq_wrapper *cc){
            _chifn=cc;
        }

        void set_associates(array_1d<int> &aa){
            int i;
            _associates.reset_preserving_room();
            for(i=0;i<aa.get_dim();i++){
                _associates.set(i,aa.get_data(i));
            }
        }

        void get_min_pt(array_1d<double> &pp){
            int i;
            for(i=0;i<_chifn->get_dim();i++){
                pp.set(i,_particles.get_data(_mindex,i));
            }
        }

        int get_n_particles(){
            return _particles.get_rows();
        }

        double get_mu(int dex){
            return _mu_arr.get_data(dex);
        }

        const array_1d<double> get_pt(int dex){
            return _particles(dex);
        }

        double get_pt(int dex, int i){
            return _particles.get_data(dex,i);
        }

        void get_pt(int dex, array_1d<double> &pp){
            int i;
            for(i=0;i<_chifn->get_dim();i++){
                pp.set(i,_particles.get_data(dex,i));
            }
        }

        void add_particle(const array_1d<double> &pt){
            _particles.add_row(pt);
        }

        void get_seed(array_2d<double>&);

        void set_norm();
        void reset();
        void initialize();
        void bump_particles();
        void kick(int);
        void sample(int,int);
        void sample(int);

    private:
        chisq_wrapper *_chifn;
        array_1d<int> _associates;
        array_1d<double> _median_associate;
        array_2d<double> _bases;
        array_1d<double> _norm,_min,_max;
        array_1d<double> _mu_arr;
        array_1d<double> _req_temp;
        int _mindex;
        double _mu_min;
        double _temp;
        int _attempted;
        array_2d<double> _particles;
        double _scalar_acceptance;
        double _scalar_steps;
        double _envelope;
        double _target_rate;
};

#endif
