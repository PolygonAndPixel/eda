#include "explorers.h"

void explorers::get_seed(array_2d<double> &seed){

    array_1d<int> mu_dex;
    mu_dex.set_name("exp_seed_mu_dex");
    array_1d<double> mu_sorted;
    mu_sorted.set_name("exp_mu_sorted");

    int i;
    for(i=0;i<_mu_arr.get_dim();i++){
        mu_dex.set(i,i);
    }
    sort(_mu_arr,mu_sorted,mu_dex);

    int j,k;
    double rr;

    for(i=0;i<_chifn->get_dim()+1;i++){
        seed.add_row(_particles(mu_dex.get_data(i)));
        for(j=0;j<_chifn->get_dim();j++){
            _particles.set(mu_dex.get_data(i),j,_median_associate.get_data(j));
        }
        for(j=0;j<_chifn->get_dim();j++){
            rr=normal_deviate(_chifn->get_dice(),0.0,2.0);
            for(k=0;k<_chifn->get_dim();k++){
                _particles.add_val(mu_dex.get_data(i),k,rr*_norm.get_data(j)*_bases.get_data(j,k));
            }
        }
    }
}

void explorers::set_norm(){

    int ip,ia;
    double component;
    int i,j;

    _min.reset_preserving_room();
    _max.reset_preserving_room();
    for(ip=0;ip<_associates.get_dim();ip++){
        ia=_associates.get_data(ip);
        for(i=0;i<_chifn->get_dim();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=_chifn->get_pt(ia,j)*_bases.get_data(i,j);
            }
            if(i>=_min.get_dim() || component<_min.get_data(i)){
                _min.set(i,component);
            }
            if(i>=_max.get_dim() || component>_max.get_data(i)){
                _max.set(i,component);
            }
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
       _median_associate.set(i,0.0);
    }

    for(i=0;i<_chifn->get_dim();i++){
        _norm.set(i,0.2*(_max.get_data(i)-_min.get_data(i)));
        for(j=0;j<_chifn->get_dim();j++){
             _median_associate.add_val(j,_bases.get_data(i,j)*0.5*(_max.get_data(i)+_min.get_data(i)));
        }
    }

}

void explorers::reset(){
    _particles.reset_preserving_room();
    _attempted=0;
    _scalar_acceptance=0.0;
    _scalar_steps=0.0;
    _req_temp.reset_preserving_room();
    _temp=1.0;
}

void explorers::initialize(){
    _particles.reset_preserving_room();
    _attempted=0;
    _req_temp.reset_preserving_room();
}


void explorers::bump_particles(){
    array_1d<double> dir, trial;
    dir.set_name("exp_bump_dir");
    trial.set_name("exp_bump_trial");
    _particles.reset_preserving_room();

    int i,j;
    while(_particles.get_rows()<get_n_particles()){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(_chifn->mindex(),i));
        }
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.add_val(j,1.0*(_max.get_data(i)-_min.get_data(i))
                               *dir.get_data(i)*_bases.get_data(i,j));
            }
        }
        _particles.add_row(trial);
    }
}


void explorers::kick(int dex){
     array_1d<double> dir;
     dir.set_name("explorers_kick_dir");
     int i;
     for(i=0;i<_chifn->get_dim();i++){
         dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
     }
     dir.normalize();
     for(i=0;i<_chifn->get_dim();i++){
         _particles.set(dex,i,_chifn->get_pt(_chifn->mindex(),i));
     }
     int j;
     for(i=0;i<_chifn->get_dim();i++){
         for(j=0;j<_chifn->get_dim();j++){
             _particles.add_val(dex,j,0.5*dir.get_data(i)*_bases.get_data(i,j)*(_max.get_data(i)-_min.get_data(i)));
         }
     }
}

void explorers::sample(int n_steps){
    sample(n_steps, 0);
}

void explorers::sample(int n_steps, int with_kick){

    cost_fn dchifn(_chifn, _associates);
    dchifn.set_envelope(_envelope);

    ellipse local_ellipse;
    array_2d<double> ellipse_pts;
    ellipse_pts.set_name("explorers_sample_ellipse_pts");
    int i,j;
    for(i=0;i<_associates.get_dim();i++){
        ellipse_pts.add_row(_chifn->get_pt(_associates.get_data(i)));
    }
    local_ellipse.build(ellipse_pts);
    _bases.reset_preserving_room();
    _bases.set_dim(_chifn->get_dim(),_chifn->get_dim());
    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            _bases.set(i,j,local_ellipse.bases(i,j));
        }
    }

    set_norm();


    if(with_kick==1){
        printf("\nkicking explorers\n");
        for(i=0;i<get_n_particles();i++){
            kick(i);
        }
    }

    _mindex=-1;

    double mu;
    int i_step;
    int i_dim,ip;
    int accept_it;
    double rr,roll,ratio;
    array_1d<double> trial;
    trial.set_name("exp_sample_trial");

    for(ip=0;ip<get_n_particles();ip++){
        mu=dchifn(_particles(ip));
        if(ip==0 || mu<_mu_min){
            _mindex=ip;
            _mu_min=mu;
        }
        _mu_arr.set(ip,mu);
    }


    double needed_temp;
    double old_temp;
    array_1d<int> req_temp_dex;
    array_1d<double> req_temp_sorted;
    req_temp_dex.set_name("exp_sample_req_temp_dex");
    req_temp_sorted.set_name("exp_sample_req_temp_sorted");

    int has_been_adjusted;

    printf("    starting sampling with %e -- %d\n",_mu_min,get_n_particles());

    for(i_step=0;i_step<n_steps;i_step++){

        for(ip=0;ip<get_n_particles();ip++){
            i_dim=_chifn->random_int()%_chifn->get_dim();
            rr=normal_deviate(_chifn->get_dice(),0.0,1.0);
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,_particles.get_data(ip,i)+
                          rr*_norm.get_data(i_dim)*_bases.get_data(i_dim,i));
            }
            accept_it=0;
            mu=dchifn(trial);
            if(mu<_mu_arr.get_data(ip)){
                accept_it=1;
                needed_temp=1.0;
            }
            else{
                ratio=exp(-0.5*(mu-_mu_arr.get_data(ip))/_temp);
                roll=_chifn->random_double();
                needed_temp=-0.5*(mu-_mu_arr.get_data(ip))/log(roll);
                if(needed_temp<1.0){
                    needed_temp=1.0;
                }
                if(roll<ratio){
                    accept_it=1;
                }
            }

            _req_temp.add(needed_temp);
            _scalar_steps+=1.0;
            if(accept_it==1){
                _scalar_acceptance+=1.0;
                _mu_arr.set(ip,mu);
                for(i=0;i<_chifn->get_dim();i++){
                    _particles.set(ip,i,trial.get_data(i));
                }
                if(mu<_mu_min){
                     _mindex=ip;
                     _mu_min=mu;
                }
            }

        }
        _attempted++;

        if(_attempted>0 && _attempted%(2*_chifn->get_dim())==0){
            printf("assessing temp\n");
            if(_scalar_acceptance>(_target_rate+0.1)*_scalar_steps ||
               _scalar_acceptance<(_target_rate-0.1)*_scalar_steps){
                old_temp=_temp;
                req_temp_sorted.reset_preserving_room();
                req_temp_dex.reset_preserving_room();
                for(i=0;i<_req_temp.get_dim();i++){
                    req_temp_dex.set(i,i);
                }
                sort(_req_temp, req_temp_sorted, req_temp_dex);
                _temp=req_temp_sorted.get_data(req_temp_dex.get_dim()/2);
                if(fabs(1.0-old_temp/_temp)<0.01){
                    if(_scalar_acceptance>(_target_rate+0.1)*_scalar_steps){
                        _temp*=0.8;
                    }
                    else{
                        _temp*=1.25;
                    }
                    if(_temp<1.0){
                         _temp=1.0;
                    }
                }
                if(fabs(1.0-old_temp/_temp)>0.01){
                    _req_temp.reset();
                    _attempted=0;
                    _scalar_acceptance=0.0;
                    _scalar_steps=0.0;
                }
            }
        }

    }

    printf("    sampling min %e\n",_mu_min);
    printf("    temp %e\n",_temp);
    printf("    accepted %.1f steps %.1f\n",_scalar_acceptance,
    _scalar_steps);
    printf("    ct %d\n",_chifn->get_pts());
}
