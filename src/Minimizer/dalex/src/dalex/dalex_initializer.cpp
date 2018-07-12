#include "dalex_initializer.h"

void dalex_initializer::search(){
    safety_check();
    int total_per=100*_chifn->get_dim();
    int adjust_every=total_per/10;
    int n_groups=2;
    int n_particles=n_groups*(_chifn->get_dim()+1)+_chifn->get_dim()/2;
    double global_norm=0.05;

    array_1d<int> accepted,accepted_sorted,accepted_dex;
    array_1d<int> total_accepted;
    accepted.set_name("mcmc_init_accepted");
    accepted_sorted.set_name("mcmc_init_acc_sorted");
    accepted_dex.set_name("mcmc_init_acc_dex");
    total_accepted.set_name("total_accepted");

    double _temp=1.0;

    array_1d<double> trial,dir,norm;
    trial.set_name("mcmc_init_trial");
    dir.set_name("mcmc_init_dir");
    norm.set_name("mcmc_init_norm");

    double rr;

    int ip,i,j,k,i_step,i_found;

    for(i=0;i<_chifn->get_dim();i++){
        norm.set(i,_chifn->get_characteristic_length(i));
    }

    int adjusted=0;
    double mu,mu_true,roll,ratio;
    int accept_it;
    int min_acc,max_acc,med_acc;

    array_1d<double> min_vals,min_val_sorted;
    array_1d<int> min_dexes;
    min_vals.set_name("min_vals");
    min_val_sorted.set_name("min_val_sorted");
    min_dexes.set_name("min_dexes");

    array_1d<int> current_particles;
    current_particles.set_name("mcmc_init_current_particles");

    asymm_array_2d<int> trails;
    trails.set_name("mcmc_init_trails");

    for(i=0;i<n_particles;i++){
        accepted.set(i,0);
        total_accepted.set(i,0);
    }

    int i_best;
    int n_cand=100;

    array_1d<double> c_v,c_v_s;
    array_1d<int> c_v_d;

    for(ip=0;ip<n_particles;ip++){
        _chifn->set_search_type(-1*(ip+2));
        i_best=-1;

        while(i_best<0){
            for(i=0;i<_chifn->get_dim();i++){
                dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
            dir.normalize();
            for(i=0;i<_chifn->get_dim();i++){
                 mu=_chifn->get_max(i)-_chifn->get_min(i);
                 trial.set(i,0.5*(_chifn->get_max(i)+_chifn->get_min(i)));
                 trial.add_val(i,0.5*mu*dir.get_data(i));
            }
            mu=evaluate(trial,&i_best,ip,&mu_true);

        }
        _particles.set(ip,i_best);
        trails.set(ip,0,i_best);
    }


    double needed_temp;
    array_1d<double> needed_temp_arr,needed_temp_sorted;
    array_1d<int> needed_temp_dex;
    double old_temp;
    int has_been_adjusted;
    int step_ct=0;
    int needs_adjustment;

    array_1d<double> geo_center,local_min,local_max;
    int min_pt_connected;
    double dd_min,dd_best,dd;

    int n_jumps=0;
    int n_opt_out=0;

    array_2d<double> bases;
    array_1d<double> vv;

    int i_dim,i_half;
    double sgn;

    printf("starting steps with min %e\n",_chifn->chimin());
    for(i_step=0;i_step<total_per;i_step++){
        if(i_step%(4*_chifn->get_dim())==0){
            bases.reset_preserving_room();
            while(bases.get_rows()!=_chifn->get_dim()){
                for(i=0;i<_chifn->get_dim();i++){
                    vv.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
                }
                for(i=0;i<bases.get_rows();i++){
                    mu=0.0;
                    for(j=0;j<_chifn->get_dim();j++){
                        mu+=vv.get_data(j)*bases.get_data(i,j);
                    }
                    for(j=0;j<_chifn->get_dim();j++){
                        vv.subtract_val(j,mu*bases.get_data(i,j));
                    }
               }
               mu=vv.normalize();
                if(mu>1.0e-10){
                    bases.add_row(vv);
                }
            }


            for(i=0;i<_chifn->get_dim();i++){
                local_min.set(i,2.0*exception_value);
                local_max.set(i,-2.0*exception_value);
                for(ip=0;ip<_particles.get_dim();ip++){
                    mu=0.0;
                    i_half=trails.get_data(ip,trails.get_cols(ip)/2);
                    for(j=0;j<_chifn->get_dim();j++){
                       mu+=_chifn->get_pt(i_half,j)*bases.get_data(i,j);
                    }
                    if(mu<local_min.get_data(i)){
                        local_min.set(i,mu);
                    }
                    if(mu>local_max.get_data(i)){
                        local_max.set(i,mu);
                    }
                }
            }

            for(i=0;i<_chifn->get_dim();i++){
                norm.set(i,global_norm*(local_max.get_data(i)-local_min.get_data(i)));
            }

        }

        for(ip=0;ip<n_particles;ip++){
            _chifn->set_search_type(-1*(ip+2));

            i_dim=_chifn->random_int()%_chifn->get_dim();

            rr=normal_deviate(_chifn->get_dice(),0.0,1.0);

            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,_chifn->get_pt(_particles.get_data(ip),i)+
                            rr*norm.get_data(i_dim)*bases.get_data(i_dim,i));
            }
            mu=evaluate(trial,&i_found,ip,&mu_true);


            accept_it=0;


            if(ip>=_particles.get_dim() || mu<_chifn->get_fn(_particles.get_data(ip))){
                accept_it=1;
                needed_temp=1.0;
            }
            else{
                roll=_chifn->random_double();
                ratio=exp(-0.5*(mu-_chifn->get_fn(_particles.get_data(ip)))/_temp);
                needed_temp=-0.5*(mu-_chifn->get_fn(_particles.get_data(ip)))/log(roll);
                if(roll<ratio){
                    accept_it=1;
                }
            }

            needed_temp_dex.add(needed_temp_arr.get_dim());
            if(needed_temp<1.0){
                needed_temp_arr.add(1.0);
            }
            else{
                needed_temp_arr.add(needed_temp);
            }

            if(accept_it==1){
                _particles.set(ip,i_found);
                trails.add(ip,i_found);
                accepted.add_val(ip,1);
                total_accepted.add_val(ip,1);
            }
            else{
                trails.add(ip,_particles.get_data(ip));
            }

        }
        step_ct++;

        if(i_step>0 && i_step%adjust_every==0){

            accepted_sorted.reset_preserving_room();
            accepted_dex.reset_preserving_room();
            for(i=0;i<n_particles;i++){
                accepted_dex.set(i,i);
            }
            sort(accepted,accepted_sorted,accepted_dex);

            med_acc=accepted_sorted.get_data(accepted_dex.get_dim()/2);
            min_acc=accepted_sorted.get_data(0);
            max_acc=accepted_sorted.get_data(accepted_dex.get_dim()-1);

            has_been_adjusted=0;

            needs_adjustment=0;
            if(med_acc<step_ct/3){
                needs_adjustment=1;
            }
            else if(med_acc>(3*step_ct)/4){
                needs_adjustment=-1;
            }

            if(needs_adjustment!=0){
                sort(needed_temp_arr, needed_temp_sorted, needed_temp_dex);
                old_temp=_temp;
                _temp=needed_temp_sorted.get_data(needed_temp_dex.get_dim()/2);

                if(fabs(1.0-(_temp/old_temp))>0.01){
                    has_been_adjusted=1;
                }
            }

            for(ip=0;ip<_particles.get_dim();ip++){
                if(_since_min.get_data(ip)>=adjust_every){
                    min_pt_connected=0;
                    for(i=0;i<_chifn->get_dim();i++){
                        trial.set(i,0.5*(_chifn->get_pt(_chifn->mindex(),i)+
                                   _chifn->get_pt(_local_min.get_data(ip),i)));
                    }
                    mu=evaluate(trial,&i_found,ip,&mu_true);
                    if(mu<_chifn->get_fn(_local_min.get_data(ip)) || _local_min.get_data(ip)==_chifn->mindex()){
                        min_pt_connected=1;
                    }

                    if(min_pt_connected==0){
                        n_opt_out++;
                    }

                    if(min_pt_connected==1){
                        n_jumps++;

                        i_best=-1;
                        for(k=0;k<100 || i_best<0;k++){
                            for(i=0;i<_chifn->get_dim();i++){
                                dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
                            }
                            dir.normalize();
                            for(i=0;i<_chifn->get_dim();i++){
                                trial.set(i,_chifn->get_pt(_chifn->mindex(),i)+
                                          0.5*(_chifn->get_max(i)-_chifn->get_min(i))*dir.get_data(i));
                            }

                            for(i=0;i<_particles.get_dim();i++){
                                dd=0.0;
                                for(j=0;j<_chifn->get_dim();j++){
                                    dd+=power((trial.get_data(j)-_chifn->get_pt(_abs_min.get_data(i),j))/
                                               (_chifn->get_max(j)-_chifn->get_min(j)),2);
                                }
                                if(i==0 || dd<dd_min){
                                     dd_min=dd;
                                }
                            }

                            mu=evaluate(trial,&i_found,ip,&mu_true);

                            if(i_best<0 || dd_min>dd_best){
                                dd_best=dd_min;
                                i_best=i_found;
                            }

                        }

                        _particles.set(ip,i_best);
                        _local_min.set(ip,i_best);
                        trails.add(ip, i_best);
                        _since_min.set(ip,0);
                    }
                }
            }



            needed_temp_sorted.reset_preserving_room();

            for(ip=0;ip<_particles.get_dim();ip++){
                current_particles.set(ip,_particles.get_data(ip));
            }

            if(has_been_adjusted==1){
                needed_temp_arr.reset_preserving_room();
                needed_temp_dex.reset_preserving_room();
               for(i=0;i<n_particles;i++){
                   accepted.set(i,0);
                }
                step_ct=0;
            }

            if(has_been_adjusted==1){
                adjusted++;
            }
        }
    }

    _chifn->set_search_type(_type_init);

    array_1d<int> connected;
    double mu_min;
    int dex_min;

    for(i=0;i<_abs_min.get_dim();i++){
        if(i==0 || _chifn->get_fn(_abs_min.get_data(i))<mu_min){
             mu_min=_chifn->get_fn(_abs_min.get_data(i));
             dex_min=_abs_min.get_data(i);
        }
    }

    double min_disconnected=2.0*exception_value;
    int n_disconnected=0;

    for(i=0;i<_abs_min.get_dim();i++){
        if(_abs_min.get_data(i)==dex_min){
            connected.set(i,1);
        }
        else{
            for(j=0;j<_chifn->get_dim();j++){
                trial.set(j,0.5*(_chifn->get_pt(dex_min,j)+
                            _chifn->get_pt(_abs_min.get_data(i),j)));
            }
            mu=evaluate(trial,&j,-1,&mu_true);
            if(mu<_chifn->get_fn(_abs_min.get_data(i))){
                connected.set(i,1);
            }
            else{
                connected.set(i,0);
                n_disconnected++;
                if(_chifn->get_fn(_abs_min.get_data(i))<min_disconnected){
                    min_disconnected=_chifn->get_fn(_abs_min.get_data(i));
                }
            }
        }
    }

    printf("called %d -- %e\n",_chifn->get_pts(),_chifn->chimin());
    printf("min disconnected %e - %d\n",min_disconnected,n_disconnected);
    printf("jumped %d vs %d\n",n_jumps,n_opt_out);

    array_1d<double> smin,smax;
    for(i=0;i<_chifn->get_dim();i++){
        smin.set(i,0.0);
        smax.set(i,_chifn->get_characteristic_length(i));
    }
    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(smin,smax);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();

    min_vals.reset();
    min_val_sorted.reset();
    min_dexes.reset();
    for(i=0;i<_abs_min.get_dim();i++){
        min_vals.add(_chifn->get_fn(_abs_min.get_data(i)));
        min_dexes.add(_abs_min.get_data(i));
    }
    sort(min_vals, min_val_sorted, min_dexes);
    array_2d<double> seed;
    array_1d<int> chosen;
    for(i=0;i<_chifn->get_dim()+1;i++){
        seed.add_row(_chifn->get_pt(min_dexes.get_data(i)));
        chosen.add(min_dexes.get_data(i));
    }
    ffmin.find_minimum(seed,trial);

    seed.reset_preserving_room();
    min_vals.reset();
    min_val_sorted.reset();
    min_dexes.reset();
    for(i=0;i<_abs_min.get_dim();i++){
        if(connected.get_data(i)==0){
            min_vals.add(_chifn->get_fn(_abs_min.get_data(i)));
            min_dexes.add(_abs_min.get_data(i));
        }
    }
    sort(min_vals, min_val_sorted, min_dexes);
    for(i=0;i<min_dexes.get_dim() && seed.get_rows()!=_chifn->get_dim()+1;i++){
        seed.add_row(_chifn->get_pt(min_dexes.get_data(i)));
        chosen.add(min_dexes.get_data(i));
    }

    while(seed.get_rows()!=_chifn->get_dim()+1){
        i=_chifn->random_int()%_abs_min.get_dim();
        if(chosen.contains(_abs_min.get_data(i))==0 && _abs_min.get_data(i)!=dex_min){
            seed.add_row(_chifn->get_pt(_abs_min.get_data(i)));
            chosen.add(_abs_min.get_data(i));
        }
    }

    ffmin.find_minimum(seed,trial);

}
