#include "dalex.h"

void dalex::build(chisq_wrapper *cc){

    _chifn=cc;

    _min_explorers.set_chifn(_chifn);

    index_t i,j;
    _basis_vectors.set_cols(_chifn->get_dim());
    for(i=0;i<_chifn->get_dim();i++){
        _basis_norm.set(i,_chifn->get_characteristic_length(i));
        for(j=0;j<_chifn->get_dim();j++){
            if(i==j){
                _basis_vectors.set(i,j,1.0);
            }
            else{
                _basis_vectors.set(i,j,0.0);
            }
        }
    }

}

value_t dalex::get_norm(index_t dex){

    if(_good_points.get_dim()==0){
        return _chifn->get_characteristic_length(dex);
    }

    value_t min,max;
    index_t i,ip;
    for(i=0;i<_good_points.get_dim();i++){

        ip=_good_points.get_data(i);

        if(i==0 || _chifn->get_pt(ip,dex)<min){
            min=_chifn->get_pt(ip,dex);
        }

        if(i==0 || _chifn->get_pt(ip,dex)>max){
            max=_chifn->get_pt(ip,dex);
        }
    }

    if(max-min<1.0e-20){
        return _chifn->get_characteristic_length(dex);
    }

    return max-min;
}

void dalex::search(){
    safety_check("search");
    index_t i,j,i_found;
    value_t mu;
    index_t pts_0=_chifn->get_pts();
    array_1d<value_t> pt;
    array_1d<index_t> to_use,explorer_dex;
    assess_good_points();

    index_t has_explored=0;

    while(_chifn->get_pts()<_limit){
        _chifn->set_search_type(_type_refine);
        iterate_on_minimum();

        _chifn->set_search_type(_type_tendril);
        tendril_search();
        _update_good_points();
    }

}


void dalex::simplex_search(){
    array_1d<index_t> empty;
    simplex_search(empty);
}

void dalex::simplex_search(index_t ii){
    array_1d<index_t> specified;
    specified.add(ii);
    simplex_search(specified);
}

void dalex::simplex_search(array_1d<index_t> &specified){

    // printf("\n    doing dalex_simplex %e -- %e %e\n",chimin(),target(),_target_factor);
    safety_check("simplex_search");
    array_1d<value_t> min,max;
    min.set_name("dalex_simplex_min");
    max.set_name("dalex_simplex_max");
    array_2d<value_t> seed;
    seed.set_name("dalex_simplex_seed");
    array_1d<value_t> trial,grad;

    index_t i,j;
    index_t pt_0=_chifn->get_pts();

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    for(i=0;i<specified.get_dim();i++){
        seed.add_row(_chifn->get_pt(specified.get_data(i)));
    }

    array_1d<value_t> mu_arr,mu_arr_sorted;
    array_1d<index_t> mu_dex;
    for(i=0;i<_min_explorers.get_n_particles();i++){
        mu_dex.set(i,i);
        mu_arr.set(i,_min_explorers.get_mu(i));
    }
    sort(mu_arr, mu_arr_sorted,mu_dex);

    for(i=0;i<_chifn->get_dim();i++){
        j=mu_dex.get_data(i);
        _min_explorers.get_pt(j,trial);
            seed.add_row(trial);
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(_chifn);
    ffmin.set_minmax(min,max);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.use_gradient();
    ffmin.find_minimum(seed,trial);

    // printf("    after dalex_simplex chimin %e\n",chimin());
    _simplex_mindex=mindex();

    array_1d<value_t> min_pt;
    ffmin.get_minpt(min_pt);
    value_t mu;
    index_t i_found;
    evaluate(min_pt, &mu, &i_found);

    _update_good_points(pt_0);

}


index_t dalex::bisection(index_t ilow, const array_1d<value_t> &dir, value_t local_target, value_t tol){
    safety_check("bisection(index_t, arr)");
    array_1d<value_t> trial_high;
    trial_high.set_name("dalex_bisection_trial_high");
    value_t rr=1.0;
    index_t ii,i_found;
    value_t mu=-2.0*exception_value;
    for(ii=0;ii<_chifn->get_dim();ii++){
        trial_high.set(ii,_chifn->get_pt(ilow,ii));
    }

    while(mu<=local_target){
        for(ii=0;ii<_chifn->get_dim();ii++){
            trial_high.add_val(ii,rr*dir.get_data(ii));
        }
        evaluate(trial_high, &mu, &i_found);
        rr*=2.0;
    }

    return bisection(_chifn->get_pt(ilow), trial_high, local_target, tol);
}


index_t dalex::bisection(index_t ilow, index_t ihigh, value_t local_target, value_t tol){
    safety_check("bisection(index_t, index_t)");
    return bisection(_chifn->get_pt(ilow), _chifn->get_pt(ihigh),
                     local_target, tol);
}


index_t dalex::bisection(const array_1d<value_t>& lowball_in, const array_1d<value_t>& highball_in,
                     value_t local_target, value_t tol){

    safety_check("bisection(arr, arr)");
    array_1d<value_t> trial,lowball,highball;
    trial.set_name("dalex_bisection_trial");
    lowball.set_name("dalex_bisection_lowball");
    highball.set_name("dalex_bisection_highball");
    index_t i_found,ii,jj,i_trial;

    value_t mu,flow,fhigh;
    evaluate(lowball_in, &flow, &ii);
    evaluate(highball_in, &fhigh, &jj);

    if(flow>local_target){
        // printf("WARNING flow is greater than target %e %e\n",flow,local_target);
        // printf("min %e target %e\n",chimin(),target());
        exit(1);
    }

    i_found=ii;

    for(ii=0;ii<_chifn->get_dim();ii++){
        lowball.set(ii,lowball_in.get_data(ii));
        highball.set(ii,highball_in.get_data(ii));
    }

    index_t ct;
    value_t wgt_low;
    for(ct=0;ct<20 &&
       (ct<5 || fabs(_chifn->get_fn(i_found)-local_target)>tol); ct++){

        wgt_low=(fhigh-local_target)/(fhigh-flow);
        if(wgt_low<0.25 || wgt_low>0.75){
            wgt_low=0.5;
        }

        for(ii=0;ii<_chifn->get_dim();ii++){
            trial.set(ii, wgt_low*lowball.get_data(ii)+(1.0-wgt_low)*highball.get_data(ii));
        }

        evaluate(trial,&mu,&i_trial);

        if(mu<local_target){
            flow=mu;
            for(ii=0;ii<_chifn->get_dim();ii++){
                lowball.set(ii,trial.get_data(ii));
            }
        }
        else{
            fhigh=mu;
            for(ii=0;ii<_chifn->get_dim();ii++){
                highball.set(ii,trial.get_data(ii));
            }
        }

        if(fabs(mu-local_target)<fabs(_chifn->get_fn(i_found)-local_target)){
            i_found=i_trial;
        }
    }

    return i_found;

}


void dalex::find_trial_bases(index_t idim, array_1d<value_t> &dx, array_2d<value_t> &bases_out){

    index_t i,j;

    array_1d<value_t> trial_dir;
    trial_dir.set_name("dalex_find_trial_bases_trial_dir");

    if(_basis_vectors.get_rows()==0){
        _basis_vectors.reset();
        _basis_vectors.set_cols(_chifn->get_dim());
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                _basis_vectors.set(i,j,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
        }
    }

    if(dx.get_dim()==0){
        for(i=0;i<_chifn->get_dim();i++){
            dx.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dx.normalize();
    }

    bases_out.set_cols(_basis_vectors.get_cols());
    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            bases_out.set(i,j,_basis_vectors.get_data(i,j));
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        bases_out.add_val(idim,i,dx.get_data(i));
    }
    bases_out.normalize_row(idim);

    index_t ix,jx;
    value_t mu;

    for(ix=idim+1;ix!=idim;){
        if(ix>=_chifn->get_dim()){
            ix=0;
        }

        for(jx=idim;jx!=ix;){
            if(ix>=_chifn->get_dim()){
                jx=0;
            }

            mu=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                mu+=bases_out.get_data(ix,i)*bases_out.get_data(jx,i);
            }
            for(i=0;i<_chifn->get_dim();i++){
                bases_out.subtract_val(ix,i,mu*bases_out.get_data(jx,i));
            }

            if(jx<_chifn->get_dim()-1)jx++;
            else jx=0;
        }

        bases_out.normalize_row(ix);

        if(ix<_chifn->get_dim()-1)ix++;
        else ix=0;

    }

    validate_bases(bases_out,"dalex_find_trial_bases");


}

void dalex::validate_bases(array_2d<value_t> &bases, char *whereami){
    index_t ix,i,jx;
    value_t mu;
    /////////////////testing
    for(ix=0;ix<_chifn->get_dim();ix++){
        bases.normalize_row(ix);
        mu=0.0;
        for(i=0;i<_chifn->get_dim();i++){
            mu+=bases.get_data(ix,i)*bases.get_data(ix,i);
        }
        if(fabs(mu-1.0)>1.0e-6){
            // printf("WARNING in %s, square norm %e\n",whereami,mu);
            exit(1);
        }

        for(jx=ix+1;jx<_chifn->get_dim();jx++){
            mu=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                mu+=bases.get_data(ix,i)*bases.get_data(jx,i);
            }

            if(fabs(mu)>1.0e-6){
                // printf("WARNING in %s, dot product %e\n",whereami,mu);
                exit(1);
            }
        }
    }

}

void dalex::guess_bases(array_2d<value_t> &bases){
    safety_check("guess_bases");

    array_2d<value_t> covar;
    covar.set_name("dalex_guess_bases_covar");
    covar.set_cols(_chifn->get_dim());
    bases.set_cols(_chifn->get_dim());
    index_t ix,iy;
    value_t covarmax=-1.0;

    find_covariance_matrix(mindex(),covar);

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=ix;iy<_chifn->get_dim();iy++){
            if(fabs(covar.get_data(ix,iy))>covarmax){
                covarmax=fabs(covar.get_data(ix,iy));
            }
        }
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=0;iy<_chifn->get_dim();iy++){
            covar.divide_val(ix,iy,covarmax);
        }
    }

    // printf("assembled covariance matrix\n");

    array_2d<value_t> evecs;
    evecs.set_name("dalex_guess_bases_evecs");
    array_1d<value_t> evals;
    evals.set_name("dalex_guess_bases_evals");

    evecs.set_cols(_chifn->get_dim());

    try{
        eval_symm(covar,evecs,evals,0.1);
    }
    catch(index_t iex){
        // printf("Guess failed on of eigen vectors\n");
        throw -1;
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=0;iy<_chifn->get_dim();iy++){
            bases.set(ix,iy,evecs.get_data(ix,iy));
        }
        bases.normalize_row(ix);
    }

    // printf("validated guessed bases\n");
}



void dalex::find_bases(){

    safety_check("find_bases");

    _minimizers.reset_preserving_room();

    array_1d<value_t> dir;
    index_t i,j;
    for(i=0;i<_basis_associates.get_dim();i++){
        if(_chifn->get_fn(_basis_associates.get_data(i))>target()){
            _basis_associates.remove(i);
            i--;
        }
    }
    index_t i_pt;
    value_t mu;
    index_t n_0=_basis_associates.get_dim();
    index_t old_search_type = _chifn->get_search_type();
    _chifn->set_search_type(_type_find_bases);
    while(_basis_associates.get_dim()<4*_chifn->get_dim()+n_0){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        dir.normalize();
        i_pt=bisection(mindex(),dir,0.5*(target()+chimin()),0.1);

        /*if(fabs(_chifn->get_fn(i_pt)-0.5*(target()+chimin()))>1.0 && chimin()<500.0){
            // printf("WARNING failed to get associate within tol %e %e-- %e %e\n",
            0.5*(target()+chimin()),_chifn->get_fn(i_pt),chimin(),target());
            // printf("first\n");
            exit(1);
        }*/

        if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
            _basis_associates.add(i_pt);
            i_pt=bisection(mindex(),dir,0.25*chimin()+0.75*target(),0.1);
            if(i_pt!=mindex() && _basis_associates.contains(i_pt)==0){
                _basis_associates.add(i_pt);

                /*if(fabs(_chifn->get_fn(i_pt)-(0.75*target()+0.25*chimin()))>1.0 && chimin()<500.0){
                    // printf("WARNING failed to get associate within tol %e %e -- %e %e\n",
                    0.75*target()+0.25*chimin(),_chifn->get_fn(i_pt),chimin(),target());
                    // printf("second\n");
                    exit(1);
                }*/

            }
        }
    }
    _chifn->set_search_type(old_search_type);

    array_2d<value_t> trial_bases;
    array_1d<value_t> trial_model,dx;

    trial_bases.set_name("dalex_find_bases_trial_bases");
    trial_model.set_name("dalex_find_bases_trial_model");
    dx.set_name("dalex_find_bases_dx");

    if(_basis_vectors.get_rows()==0){
        find_trial_bases(0,dx,trial_bases);
        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                _basis_vectors.set(i,j,trial_bases.get_data(i,j));
            }
        }
    }

    index_t dimsq=_chifn->get_dim()*_chifn->get_dim();

    if(_basis_associates.get_dim()==0){
        // printf("WARNING _basis associates is empty\n");
        exit(1);
    }


    index_t ct,idim,aborted,max_abort,total_aborted,changed_bases;
    value_t error0,error,errorBest,stdev,stdevlim,error1;

    stdev=0.1/sqrt(value_t(_chifn->get_dim()));
    stdevlim=1.0e-5/sqrt(value_t(_chifn->get_dim()));
    max_abort=_chifn->get_dim()*100;

    error=basis_error(_basis_vectors,_basis_model);
    error0=error;
    error1=error;
    errorBest=error;
    aborted=0;
    total_aborted=0;
    changed_bases=0;
    ct=0;

    // printf("center %e min %e\n",_chifn->get_fn(mindex()),chimin());
    // printf("error0 %e\n", error0);

    if(1==1){
        // printf("guessing basis from second derivative\n");
        try{
            guess_bases(trial_bases);
            error=basis_error(trial_bases,trial_model);
            // printf("guess got error %e\n",error);
            if(error<error0){
                changed_bases=1;
                for(i=0;i<_chifn->get_dim();i++){
                    _basis_model.set(i,trial_model.get_data(i));
                    for(j=0;j<_chifn->get_dim();j++){
                        _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                    }
                }
                error1=error;
                errorBest=error;
            }
        }
        catch(index_t iex){
            // printf("never mind; guess did not work\n");
        }
    }

    value_t penalty=0.0;
    for(i=0;i<_chifn->get_dim();i++){
        if(_basis_model.get_data(i)<0.0){
            penalty+=1.0;
        }
    }

    value_t pp0=penalty;

    while(stdev>stdevlim && aborted<max_abort){
        ct++;
        idim=-1;
        while(idim>=_chifn->get_dim() || idim<0){
            idim=_chifn->random_int()%_chifn->get_dim();
        }

        for(i=0;i<_chifn->get_dim();i++){
            dx.set(i,normal_deviate(_chifn->get_dice(),0.0,stdev));
        }

        find_trial_bases(idim,dx,trial_bases);
        error=basis_error(trial_bases,trial_model);

        if(error<errorBest){
            penalty=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                if(trial_model.get_data(i)<0.0){
                    penalty+=1.0;
                }
            }
            if(error1-error>1.0e-5*error){
                aborted=0;
            }
            else{
                aborted++;
                total_aborted++;
            }

            changed_bases=1;
            for(i=0;i<_chifn->get_dim();i++){
                _basis_model.set(i,trial_model.get_data(i));
                for(j=0;j<_chifn->get_dim();j++){
                    _basis_vectors.set(i,j,trial_bases.get_data(i,j));
                }
            }
            errorBest=error;
        }
        else{
            aborted++;
            total_aborted++;
        }

        if(ct%(max_abort/2)==0){
            if(total_aborted<(3*ct)/4)stdev*=1.5;
            else if(total_aborted>(3*ct)/4)stdev*=0.5;
        }

        if(ct%1000==0){
            error1=errorBest;
            // printf("    ct %d error %.2e pp %.1e from %.2e %.1e min %.3e pts %d\n",
            // ct,errorBest,penalty,error0,pp0,chimin(),_basis_associates.get_dim());
        }
    }

    index_t k;
    array_1d<value_t> min,max,focused_min,focused_max;
    min.set_name("find_bases_min");
    max.set_name("find_bases_max");
    focused_min.set_name("find_bases_focused_min");
    focused_max.set_name("find_bases_focused_max");
    if(changed_bases==1){
        for(i=0;i<_good_points.get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                mu=0.0;
                for(k=0;k<_chifn->get_dim();k++){
                    mu+=_chifn->get_pt(_good_points.get_data(i),k)*_basis_vectors.get_data(j,k);
                }
                if(i==0 || mu<min.get_data(j)){
                    min.set(j,mu);
                }
                if(i==0 || mu>max.get_data(j)){
                    max.set(j,mu);
                }

                if(_basis_associates.contains(_good_points.get_data(i))==1){
                    if(j>=focused_min.get_dim() || mu<focused_min.get_data(j)){
                        focused_min.set(j,mu);
                    }
                    if(j>=focused_max.get_dim() || mu>focused_max.get_data(j)){
                        focused_max.set(j,mu);
                    }
                }
            }
        }
        for(i=0;i<_chifn->get_dim();i++){
            _basis_norm.set(i,max.get_data(i)-min.get_data(i));
            if(_basis_norm.get_data(i)<1.0e-20){
                _basis_norm.set(i,1.0);
            }

            _basis_associate_norm.set(i,focused_max.get_data(i)-focused_min.get_data(i));
            if(_basis_associate_norm.get_data(i)<1.0e-20){
                // printf("WARNING basis associate norm %d %e\n",
                // i,_basis_associate_norm.get_data(i));
            }
        }

    }

    _basis_chimin=chimin();
    // printf("done finding bases\n");


}


value_t dalex::basis_error(array_2d<value_t> &trial_bases, array_1d<value_t> &trial_model){

    if(_basis_associates.get_dim()<=0){
        // printf("WARNING cannot calculate basis error there are only %d associates\n",
        // _basis_associates.get_dim());

        exit(1);
    }

    safety_check("basis_error");

    trial_model.zero();
    if(_basis_ddsq.get_rows()>_basis_associates.get_dim()){
        _basis_ddsq.reset();
    }

    if(_basis_ddsq.get_cols()!=_chifn->get_dim()){
        _basis_ddsq.set_cols(_chifn->get_dim());
    }

    if(_basis_mm.get_dim()!=_chifn->get_dim()*_chifn->get_dim()){
        _basis_mm.set_dim(_chifn->get_dim()*_chifn->get_dim());
    }

    if(_basis_bb.get_dim()!=_chifn->get_dim()){
        _basis_bb.set_dim(_chifn->get_dim());
    }

    if(_basis_vv.get_dim()!=_chifn->get_dim()){
        _basis_vv.set_dim(_chifn->get_dim());
    }

    _basis_mm.zero();
    _basis_bb.zero();
    _basis_vv.zero();
    _basis_ddsq.zero();

    index_t i,j,ix;
    value_t mu;
    for(ix=0;ix<_basis_associates.get_dim();ix++){
        for(i=0;i<_chifn->get_dim();i++){
            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=(_chifn->get_pt(_basis_associates.get_data(ix),j)-_chifn->get_pt(mindex(),j))*trial_bases.get_data(i,j);
            }
            _basis_ddsq.set(ix,i,mu*mu);
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        for(j=0;j<_basis_associates.get_dim();j++){
            _basis_bb.add_val(i,_basis_ddsq.get_data(j,i)*(_chifn->get_fn(_basis_associates.get_data(j))-chimin()));
        }
    }

    index_t k;
    for(i=0;i<_chifn->get_dim();i++){
        for(j=i;j<_chifn->get_dim();j++){
            ix=i*_chifn->get_dim()+j;
            for(k=0;k<_basis_associates.get_dim();k++){
                _basis_mm.add_val(ix,_basis_ddsq.get_data(k,i)*_basis_ddsq.get_data(k,j));
            }
            if(j!=i){
                _basis_mm.set(j*_chifn->get_dim()+i,_basis_mm.get_data(ix));
            }
        }
    }

    try{
        naive_gaussian_solver(_basis_mm,_basis_bb,trial_model,_chifn->get_dim());
    }
    catch(index_t iex){
        // printf("WARNING basis_error was no good\n");
        return 2.0*exception_value;
    }

    value_t error=0.0,chi_model;
    for(i=0;i<_basis_associates.get_dim();i++){
        chi_model=chimin();
        for(j=0;j<_chifn->get_dim();j++){
            chi_model+=trial_model.get_data(j)*_basis_ddsq.get_data(i,j);
        }
        error+=power(_chifn->get_fn(_basis_associates.get_data(i))-chi_model,2);
    }

    return error/value_t(_basis_associates.get_dim());

}


void dalex::find_covariance_matrix(index_t iCenter, array_2d<value_t> &covar){
    safety_check("find_covariance_matrix");


    array_1d<value_t> trial,center,norm;
    trial.set_name("dalex_findCovar_trial");
    norm.set_name("dalex_findCovar_norm");
    center.set_name("dalex_findCovar_center");

    value_t fcenter;
    index_t i;

    fcenter=_chifn->get_fn(iCenter);

    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,_chifn->get_pt(iCenter,i));
        norm.set(i,1.0);
    }

    array_2d<value_t> fpp,fpm,fmp,fmm;
    fpp.set_name("dalex_findCovar_fpp");
    fpm.set_name("dalex_findCovar_fpm");
    fmp.set_name("dalex_findCovar_fmp");
    fmm.set_name("dalex_findCovar_fmm");

    fpp.set_cols(_chifn->get_dim());
    fpm.set_cols(_chifn->get_dim());
    fmp.set_cols(_chifn->get_dim());
    fmm.set_cols(_chifn->get_dim());

    array_1d<value_t> f2p,f2m;
    f2p.set_name("dalex_findCovar_f2p");
    f2m.set_name("dalex_findCovar_f2m");

    index_t ifpp,ifpm,ifmp,ifmm,if2p,if2m;

    value_t mu;
    array_1d<value_t> dx;
    dx.set_name("dalex_findCovar_dx");
    for(i=0;i<_chifn->get_dim();i++){
        dx.set(i,1.0e-2);
    }

    index_t ix,iy,keepGoing,ctAbort,ctAbortMax,calledMax;
    index_t ibefore=_chifn->get_called();

    ctAbort=0;
    ctAbortMax=100;
    calledMax = 10*_chifn->get_dim()*_chifn->get_dim();

    for(i=0;i<_chifn->get_dim();i++){
        trial.set(i,center.get_data(i));
    }

    for(ix=0;ix<_chifn->get_dim();ix++){
        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,center.get_data(i));
        }

        if(_chifn->get_called()-ibefore>calledMax ||
           ctAbort>=ctAbortMax){
                // printf("Could not find CoVar; aborting\n");
                // printf("ctAbort %d\n",ctAbort);
                // printf("called %d\n",_chifn->get_called()-ibefore);
                throw -1;
        }

        if(iCenter<0){
            // printf("Center is invalid; aborting\n");
            throw -1;
        }

        keepGoing=1;
        trial.set(ix,center.get_data(ix)+2.0*dx.get_data(ix)*norm.get_data(ix));
        evaluate(trial,&mu,&if2p);

        if(if2p>=0 && if2p!=iCenter){
            f2p.set(ix,mu);
        }
        else if(if2p==iCenter){
            dx.multiply_val(ix,1.5);
            keepGoing=0;
            ix--;
            ctAbort++;
        }
        else if(if2p<0){
            center.subtract_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
            evaluate(center,&fcenter,&iCenter);
            keepGoing=0;
            ix--;
            ctAbort++;
        }

        if(keepGoing==1){
            trial.set(ix,center.get_data(ix)-2.0*dx.get_data(ix)*norm.get_data(ix));
            evaluate(trial,&mu,&if2m);

            if(if2m>=0 && if2m!=iCenter){
                f2m.set(ix,mu);
            }
            else if(if2m==iCenter){
                dx.multiply_val(ix,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(if2m<0){
                center.add_val(ix,2.5*dx.get_data(ix)*norm.get_data(ix));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
        }

        for(iy=ix-1;iy>=0 && keepGoing==1;iy--){
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,center.get_data(i));
            }

            if(_chifn->get_called()-ibefore>calledMax ||
               ctAbort>=ctAbortMax){
                // printf("Could not find CoVar; aborting\n");
                // printf("ctAbort %d\n",ctAbort);
                // printf("called %d\n",_chifn->get_called()-ibefore);
                throw -1;
            }

            if(iCenter<0){
                // printf("center is invalid; aborting\n");
                throw -1;
            }

            trial.set(ix,center.get_data(ix)+dx.get_data(ix)*norm.get_data(ix));
            trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
            evaluate(trial,&mu,&ifpp);
            if(ifpp>=0 && ifpp!=iCenter){
                fpp.set(ix,iy,mu);
            }
            else if(ifpp==iCenter){
                dx.multiply_val(ix,1.5);
                dx.multiply_val(iy,1.5);
                keepGoing=0;
                ix--;
                ctAbort++;
            }
            else if(ifpp<0){
                center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                evaluate(center,&fcenter,&iCenter);
                keepGoing=0;
                ix--;
                ctAbort++;
            }

            if(keepGoing==1){
               trial.set(iy,center.get_data(iy)-dx.get_data(iy)*norm.get_data(iy));
               evaluate(trial,&mu,&ifpm);
               if(ifpm>=0 && ifpm!=iCenter){
                   fpm.set(ix,iy,mu);
               }
               else if(ifpm==iCenter){
                   dx.multiply_val(ix,1.5);
                   dx.multiply_val(iy,1.5);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
               else if(ifpm<0){
                   center.subtract_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                   center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                   evaluate(center,&fcenter,&iCenter);
                   keepGoing=0;
                   ix--;
                   ctAbort++;
               }
            }

            if(keepGoing==1){
                trial.set(ix,center.get_data(ix)-dx.get_data(ix)*norm.get_data(ix));
                evaluate(trial,&mu,&ifmm);
                if(ifmm>=0 && ifmm!=iCenter){
                    fmm.set(ix,iy,mu);
                }
                else if(ifmm==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmm<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.add_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

            if(keepGoing==1){
                trial.set(iy,center.get_data(iy)+dx.get_data(iy)*norm.get_data(iy));
                evaluate(trial,&mu,&ifmp);
                if(ifmp>=0 && ifmp!=iCenter){
                    fmp.set(ix,iy,mu);
                }
                else if(ifmp==iCenter){
                    dx.multiply_val(ix,1.5);
                    dx.multiply_val(iy,1.5);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
                else if(ifmp<0){
                    center.add_val(ix,1.5*dx.get_data(ix)*norm.get_data(ix));
                    center.subtract_val(iy,1.5*dx.get_data(iy)*norm.get_data(iy));
                    evaluate(center,&fcenter,&iCenter);
                    keepGoing=0;
                    ix--;
                    ctAbort++;
                }
            }

        }
    }

    covar.set_cols(_chifn->get_dim());
    for(ix=0;ix<_chifn->get_dim();ix++){
        covar.set(ix,ix,0.25*(f2p.get_data(ix)+f2m.get_data(ix)-2.0*fcenter)/power(dx.get_data(ix)*norm.get_data(ix),2));
    }

    value_t num,denom;
    for(ix=0;ix<_chifn->get_dim();ix++){
        for(iy=ix-1;iy>=0;iy--){
            num=0.25*(fpp.get_data(ix,iy)+fmm.get_data(ix,iy)-fmp.get_data(ix,iy)-fpm.get_data(ix,iy));
            denom=dx.get_data(ix)*norm.get_data(ix)*dx.get_data(iy)*norm.get_data(iy);
            covar.set(ix,iy,num/denom);
            covar.set(iy,ix,num/denom);
        }
    }


}

void dalex::get_negative_gradient(index_t i_origin, cost_fn &cost, ellipse &dummy_ellipse, array_1d<value_t> &out_dir){
    value_t mu_pos, mu_neg;
    array_1d<value_t> trial;
    trial.set_name("get_neg_grad_trial");
    index_t i,j;
    out_dir.reset_preserving_room();
    for(i=0;i<_chifn->get_dim();i++){
        out_dir.set(i,0.0);
    }
    value_t rat=0.01;
    value_t delta;
    for(i=0;i<_chifn->get_dim();i++){
        delta=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            delta+=(_chifn->get_pt(i_origin,j)-_chifn->get_pt(mindex(),j))*dummy_ellipse.bases(i,j);
        }
        delta*=rat;

        for(j=0;j<_chifn->get_dim();j++){
            trial.set(j,_chifn->get_pt(i_origin,j)+rat*dummy_ellipse.radii(i)*dummy_ellipse.bases(i,j));
        }
        mu_pos=cost(trial);
        for(j=0;j<_chifn->get_dim();j++){
            trial.set(j,_chifn->get_pt(i_origin,j)-rat*dummy_ellipse.radii(i)*dummy_ellipse.bases(i,j));
        }
        mu_neg=cost(trial);
        for(j=0;j<_chifn->get_dim();j++){
            out_dir.add_val(j,(mu_neg-mu_pos)*dummy_ellipse.bases(i,j)/(2.0*rat*dummy_ellipse.radii(i)));
        }
    }
}

index_t dalex::simplex_boundary_search(const index_t specified, const index_t i_origin,
                                   ellipse_list &exclusion_zones, index_t *i_next){

    safety_check("simplex_boundary_search");
    // printf("\ndoing dalex.simplex_boundary_search() %d\n",_chifn->get_pts());
    index_t pt_start=_chifn->get_pts();
    assess_good_points();
    // printf("specified %d origin %d\n",specified,i_origin);

    if(_chifn->get_dim()>9 && specified>=0){
        // printf("    starting from %e -- %e %e\n",
        // _chifn->get_fn(specified),_chifn->get_pt(specified,6), _chifn->get_pt(specified,9));
    }

    index_t i_node,i_pt;
    index_t i,j,k;
    value_t xmin,xmax,xx;

    index_t n_good_0=_good_points.get_dim();

    array_1d<index_t> associates;
    associates.set_name("dalex_simplex_boundary_associates");
    index_t i_start;
    index_t n_thin=-1;
    index_t thin_ct=0;
    index_t ip,io;

    value_t mu;
    array_1d<value_t> trial;
    trial.set_name("simplex_boundary_trial");

    array_1d<index_t> end_points_checked;
    array_1d<value_t> end_points_chi;
    end_points_checked.set_name("end_points_checked");
    end_points_chi.set_name("end_points_chi");

    index_t accept_it;
    index_t associate_iterations=0;
    index_t thin_to=20000;

    accept_it=0;
    while(accept_it==0){
        associate_iterations++;
        if(associate_iterations>2){
            // printf("WARNING associate_iterations %d\n",associate_iterations);
            exit(1);
        }
        for(i=0;i<_tendril_path.get_rows();i++){
            ip=_tendril_path.get_data(i,0);
            io=_tendril_path.get_data(i,1);
            if(end_points_checked.contains(io)==0){
                for(j=0;j<_chifn->get_dim();j++){
                    trial.set(j,0.5*(_chifn->get_pt(specified,j)+_chifn->get_pt(io,j)));
                }

                evaluate(trial,&mu,&k);
                end_points_checked.add(io);
                end_points_chi.add(mu);
            }
            else{
                for(j=0;j<end_points_checked.get_dim();j++){
                    if(end_points_checked.get_data(j)==io){
                        mu=end_points_chi.get_data(j);
                        break;
                    }
                }
            }

            if(mu<target() || specified<=0){
                if(associates.contains(io)==0){
                    associates.add(io);
                }
            }

            if(associates.contains(ip)==0){
                if(specified<=0 || _has_struck==0 || mu<target()){
                    thin_ct++;
                    if(n_thin<0 || thin_ct%n_thin==0){
                        associates.add(ip);
                    }
                }
            }
        }

        if(associates.get_dim()-end_points_checked.get_dim()<2*thin_to){
            accept_it=1;
        }
        else{
            n_thin=associates.get_dim()/thin_to;
            if(n_thin==1){
                n_thin=2;
            }
            thin_ct=0;
            associates.reset();
        }
    }

    index_t associates_from_tendril=1;
    if(associates.get_dim()==0){
        associates_from_tendril=0;

        if(_good_points.get_dim()>thin_to){
            n_thin=_good_points.get_dim()/thin_to;
            if(n_thin==1){
                n_thin=2;
            }
        }
        for(i=0;i<_good_points.get_dim();i++){
            if(i%n_thin==0){
                associates.add(_good_points.get_data(i));
            }
        }
    }

    cost_fn dchifn(_chifn,associates);

    // printf("    associates %d path %d\n", associates.get_dim(),_tendril_path.get_rows());

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    array_1d<value_t> min,max;
    min.set_name("dalex_simplex_search_min");
    max.set_name("dalex_simplex_search_min");

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    ffmin.set_minmax(min,max);

    index_t old_type=_chifn->get_search_type();
    _chifn->set_search_type(_type_tendril_seed);

    array_2d<value_t> seed;
    seed.set_name("dalex_simplex_search_seed");
    array_1d<value_t> trial1,trial2;
    trial1.set_name("dalex_simplex_trial1");
    trial2.set_name("dalex_simplex_trial2");

    index_t d_step = _strikes+1;

    array_2d<value_t> ellipse_pts;
    ellipse_pts.set_name("dalex_simplex_boundary_ellipse_pts");
    ellipse dummy_ellipse;
    if(i_origin>=0 && associates_from_tendril==1){
        for(i=0;i<associates.get_dim();i++){
            ellipse_pts.add_row(_chifn->get_pt(associates.get_data(i)));
        }
    }

    array_1d<value_t> distances,distances_sorted;
    array_1d<index_t> distances_dex;
    distances.set_name("simp_bou_dist");
    distances_sorted.set_name("simp_bou_dist_sorted");
    distances_dex.set_name("simp_bou_dist_dex");
    if(ellipse_pts.get_rows()<2*_chifn->get_dim()){
        // printf("    need to build ellipse out of nearest good points\n");
        for(i=0;i<_good_points.get_dim();i++){
            distances.add(cardinal_distance(specified,_good_points.get_data(i)));
            distances_dex.add(_good_points.get_data(i));
        }
        sort(distances,distances_sorted,distances_dex);
        ellipse_pts.reset_preserving_room();
        for(i=0;i<2*_chifn->get_dim();i++){
            ellipse_pts.add_row(_chifn->get_pt(distances_dex.get_data(i)));
        }
    }

    dummy_ellipse.build(ellipse_pts);

    index_t i_anchor;
    array_1d<value_t> base_dir;
    base_dir.set_name("dalex_simplex_boundary_base_dir");
    if(i_origin>=0){
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(specified,i)-_chifn->get_pt(i_origin,i));
        }
    }
    else{
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(specified,i)-_chifn->get_pt(mindex(),i));
        }
        if(_chifn->get_fn(specified)>target()){
            for(i=0;i<_chifn->get_dim();i++){
                base_dir.multiply_val(i,-1.0);
            }
        }
    }
    base_dir.normalize();

    i_anchor=specified;

    array_1d<value_t> anchor;
    anchor.set_name("simplex_boundary_anchor");
    for(i=0;i<_chifn->get_dim();i++){
        anchor.set(i,_chifn->get_pt(i_anchor,i));
    }
    array_1d<value_t> bisect_dir,epsilon;
    bisect_dir.set_name("simplex_boundary_bisect_dir");
    epsilon.set_name("simplex_boundary_bisect_dir");
    index_t i_bisect1,i_bisect2,i_chosen;
    value_t mu1,mu2;
    value_t component,rat;
    array_1d<index_t> kept_dex;
    kept_dex.set_name("simplex_boundary_kept_dex");
    array_1d<value_t> avg_pt;
    avg_pt.set_name("simplex_boundary_avg_pt");

    value_t local_target,tt;
    value_t dd,avg_dd;

    if(specified>=0){
        i_anchor=specified;
        if(_chifn->get_fn(i_anchor)<target()){
            local_target=target();
        }
        else{
            local_target=_chifn->get_fn(i_anchor)+0.5*(target()-chimin());
        }
        while(seed.get_rows()!=_chifn->get_dim()+1){
            for(i=0;i<_chifn->get_dim();i++){
                for(j=0;j<_chifn->get_dim();j++){
                    bisect_dir.set(j,dummy_ellipse.bases(i,j)+base_dir.get_data(j));
                }
                i_bisect1=i_anchor;
                tt=local_target;
                while(i_bisect1==i_anchor){
                    i_bisect1=bisection(i_anchor,bisect_dir,tt,0.001);
                    tt+=0.5*(target()-chimin());
                }
                for(j=0;j<_chifn->get_dim();j++){
                    bisect_dir.set(j,-1.0*dummy_ellipse.bases(i,j)+base_dir.get_data(j));
                }

                // note the bug below; i_bisect2 does not get used
                // because it will never be set to i_bisect2 != i_anchor
                // this seems to have not hurt performance
                i_bisect2=i_anchor;
                tt=local_target;
                while(i_bisect1==i_anchor){
                    i_bisect2=bisection(i_anchor,bisect_dir,tt,0.001);
                    tt+=0.5*(target()-chimin());
                }
                mu1=dchifn(_chifn->get_pt(i_bisect1));
                mu2=dchifn(_chifn->get_pt(i_bisect2));

                i_chosen=-1;
                if(i_bisect1==i_anchor && i_bisect2!=i_anchor){
                    i_chosen=i_bisect2;
                }
                else if(i_bisect1!=i_anchor && i_bisect2==i_anchor){
                    i_chosen=i_bisect1;
                }
                else if(i_bisect1!=i_anchor && i_bisect2!=i_anchor){
                    if(mu1<mu2){
                        i_chosen=i_bisect1;
                    }
                    else{
                        i_chosen=i_bisect2;
                    }
                }
                else{
                    // printf("could not find good i_bisect\n");
                    exit(1);
                }

                if(i_chosen>=0){
                    for(j=0;j<_chifn->get_dim();j++){
                        bisect_dir.set(j,0.5*(_chifn->get_pt(i_anchor,j)+_chifn->get_pt(i_chosen,j)));
                    }
                    evaluate(bisect_dir,&mu1,&j);
                    seed.add_row(bisect_dir);
                    kept_dex.add(i_chosen);
                }
            }

            if(seed.get_rows()==_chifn->get_dim()){
                avg_dd=0.0;
                for(i=0;i<_chifn->get_dim();i++){
                    dd=0.0;
                    for(j=0;j<_chifn->get_dim();j++){
                        dd+=power(_chifn->get_pt(i_anchor,j)-seed.get_data(i,j),2);
                    }
                    avg_dd+=sqrt(dd);
                }
                avg_dd/=value_t(_chifn->get_dim());
                for(i=0;i<_chifn->get_dim();i++){
                    avg_pt.set(i,_chifn->get_pt(specified,i)+avg_dd*base_dir.get_data(i));
                }
                evaluate(avg_pt,&mu1,&i);
                // printf("    fn_avg %e\n",dchifn(avg_pt));
                seed.add_row(avg_pt);
            }
        }
    }
    else{
        // printf("calling _explorers.get_seed(); did not expect that\n");
        exit(1);
    }

    // printf("fn anchor %e; %d %d %d\n",
    // _chifn->get_fn(i_anchor),specified,i_origin,i_anchor);

    _chifn->set_search_type(old_type);

    value_t mu_min;
    value_t start_min,start_max;
    index_t i_start_min;
    array_1d<value_t> start_vals,start_vals_sorted;
    array_1d<index_t> start_vals_dex;
    start_vals.set_name("simplex_boundary_start_vals");
    start_vals_sorted.set_name("simplex_boundary_start_vals_sorted");
    start_vals_dex.set_name("simplex_boundary_start_vals_dex");
    for(i=0;i<seed.get_rows();i++){
        mu=dchifn(seed(i));
        start_vals.add(mu);
        start_vals_dex.add(i);
        if(i==0 || mu<start_min){
            start_min=mu;
            _chifn->evaluate(seed(i),&mu,&i_start_min);
        }
        if(i==0 || mu>start_max){
            start_max=mu;
        }
    }
    sort(start_vals,start_vals_sorted,start_vals_dex);
    mu=start_vals_sorted.get_data(start_vals_dex.get_dim()/2);
    // printf("    starting from %e; %e; %e\n",start_min,mu,start_max);
    // printf("    starting from (delta) %e; %e; %e\n",start_min-target(),mu-target(),start_max-target());

    array_1d<value_t> minpt;
    minpt.set_name("dalex_simplex_search_minpt");

    // loop over points, calling dchifn so that they
    // get stored in dchifn's cache
    for(i=pt_start;i<_chifn->get_pts();i++){
        mu=dchifn(_chifn->get_pt(i));
    }

    ffmin.find_minimum(seed,minpt);

    array_1d<index_t> path_row;
    path_row.set_name("path_row");

    i_next[0]=-1;
    value_t fn,cost_min;
    index_t valid_cache;
    cost_min=2.0*exception_value;
    for(i=pt_start;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            valid_cache=dchifn.get_cached_values(i,&fn);
            if(valid_cache==0){
                // printf("bad cache %d %d\n",i,pt_start);
                exit(1);
            }
            if(valid_cache==1 && fn<cost_min){
                i_next[0]=i;
                cost_min=fn;
            }
        }
    }
    if(i_next[0]<0){
        evaluate(minpt,&mu,i_next);
    }

    index_t pre_fill=_chifn->get_pts();
    array_1d<value_t> perp_dir;
    array_1d<value_t> trial_dir;
    perp_dir.set_name("dalex_simp_bou_perp_dir");
    trial_dir.set_name("dalex_simp_bou_trial_dir");
    value_t base_norm;
    index_t iteration;
    value_t step;
    value_t rr;
    index_t n_good_fill=0;
    if(i_next[0]!=specified){
        for(i=0;i<_chifn->get_dim();i++){
            base_dir.set(i,_chifn->get_pt(i_next[0],i)-_chifn->get_pt(specified,i));
        }
        base_norm=base_dir.normalize();
        for(iteration=0;iteration<2*_chifn->get_dim();iteration++){
            for(i=0;i<_chifn->get_dim();i++){
                perp_dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
            }
            component=0.0;
            for(i=0;i<_chifn->get_dim();i++){
                component+=perp_dir.get_data(i)*base_dir.get_data(i);
            }
            for(i=0;i<_chifn->get_dim();i++){
                perp_dir.subtract_val(i,component*base_dir.get_data(i));
            }
            perp_dir.normalize();
            rr=_chifn->random_double();
            for(i=0;i<_chifn->get_dim();i++){
                trial_dir.set(i,base_dir.get_data(i)+rr*perp_dir.get_data(i));
            }
            trial_dir.normalize();
            for(step=0.1*base_norm;step<1.05*base_norm;step+=0.1*base_norm){
                for(i=0;i<_chifn->get_dim();i++){
                    trial.set(i,_chifn->get_pt(specified,i)
                               +step*trial_dir.get_data(i));
                }
                evaluate(trial,&mu,&i);
                if(mu<target()){
                    n_good_fill++;
                }
            }
        }
    }

    // printf("    filling called %d good %d\n",_chifn->get_pts()-pre_fill,n_good_fill);

    _update_good_points(pt_start);

    for(i=pre_fill;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            fn=dchifn(_chifn->get_pt(i));
            if(fn<cost_min){
                i_next[0]=i;
                cost_min=fn;
            }
        }
    }
    if(i_next[0]<0){
        // printf("WARNING could not set i_next\n");
        exit(1);
    }

    if(_chifn->get_dim()>9){
        // printf("    actually found %e -- %.3e %.3e; %.3e %.3e\n",
        // _chifn->get_fn(i_next[0]),_chifn->get_pt(i_next[0],6), _chifn->get_pt(i_next[0],9),
        // _chifn->get_pt(i_next[0],0), _chifn->get_pt(i_next[0], 1));
    }

    // printf("    adjusted %e from %e\n",
    // dchifn(_chifn->get_pt(i_next[0])),_chifn->get_fn(i_next[0]));

    // printf("    min is %e target %e\n",chimin(),target());
    if(_chifn->get_dim()>9){
       // printf("    minpt at %e %e\n",
    //    _chifn->get_pt(mindex(),6),
    //    _chifn->get_pt(mindex(),9));
    }


    index_t path_start,path_end,path_mid;
    path_start=-1;
    path_end=-1;
    path_mid=-1;
    array_1d<value_t> path_dist,path_dist_sorted;
    array_1d<index_t> path_dist_dex;
    path_dist.set_name("dalex_simp_bou_path_dist");
    path_dist_sorted.set_name("dalex_simp_bou_path_dist_sorted");
    path_dist_dex.set_name("dalex_simp_bou_path_dist_dex");

    for(i=pt_start;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            mu=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                mu+=(_chifn->get_pt(i,j)-_chifn->get_pt(specified,j))*base_dir.get_data(j);
            }
            path_dist.add(mu);
            path_dist_dex.add(i);
        }
    }
    if(path_dist.get_dim()>1){
        sort(path_dist, path_dist_sorted, path_dist_dex);
    }

    if(path_dist.get_dim()>0){
        path_start=path_dist_dex.get_data(0);
    }
    if(path_dist.get_dim()>1){
        path_end=path_dist_dex.get_data(path_dist.get_dim()-1);
    }
    if(path_dist.get_dim()>2){
        path_mid=path_dist_dex.get_data(path_dist.get_dim()/2);
    }

    value_t dd_start,dd_mid,dd_end;

    if(path_start>=0 || path_end>=0 || path_mid>0){
        // printf("going to try to add to path %d %d\n",
        // path_start,_chifn->get_pts());
        // printf("start %d end %d\n",path_start,path_end);

        for(i=path_start;i<_chifn->get_pts();i++){
            dd_start=-1.0;
            dd_end=-1.0;
            dd_mid=-1.0;
            if(_tendril_path.get_rows()==0){
                _tendril_path.set_cols(2);
            }

            if(_chifn->get_fn(i)<_chifn->target() &&
               _chifn->get_search_type_log(i)==_type_tendril){

                path_row.set(0,i);

                if(path_start>=0){
                    dd_start=distance(i,path_start);
                }
                if(path_end>=0){
                    dd_end=distance(i,path_end);
                }
                if(path_mid>=0){
                    dd_mid=distance(i,path_mid);
                }

                if(dd_start<0.0 && dd_end<0.0 && dd_mid<0.0){
                    // printf("WARNING constructing path but all dd are <0\n");
                    exit(1);
                }

                if(dd_start<0.0){
                    if(dd_end<0.0){
                        path_row.set(1,path_mid);
                    }
                    else if(dd_mid<0.0){
                        path_row.set(1,path_end);
                    }
                    else if(dd_mid<dd_end){
                        path_row.set(1,path_mid);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }
                else if(dd_mid<0.0){
                    if(dd_start<0.0){
                        path_row.set(1,path_end);
                    }
                    else if(dd_end<0.0){
                        path_row.set(1,path_start);
                    }
                    else if(dd_start<dd_end){
                        path_row.set(1,path_start);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }
                else if(dd_end<0.0){
                    if(dd_start<0.0){
                        path_row.set(1,path_mid);
                    }
                    else if(dd_mid<0.0){
                        path_row.set(1,path_start);
                    }
                    else if(dd_start<dd_mid){
                        path_row.set(1,path_start);
                    }
                    else{
                        path_row.set(1,path_mid);
                    }
                }
                else{
                    if(dd_start<=dd_end && dd_start<=dd_mid){
                        path_row.set(1,path_start);
                    }
                    else if(dd_mid<=dd_end && dd_mid<=dd_start){
                        path_row.set(1,path_mid);
                    }
                    else{
                        path_row.set(1,path_end);
                    }
                }

                if(path_row.get_data(1)<0){
                    // printf("WARNING somehow path_row(1)<0\n");
                    exit(1);
                }
                _tendril_path.add_row(path_row);
            }
        }
    }

    index_t is_a_strike=0;

    for(i=0;i<exclusion_zones.ct() && is_a_strike==0;i++){
        if(exclusion_zones(i)->contains(_chifn->get_pt(i_next[0]))==1){
            is_a_strike=1;
        }
    }

    if(is_a_strike==1){
        return 1;
    }
    else{
        return 0;
    }

}

index_t dalex::_exploration_simplex(index_t i1, index_t i0, array_1d<index_t> &associates){

    array_2d<value_t> associate_pts;
    associate_pts.set_name("dalex_exp_simp_associate_pts");
    index_t i;
    for(i=0;i<associates.get_dim();i++){
        associate_pts.add_row(_chifn->get_pt(associates.get_data(i)));
    }
    cost_fn dchifn(_chifn, associates);

    array_2d<value_t> seed;
    seed.set_name("dalex_exp_sim_seed");
    array_1d<value_t> min,max;
    min.set_name("dalex_exp_simp_min");
    max.set_name("dalex_exp_simp_max");

    array_1d<value_t> dir;
    dir.set_name("dalex_exp_sim_dir");
    array_2d<value_t> dir_log;
    dir_log.set_name("dalex_exp_sim_dir_log");
    for(i=0;i<_chifn->get_dim();i++){
        dir.set(i,_chifn->get_pt(i1,i)-_chifn->get_pt(i0,i));
    }
    value_t dir_norm;
    dir_norm=dir.normalize();
    array_1d<value_t> newpt,newdir;
    newpt.set_name("dalex_exp_simp_newpt");
    newdir.set_name("dalex_exp_simp_newdir");
    for(i=0;i<_chifn->get_dim();i++){
        newpt.set(i,_chifn->get_pt(i1,i)+0.1*dir_norm*dir.get_data(i));
    }
    seed.add_row(newpt);
    dir_log.add_row(dir);
    value_t component;
    index_t j;
    while(seed.get_rows()<_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(),0.0,1.0));
        }
        for(i=0;i<dir_log.get_rows();i++){
            component=0.0;
            for(j=0;j<_chifn->get_dim();j++){
                component+=dir.get_data(j)*dir_log.get_data(i,j);
            }
            for(j=0;j<_chifn->get_dim();j++){
                dir.subtract_val(j,component*dir_log.get_data(i,j));
            }
        }
        component=dir.normalize();
        if(component>1.0e-20){
            for(i=0;i<_chifn->get_dim();i++){
                newdir.set(i,dir.get_data(i)+dir_log.get_data(0,i));
            }
            newdir.normalize();
            for(i=0;i<_chifn->get_dim();i++){
                newpt.set(i,_chifn->get_pt(i1,i)+0.1*dir_norm*newdir.get_data(i));
            }
            seed.add_row(newpt);
            dir_log.add_row(dir);
        }
    }
    for(i=0;i<_chifn->get_dim();i++){
        newpt.set(i,0.0);
    }
    for(i=0;i<seed.get_rows();i++){
        for(j=0;j<_chifn->get_dim();j++){
            newpt.add_val(j,seed.get_data(i,j));
        }
    }
    for(i=0;i<_chifn->get_dim();i++){
        newpt.divide_val(i,value_t(seed.get_rows()));
    }
    seed.add_row(newpt);

    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_chisquared(&dchifn);
    ffmin.set_dice(_chifn->get_dice());
    ffmin.set_minmax(min,max);
    ffmin.use_gradient();
    ffmin.find_minimum(seed,newpt);
    value_t mu;
    index_t i_found;
    evaluate(newpt,&mu,&i_found);
    // printf("exp simplex found %e\n",mu);
    return i_found;

}


void dalex::find_tendril_candidates(value_t factor_in){

    // printf("finding tendril candidates\n");

    array_1d<index_t> associates;
    associates.set_name("find_tendrils_associates");
    array_1d<index_t> associates_raw;
    associates_raw.set_name("find_tendril_associates_raw");
    array_2d<value_t> ellipse_pts;
    ellipse_pts.set_name("find_tendrils_ellipse_pts");
    array_1d<value_t> center;
    center.set_name("find_tendrils_center");
    index_t i,j;
    for(i=0;i<_chifn->get_pts();i++){
        if(_chifn->get_fn(i)<target()){
            ellipse_pts.add_row(_chifn->get_pt(i));
            associates.add(i);
            associates_raw.add(i);
        }
    }
    if(ellipse_pts.get_rows()==0){
        // printf("CANNOT find tendrils; not good points\n");
        exit(1);
    }
    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,_chifn->get_pt(mindex(),i));
    }

    index_t n_thin;

    if(associates.get_dim()>20000){
        n_thin=associates_raw.get_dim()/20000;
        if(n_thin==1){
            n_thin=2;
        }
        associates.reset_preserving_room();
        for(i=0;i<associates_raw.get_dim();i++){
            if(i%n_thin==0){
                associates.add(associates_raw.get_data(i));
            }
        }
    }

    ellipse good_ellipse;
    good_ellipse.build(center,ellipse_pts);

    cost_fn dchifn(_chifn,associates,0);
    value_t envelope=0.25*(target()-chimin());
    if(envelope<2.0){
        envelope=2.0;
    }
    dchifn.set_envelope(envelope);

    array_1d<value_t> min,max;
    min.set_name("find_tendrils_min");
    max.set_name("find_tendrils_max");
    for(i=0;i<_chifn->get_dim();i++){
        min.set(i,0.0);
        max.set(i,_chifn->get_characteristic_length(i));
    }

    simplex_minimizer ffmin;
    ffmin.set_minmax(min,max);
    ffmin.set_chisquared(&dchifn);

    array_1d<value_t> trial,minpt;
    trial.set_name("find_tendrils_trial");
    minpt.set_name("find_tendrils_minpt");
    value_t mu;
    index_t i_found;

    array_1d<index_t> particles,origins;
    particles.set_name("find_tendrils_particles");
    origins.set_name("find_tendrils_origins");
    array_1d<value_t> fn_val;
    fn_val.set_name("find_tendrils_fn_val");
    array_1d<index_t> fn_val_dex;
    fn_val_dex.set_name("find_tendrils_fn_val_dex");

    array_2d<value_t> seed;
    seed.set_name("find_tendrils_seed");

    index_t idim,jdim;
    value_t sgn,cost;
    index_t pt_start;
    array_1d<index_t> at_least_one;
    at_least_one.set_name("find_tendril_at_least_one");
    value_t factor;
    index_t failures;

    // printf("n associates %d\n",associates.get_dim());
    for(idim=0;idim<_chifn->get_dim();idim++){
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            pt_start=_chifn->get_pts();
            seed.reset_preserving_room();
            i_found=-1;
            factor=factor_in;
            while(i_found<0){
                // // printf("    running with factor %e\n",factor);
                for(i=0;i<_chifn->get_dim();i++){
                    trial.set(i,center.get_data(i)+sgn*factor*good_ellipse.radii(idim)*good_ellipse.bases(idim,i));
                }
                evaluate(trial,&mu,&i_found);
                // // printf("    first point is %d\n",i_found);
                if(factor>2.0){
                    factor-=1.0;
                }
                else{
                    factor*=0.5;
                }
            }
            origins.add(i_found);
            seed.add_row(trial);
            for(jdim=0;jdim<_chifn->get_dim();jdim++){
                for(i=0;i<_chifn->get_dim();i++){
                   trial.set(i,seed.get_data(0,i)+0.1*good_ellipse.radii(jdim)*good_ellipse.bases(jdim,i));
                }
                seed.add_row(trial);
            }
            failures=0;
            ffmin.find_minimum(seed,minpt);
            evaluate(minpt,&mu,&i_found);
            failures++;
            while(_chifn->get_fn(i_found)>target() && failures<3){
                // // printf("    trying again because %e > %e\n",
                // _chifn->get_fn(i_found),target());
                seed.reset_preserving_room();
                seed.add_row(minpt);
                for(jdim=0;jdim<_chifn->get_dim();jdim++){
                    for(i=0;i<_chifn->get_dim();i++){
                        trial.set(i,seed.get_data(0,i)+0.05*good_ellipse.radii(jdim)*good_ellipse.bases(jdim,i));
                    }
                    seed.add_row(trial);
                }
                ffmin.set_dice(_chifn->get_dice());
                ffmin.use_gradient();
                ffmin.find_minimum(seed,minpt);
                evaluate(minpt,&mu,&i_found);
                ffmin.do_not_use_gradient();
                failures++;
            }
            if(_chifn->get_fn(i_found)<target()){
                at_least_one.add(i_found);
            }
            particles.add(i_found);
            cost=dchifn(minpt);
            fn_val.add(cost);
            // if(_chifn->get_dim()>9){
                // // printf("    got %d %e -- %e %e\n\n",particles.get_dim(),mu,minpt.get_data(6),minpt.get_data(9));
            // }

            j=0;
            for(i=pt_start;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target()){
                    associates.add(i);
                    associates_raw.add(i);
                    j++;
                }
            }
            if(j>0){
                if(associates.get_dim()>20000){
                    n_thin=associates_raw.get_dim()/20000;
                    if(n_thin==1){
                        n_thin=2;
                    }
                    associates.reset_preserving_room();
                    for(i=0;i<at_least_one.get_dim();i++){
                        associates.add(at_least_one.get_data(i));
                    }
                    for(i=0;i<associates_raw.get_dim();i++){
                        if(i%n_thin==0 && associates.contains(associates_raw.get_data(i))==0){
                            associates.add(associates_raw.get_data(i));
                        }
                    }
                }
                dchifn.build(_chifn,associates,0);
            }

            if(_limit>0 && _chifn->get_pts()>_limit){
                _chifn->write_pts();
                return;
            }
        }
    }
    for(i=0;i<fn_val.get_dim();i++){
        fn_val_dex.add(i);
    }
    array_1d<value_t> fn_val_sorted;
    fn_val_sorted.set_name("find_tendrils_fn_val_sorted");
    sort(fn_val,fn_val_sorted,fn_val_dex);

    for(i=0;i<_chifn->get_dim()/2;i++){
        _particle_candidates.set(i,particles.get_data(fn_val_dex.get_data(i)));
        _origin_candidates.set(i,-1);
    }
}


void dalex::get_new_tendril(index_t *particle, index_t *origin){
    particle[0]=-1;
    origin[0]=-1;
    index_t i,j,ip,io;
    index_t is_outside;
    index_t old_type;

    cost_fn dchifn;
    array_1d<index_t> associates;
    array_1d<value_t> cost_val,cost_val_sorted;
    array_1d<index_t> cost_val_dex;
    associates.set_name("new_tendril_associates");
    cost_val.set_name("new_tendril_cost_val");
    cost_val_sorted.set_name("new_tendril_cost_val_sorted");
    cost_val_dex.set_name("new_tendril_cost_val_dex");

    while(particle[0]<0){
        if(_limit>0 && _chifn->get_pts()>_limit){
            _chifn->write_pts();
            return;
        }
        associates.reset_preserving_room();
        cost_val.reset_preserving_room();
        cost_val_sorted.reset_preserving_room();
        cost_val_dex.reset_preserving_room();
        for(i=0;i<_chifn->get_pts();i++){
            if(_chifn->get_fn(i)<target()){
                if(_chifn->get_search_type_log(i)!=_type_init_tendril){
                    associates.add(i);
                }
            }
        }
        dchifn.build(_chifn,associates,1);
        for(i=0;i<_particle_candidates.get_dim();i++){
            if(_particle_candidates.get_data(i)<0){
                cost_val.add(2.0*exception_value);
            }
            else{
                cost_val.add(dchifn(_chifn->get_pt(_particle_candidates.get_data(i))));
            }
            cost_val_dex.add(i);
        }
        sort(cost_val, cost_val_sorted, cost_val_dex);

        for(i=0;i<_particle_candidates.get_dim();i++){
            ip=_particle_candidates.get_data(cost_val_dex.get_data(i));
            io=_origin_candidates.get_data(cost_val_dex.get_data(i));
            if(ip>=0){
                is_outside=1;
                for(j=0;j<_exclusion_zones.ct();j++){
                    if(_exclusion_zones(j)->contains(_chifn->get_pt(ip))==1){
                        is_outside=0;
                        break;
                    }
                }
                if(is_outside==1){
                    particle[0]=ip;
                    origin[0]=io;
                    _particle_candidates.set(cost_val_dex.get_data(i),-1);
                    _origin_candidates.set(cost_val_dex.get_data(i),-1);
                    // printf("returning tendril %d %d\n",particle[0],origin[0]);
                    return;
                }
                else{
                    _particle_candidates.set(cost_val_dex.get_data(i),-1);
                    _origin_candidates.set(cost_val_dex.get_data(i),-1);
                }
            }
        }
        old_type=_chifn->get_search_type();
        _chifn->set_search_type(_type_init_tendril);
        find_tendril_candidates(3.0);
        _chifn->set_search_type(old_type);
    }
}


void dalex::init_fill(){
    index_t path_len=_tendril_path.get_rows();
    index_t old_type=_chifn->get_search_type();
    _chifn->set_search_type(_type_init);
    find_tendril_candidates(1.0);
    index_t i;
    for(i=0;i<_particle_candidates.get_dim();i++){
        _particle_candidates.set(i,-1);
        _origin_candidates.set(i,-1);
    }

    _chifn->set_search_type(old_type);
}

void dalex::tendril_search(){
    index_t pt_start=_chifn->get_pts();
    index_t pt_prime=_chifn->get_pts();
    array_1d<value_t> dir,midpt;
    index_t i_found,i_next;
    value_t mu;
    dir.set_name("tendril_dir");
    midpt.set_name("tendril_mid_pt");
    index_t i,j;

    array_2d<value_t> ellipse_pts;
    ellipse_pts.set_name("tendril_ellipse_pts");
    ellipse local_ellipse;
    index_t is_a_strike;

    array_1d<index_t> new_particles,new_origins;
    new_particles.set_name("tendril new_particles");
    new_origins.set_name("new_origins");

    index_t particle, origin;

    if(_tendril_init==0){
        init_fill();
        _tendril_init=1;
        if(_limit>0 && _chifn->get_pts()>_limit){
            _chifn->write_pts();
            return;
        }
    }

    index_t strikes=0;
    _has_struck=0;

    get_new_tendril(&particle,&origin);

    while(strikes<3){

        ellipse_pts.reset_preserving_room();
        is_a_strike=0;
        pt_start=_chifn->get_pts();

        is_a_strike=simplex_boundary_search(particle,origin,_exclusion_zones,&i_next);

        if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
            if(local_ellipse.contains(_chifn->get_pt(i_next),1)==1){
                is_a_strike=1;
            }
        }

        for(j=pt_start;j<_chifn->get_pts();j++){
            if(_chifn->get_fn(j)<target()){
                ellipse_pts.add_row(_chifn->get_pt(j));
            }
        }

        if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
            local_ellipse.build(ellipse_pts);
        }

        if(is_a_strike==0){
            origin=particle;
            particle=i_next;
            if(_chifn->get_fn(i_next)<target()+1.0e-6){
                strikes=0;
            }
            else{
                is_a_strike=1;
                strikes++;
            }
        }
        else{
             strikes++;
             _has_struck=1;
        }

        // printf("pts %d lim %d\n",_chifn->get_pts(),_limit);
        if(_chifn->get_dim()>9)// printf("got to %e %e\n",_chifn->get_pt(i_next,6),_chifn->get_pt(i_next,9));
        // printf("is a strike: %d; strikes %d\n",is_a_strike,strikes);
        // printf("\n");

        if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
            local_ellipse.build(ellipse_pts);
            _exclusion_zones.add(local_ellipse);
        }

    }


    // printf("\n    pts %d limit %d\n",_chifn->get_pts(),_limit);
    // printf("    strikeouts %d\n",_strikeouts);
    _chifn->write_pts();
    _update_good_points();
}

void dalex::initialize_min_exploration(){
    _min_explorers.initialize();
    array_1d<value_t> dir,midpt;
    dir.set_name("dalex_init_min_exp_dir");
    midpt.set_name("dalex_init_min_exp_midpt");
    value_t mu;
    index_t i_found;
    index_t i,j;
    while(_min_explorers.get_n_particles()<2*_chifn->get_dim()){
        for(i=0;i<_chifn->get_dim();i++){
            dir.set(i,normal_deviate(_chifn->get_dice(), 0.0, 1.0));
        }
        dir.normalize();
        i_found=bisection(mindex(),dir,target(),0.001);
        if(i_found!=mindex()){
            for(i=0;i<_chifn->get_dim();i++){
                midpt.set(i,0.5*(_chifn->get_pt(mindex(),i)+_chifn->get_pt(i_found,i)));
            }
            evaluate(midpt,&mu,&i_found);
            if(i_found!=mindex()){
                _min_explorers.add_particle(_chifn->get_pt(i_found));
            }
        }
    }

}

void dalex::min_explore(index_t n_particles, index_t n_steps){
    // printf("\nmin exploring\n");
    index_t pt_0=_chifn->get_pts();

    if(_min_explorers.get_n_particles()==0){
        initialize_min_exploration();
    }

    array_1d<index_t> associates;
    associates.set_name("dalex_explore_associates");
    index_t skip;
    if(_good_points.get_dim()<5000){
        skip=-1;
    }
    else{
        skip=5;
    }
    index_t i;
    for(i=0;i<_good_points.get_dim();i++){
        if(skip<0 || i%skip==0){
            associates.add(_good_points.get_data(i));
        }
    }

    _min_explorers.set_associates(associates);
    _min_explorers.sample(n_steps);

    _update_good_points(pt_0);
}


void dalex::iterate_on_minimum(){
    // printf("iterating with %e\n",chimin());
    value_t min_0=chimin();
    value_t min_00=min_0;

    value_t min_1=-2.0*exception_value;
    index_t i,j;

    if(_reset_chimin>exception_value){
        _reset_chimin=chimin();
    }

    if(_good_points.get_dim()==0){
        find_bases();
    }

    while(min_1<min_0){
        min_0=chimin();
        min_explore(2*_chifn->get_dim(), 4*_chifn->get_dim());
        simplex_search(mindex());
        min_1=chimin();
    }

    array_2d<index_t> tendril_cache;
    tendril_cache.set_name("tendril_cache");
    array_1d<index_t> good_cache;
    good_cache.set_name("good_cache");
    array_1d<index_t> origins;
    origins.set_name("iterate_min_origins");
    array_2d<value_t> ellipse_pts;
    ellipse_pts.set_name("iterate_ellipse_pts");
    ellipse local_ellipse;

    if(chimin()<_reset_chimin-_reset_threshold){

        for(i=0;i<_tendril_path.get_rows();i++){
            if(_chifn->get_fn(_tendril_path.get_data(i,0))<target()+1.0e-6){
                tendril_cache.add_row(_tendril_path(i));
            }
        }

        for(i=0;i<_good_points.get_dim();i++){
            if(_chifn->get_fn(_good_points.get_data(i))<target()){
                good_cache.add(_good_points.get_data(i));
            }
        }

        _reset_chimin=chimin();
        _good_points.reset_preserving_room();
        _tendril_path.reset_preserving_room();
        _exclusion_zones.reset();

        for(i=0;i<good_cache.get_dim();i++){
            _good_points.add(good_cache.get_data(i));
        }

        for(i=0;i<tendril_cache.get_rows();i++){
            _tendril_path.add_row(tendril_cache(i));
            if(origins.contains(tendril_cache.get_data(i,1))==0){
                origins.add(tendril_cache.get_data(i,1));
            }
        }

        for(i=0;i<origins.get_dim();i++){
            ellipse_pts.reset_preserving_room();
            for(j=0;j<_tendril_path.get_rows();j++){
                if(_tendril_path.get_data(j,1)==origins.get_data(i)){
                    ellipse_pts.add_row(_chifn->get_pt(_tendril_path.get_data(j,0)));
                }
            }
            if(ellipse_pts.get_rows()>2*_chifn->get_dim()){
                local_ellipse.build(ellipse_pts);
                _exclusion_zones.add(local_ellipse);
            }
        }

        find_bases();
    }
    else if(chimin()<min_00-0.01){
        find_bases();
    }

    // printf("done iterating %e %d\n",chimin(),_chifn->get_called());

}

void dalex::refine_minimum(){

    // printf("    refining minimum %e\n",chimin());

    index_t n_particles=2*_chifn->get_dim();
    index_t n_steps=5*_chifn->get_dim();
    index_t accepted=0;
    index_t rejected=0;

    array_1d<value_t> trial,coord;
    trial.set_name("refine_minimum_trial");
    coord.set_name("refine_minimum_coord");
    value_t roll,ratio,mu;
    index_t i_found;
    index_t i,j;

    _minimizers.reset_preserving_room();
    while(_minimizers.get_dim()!=n_particles){
        for(i=0;i<_chifn->get_dim();i++){
           coord.set(i,2.0*(_chifn->random_double()-0.5)*
                       _basis_associate_norm.get_data(i));
        }

        for(i=0;i<_chifn->get_dim();i++){
            trial.set(i,_chifn->get_pt(mindex(),i));
        }

        for(i=0;i<_chifn->get_dim();i++){
            for(j=0;j<_chifn->get_dim();j++){
                trial.add_val(j,coord.get_data(i)*_basis_vectors.get_data(i,j));
            }
        }

        evaluate(trial, &mu, &i_found);
        if(i_found>=0){
            _minimizers.add(i_found);
        }
    }

    index_t i_step,ip,i0;
    index_t i_dim,accept_it;
    value_t rr;
    for(i_step=0;i_step<n_steps;i_step++){
        for(ip=0;ip<_minimizers.get_dim();ip++){
            i0=_minimizers.get_data(ip);
            i_dim=_chifn->random_int()%_chifn->get_dim();
            rr=normal_deviate(_chifn->get_dice(),0.0,1.0)*0.5*_basis_associate_norm.get_data(i_dim);
            for(i=0;i<_chifn->get_dim();i++){
                trial.set(i,_chifn->get_pt(i0,i));
            }
            for(i=0;i<_chifn->get_dim();i++){
                trial.add_val(i,rr*_basis_vectors.get_data(i_dim,i));
            }
            evaluate(trial,&mu,&i_found);
            accept_it=0;
            if(i_found>=0){
                if(mu<_chifn->get_fn(i0)){
                    accept_it=1;
                }
                else{
                    roll=_chifn->random_double();
                    ratio=exp(-0.5*(mu-_chifn->get_fn(i0)));
                    if(roll<ratio){
                        accept_it=1;
                    }
                }
            }

            if(accept_it==1){
                accepted++;
                _minimizers.set(ip,i_found);
            }
            else{
                rejected++;
            }
        }
    }

    value_t f_min=exception_value;
    for(i=0;i<_minimizers.get_dim();i++){
        mu=_chifn->get_fn(_minimizers.get_data(i));
        if(mu<f_min){
            f_min=mu;
        }
    }

    // printf("    refined to %e %d %d %e\n",chimin(),accepted,rejected,f_min);

}

void dalex::compass_search(ellipse &ee){
    array_1d<value_t> center;
    center.set_name("compass_search_center");
    array_1d<value_t> dir;
    dir.set_name("compass_search_dir");
    index_t i_found;
    value_t mu;
    index_t i;
    for(i=0;i<_chifn->get_dim();i++){
        center.set(i,ee.center(i));
    }
    evaluate(center,&mu,&i_found);
    if(mu>target()){
        // printf("   cannot do compass: %e\n",mu);
        return;
    }

    index_t i_dim;
    value_t sgn;
    for(i_dim=0;i_dim<_chifn->get_dim();i_dim++){
        for(sgn=-1.0;sgn<1.1;sgn+=2.0){
            for(i=0;i<_chifn->get_dim();i++){
                dir.set(i,sgn*ee.bases(i_dim,i));
            }
            bisection(i_found, dir, target(), 0.01);
        }
    }
}
