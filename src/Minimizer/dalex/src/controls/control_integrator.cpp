#include "controls/control_integrator.h"

control_integrator::control_integrator(function_wrapper &chisq,
                     array_1d<value_t> &min, array_1d<value_t> &max,
                     array_1d<value_t> &dx, char *name_root){

    _iter_p = NULL;
    _d_threshold=60.0;
    _min.set_name("control_min");
    _max.set_name("control_max");
    _dx.set_name("control_dx");

    chisquared_distribution distro;
    _max_chi_lim_freq=distro.confidence_limit(100.0,0.95);

    _chisq=&chisq;
    index_t i;
    for(i=0;i<min.get_dim();i++){
        _min.set(i,min.get_data(i));
        _max.set(i,max.get_data(i));
        _dx.set(i,dx.get_data(i));
    }

    for(i=0;i<letters-1 && name_root[i]!=0;i++){
        _name_root[i]=name_root[i];
    }
    _name_root[i]=0;
    _chi_min=2.0*exception_value;

}

value_t control_integrator::get_min(index_t ii){
    return _min.get_data(ii);
}

value_t control_integrator::get_max(index_t ii){
    return _max.get_data(ii);
}

void control_integrator::run_analysis(value_t cc){
    array_1d<value_t> cc_arr;
    cc_arr.add(cc);
    run_analysis(cc_arr);
}

void control_integrator::run_analysis(){
    run_analysis(0.95);
}

void control_integrator::_initialize_iterate(){
    if(_iter_p!=NULL){
        delete _iter_p;
    }

    _iter_p=new default_iteration_parameters();
    _iter_p->initialize(_min, _max, _dx);
}


void control_integrator::run_analysis(array_1d<value_t> &cc){

    value_t max_confidence_limit;
    index_t ic;
    for(ic=0;ic<cc.get_dim();ic++){
        if(ic==0 || cc.get_data(ic)>max_confidence_limit){
            max_confidence_limit=cc.get_data(ic);
        }
    }


    printf("nameroot %s\n",_name_root);
    printf("max_chi_lim_freq %e\n",_max_chi_lim_freq);

    array_1d<value_t> pt;

    array_1d<value_t> chi_vals,good_min,good_max;
    array_2d<index_t> coordinates;

    pt.set_name("pt");
    good_min.set_name("good_min");
    good_max.set_name("good_max");
    chi_vals.set_name("chi_vals");
    coordinates.set_name("coordinates");
    coordinates.set_dim(1000000,_min.get_dim());
    chi_vals.set_dim(1000000);
    coordinates.reset_preserving_room();
    chi_vals.reset_preserving_room();

    value_t xx,yy,zz;
    index_t ix,iy,iz;

    for(ix=0;ix<_min.get_dim();ix++){
        good_min.set(ix,2.0*exception_value);
        good_max.set(ix,-2.0*exception_value);
    }

    value_t mu;
    value_t start=value_t(time(NULL));

    index_t ct=0;

    value_t foundMin=exception_value;

    index_t keep_going=1;

    _initialize_iterate();

    array_1d<index_t> idx;
    idx.set_name("control_idx");

    value_t expected_time;

    value_t good_threshold=_max_chi_lim_freq+_d_threshold;
    printf("good_threshold %e\n",good_threshold);

    while(keep_going==1){

        keep_going=_iter_p->get_pt(pt,idx);

        mu=_chisq[0](pt);

        if(mu<foundMin)foundMin=mu;

        if(mu<good_threshold){
            for(ix=0;ix<_min.get_dim();ix++){
                if(pt.get_data(ix)<good_min.get_data(ix)){
                    good_min.set(ix,pt.get_data(ix));
                }
                if(pt.get_data(ix)>good_max.get_data(ix)){
                    good_max.set(ix,pt.get_data(ix));
                }
            }
            coordinates.add_row(idx);
            chi_vals.add(mu);
        }

        if(_iter_p->get_current_ct()%1000000==0 && _iter_p->get_current_ct()>0){

            expected_time = ((value_t(time(NULL))-start)*_iter_p->get_total_ct())/value_t(_iter_p->get_current_ct());

            // printf("ct %ld good %d room %d in %.4e %.4e -- %.4e -- expect %e hours\n",
            // _iter_p->get_current_ct(),chi_vals.get_dim(),chi_vals.get_room(),(value_t(time(NULL))-start)/3600.0,
            // (value_t(time(NULL))-start)/value_t(_iter_p->get_current_ct()),foundMin,
            // expected_time/3600.0);
        }


    }

    printf("foundMin %e in %ld\n",foundMin,_iter_p->get_current_ct());

    start=value_t(time(NULL));
    array_1d<value_t> chi_sorted,likelihood_sorted;
    array_1d<index_t> dexes;
    index_t i;
    value_t total=0.0;
    dexes.set_dim(chi_vals.get_dim());
    likelihood_sorted.set_dim(chi_vals.get_dim());
    for(i=0;i<chi_vals.get_dim();i++)dexes.set(i,i);
    sort(chi_vals, chi_sorted, dexes);
    printf("sorting took %e\n",value_t(time(NULL))-start);

    value_t trial_chi_min;
    trial_chi_min=chi_sorted.get_data(0);
    if(trial_chi_min<_chi_min){
        _chi_min=trial_chi_min;
    }
    printf("before simplex chi min is %e\n",_chi_min);
    find_chi_min(good_min,good_max);


    for(i=chi_sorted.get_dim()-1;i>=0;i--){
        chi_sorted.subtract_val(i,_chi_min);
        likelihood_sorted.set(i,exp(-0.5*chi_sorted.get_data(i)));
        total+=exp(-0.5*chi_sorted.get_data(i));
    }

    printf("total %e\n",total);

    char name[letters];

    for(ic=0;ic<cc.get_dim();ic++){
        for(ix=0;ix<_min.get_dim();ix++){
            for(iy=ix+1;iy<_min.get_dim();iy++){
                sprintf(name,"%s_%.2f_%d_%d",_name_root,cc.get_data(ic),ix,iy);
                write_output(ix,iy,chi_sorted,likelihood_sorted,dexes,coordinates,total,
                             cc.get_data(ic),name);

            }
        }
    }

    sprintf(name, "%s_all_good_points.txt", _name_root);
    FILE *output;
    output=fopen(name, "w");
    for(ix=0;ix<coordinates.get_rows();ix++){
        if(chi_vals.get_data(ix)<foundMin+40.0){
            for(iy=0;iy<coordinates.get_cols();iy++){
                fprintf(output,"%e ",_min.get_data(iy)+coordinates.get_data(ix,iy)*_dx.get_data(iy));
            }
            fprintf(output,"%e\n",chi_vals.get_data(ix));
        }
    }
    fclose(output);

}

index_t control_integrator::get_dex(value_t value, value_t min, value_t dx, index_t max){
    index_t base=index_t((value-min)/dx);

    value_t error=min+base*dx-value;
    if(error<-0.5*dx)base--;
    else if(error>0.5*dx)base++;
    else if(fabs(error)>dx){
        // printf("WARNING in get_dex val %e min %e dx %e base %d\n",
        // value,min,dx,base);
        // printf("error %e;   %e\n",error,min+base*dx);
        exit(1);
    }

    if(base<0){
        return 0;
    }
    else if(base>max){
        return max;
    }

    return base;
}


void control_integrator::write_output(index_t xdex, index_t ydex,
                    array_1d<value_t> &chi_vals_sorted,
                    array_1d<value_t> &likelihood_sorted,
                    array_1d<index_t> &dexes, array_2d<index_t> &coordinates,
                    value_t totalLikelihood,
                    value_t confidence_limit,
                    char *outname_root){


    chisquared_distribution distro;

    value_t delta_chi_bayes=distro.confidence_limit(value_t(_min.get_dim()),confidence_limit);
    value_t chi_lim_freq=distro.confidence_limit(100.0,confidence_limit);

    value_t delta_chi_2d=distro.confidence_limit(2.0,confidence_limit);

    printf("minchi %e\n",chi_vals_sorted.get_data(0));
    array_2d<value_t> frequentistFullD,frequentist2D,frequentist2Deff,bayesianFullD,bayesian2D;
    array_2d<value_t> frequentistFullDrelative;
    frequentistFullD.set_name("freq_FullD");
    frequentistFullDrelative.set_name("freq_FullDrel");
    frequentist2D.set_name("freq_2d");
    frequentist2Deff.set_name("freq_2deff");
    bayesianFullD.set_name("bayes_FullD");
    bayesian2D.set_name("bayes_2d");

    frequentistFullD.set_cols(2);
    frequentistFullDrelative.set_cols(2);
    frequentist2D.set_cols(2);
    frequentist2Deff.set_cols(2);
    bayesianFullD.set_cols(2);
    bayesian2D.set_cols(2);

    index_t xct,yct;
    index_t i,j;
    for(i=0;i<coordinates.get_rows();i++){
        if(i==0 || coordinates.get_data(i,xdex)+1>xct){
            xct=coordinates.get_data(i,xdex)+1;
        }

        if(i==0 || coordinates.get_data(i,ydex)+1>yct){
            yct=coordinates.get_data(i,ydex)+1;
        }
    }

    printf("xct %d yct %d\n",xct,yct);

    array_2d<value_t> minChi2Grid,marginalized_likelihood;
    minChi2Grid.set_name("prepare_output_minChi2Grid");
    marginalized_likelihood.set_name("prepare_output_marginalized_likelihood");
    minChi2Grid.set_dim(xct,yct);
    marginalized_likelihood.set_dim(xct,yct);
    for(i=0;i<xct;i++){
        for(j=0;j<yct;j++){
            minChi2Grid.set(i,j,2.0*exception_value);
            marginalized_likelihood.set(i,j,0.0);
        }
    }

    printf("made marginalized_likelihood and minChi2Grid\n");

    index_t row,true_dex,ix,iy;
    for(row=chi_vals_sorted.get_dim()-1;row>=0;row--){
        true_dex=dexes.get_data(row);
        ix=coordinates.get_data(true_dex,xdex);
        iy=coordinates.get_data(true_dex,ydex);


        if(chi_vals_sorted.get_data(row)<minChi2Grid.get_data(ix,iy)){
            minChi2Grid.set(ix,iy,chi_vals_sorted.get_data(row));
        }
        marginalized_likelihood.add_val(ix,iy,likelihood_sorted.get_data(row));
    }

    value_t max_likelihood=-2.0*exception_value;
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
            if(marginalized_likelihood.get_data(ix,iy)>max_likelihood){
                max_likelihood=marginalized_likelihood.get_data(ix,iy);
            }
        }
    }

    printf("max_likelihood %e\n",max_likelihood);

    value_t min_chi2_eff=-2.0*log(max_likelihood);

    //////assemble the raw frequentist grids
    array_1d<value_t> trial,found_min,found_max;
    trial.set_name("trial");
    found_min.set_name("found_min");
    found_max.set_name("found_max");

    for(ix=0;ix<2;ix++){
        found_min.set(ix,2.0*exception_value);
        found_max.set(ix,-2.0*exception_value);
    }

    for(ix=0;ix<xct;ix++){
        trial.set(0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,_min.get_data(ydex)+iy*_dx.get_data(ydex));

            if(minChi2Grid.get_data(ix,iy)+_chi_min<=chi_lim_freq){
                frequentistFullD.add_row(trial);

                for(i=0;i<2;i++){
                    if(trial.get_data(i)<found_min.get_data(i))found_min.set(i,trial.get_data(i));
                    if(trial.get_data(i)>found_max.get_data(i))found_max.set(i,trial.get_data(i));
                }

            }

            if(minChi2Grid.get_data(ix,iy)<=delta_chi_bayes){
                frequentistFullDrelative.add_row(trial);
            }

            if(minChi2Grid.get_data(ix,iy)<=delta_chi_2d){
                frequentist2D.add_row(trial);
            }
        }
    }

    ////assemble the raw bayesian grids
    array_1d<value_t> marginalized_likelihood_line;
    array_2d<value_t> marginalized_likelihood_dexes;
    marginalized_likelihood_line.set_name("marg_like_line");
    marginalized_likelihood_dexes.set_name("marg_like_dex");
    marginalized_likelihood_dexes.set_cols(2);

    row=0;
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
           marginalized_likelihood_line.set(row,marginalized_likelihood.get_data(ix,iy));
           marginalized_likelihood_dexes.set(row,0,ix);
           marginalized_likelihood_dexes.set(row,1,iy);
           row++;
        }
    }

    array_1d<value_t> sorted;
    array_1d<index_t> sorted_dexes;
    sorted.set_name("output_sorted");
    sorted_dexes.set_name("output_sorted_dexes");
    for(i=0;i<marginalized_likelihood_line.get_dim();i++){
        sorted_dexes.set(i,i);
    }

    sort(marginalized_likelihood_line,sorted,sorted_dexes);
    value_t sum=0.0;
    row=0;
    for(i=sorted.get_dim()-1;i>=0 && sum<confidence_limit*totalLikelihood;i--){
        sum+=sorted.get_data(i);

        ix=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),0);
        iy=marginalized_likelihood_dexes.get_data(sorted_dexes.get_data(i),1);

        bayesian2D.set(row,0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        bayesian2D.set(row,1,_min.get_data(ydex)+iy*_dx.get_data(ydex));
        row++;
    }

    printf("Bayesian2d %d %e %e\n",bayesian2D.get_rows(),totalLikelihood,sum);

    array_2d<index_t> used;
    used.set_name("output_used");
    used.set_dim(xct,yct);
    for(ix=0;ix<xct;ix++){
        for(iy=0;iy<yct;iy++){
            used.set(ix,iy,0);
        }
    }

    sum=0.0;
    row=0;
    for(i=0;i<chi_vals_sorted.get_dim() && sum<confidence_limit*totalLikelihood;i++){
        sum+=likelihood_sorted.get_data(i);
        true_dex=dexes.get_data(i);

        ix=coordinates.get_data(true_dex,xdex);
        iy=coordinates.get_data(true_dex,ydex);

        if(used.get_data(ix,iy)==0){
            bayesianFullD.set(row,0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
            bayesianFullD.set(row,1,_min.get_data(ydex)+iy*_dx.get_data(ydex));
            used.set(ix,iy,1);
            row++;
        }

    }

    printf("BayesianFullD %d %e %e\n",bayesianFullD.get_rows(),totalLikelihood,sum);

    value_t effective_chi;
    for(ix=0;ix<xct;ix++){
        trial.set(0,_min.get_data(xdex)+ix*_dx.get_data(xdex));
        for(iy=0;iy<yct;iy++){
            trial.set(1,_min.get_data(ydex)+iy*_dx.get_data(ydex));

            effective_chi=-2.0*log(marginalized_likelihood.get_data(ix,iy))-min_chi2_eff;
            if(effective_chi<0.0){
                printf("WARNING eff %e min %e\n",-2.0*log(marginalized_likelihood.get_data(ix,iy)),min_chi2_eff);
                exit(1);
            }
            if(effective_chi<=delta_chi_2d){
                frequentist2Deff.add_row(trial);
            }


        }
    }

    printf("freq2deff %d %e\n",frequentist2Deff.get_rows(),min_chi2_eff);

    ///write outputs
    FILE *output;
    char outname[2*letters];

    //heat map
    sprintf(outname,"%s_heatmap.txt",outname_root);
    output=fopen(outname,"w");
    fprintf(output,"#x y min_chi^2 marginalized_likelihood\n");
    for(ix=0;ix<minChi2Grid.get_rows();ix++){
        for(iy=0;iy<minChi2Grid.get_cols();iy++){
            if(marginalized_likelihood.get_data(ix,iy)>1.0e-10*max_likelihood ||
               minChi2Grid.get_data(ix,iy)<2.0*exception_value){

                fprintf(output,"%e %e %e %e\n",
                        _min.get_data(xdex)+ix*_dx.get_data(xdex),
                        _min.get_data(ydex)+iy*_dx.get_data(ydex),
                        minChi2Grid.get_data(ix,iy),
                        marginalized_likelihood.get_data(ix,iy));
           }
        }
    }
    fclose(output);


    array_2d<value_t> boundary;
    boundary.set_name("output_boundary");

    //points wwith chisq < a priori chi^2_lim
    convert_to_boundary(frequentistFullD,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentistFullD.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open frequentistFullD\n");
        exit(1);
    }
    fprintf(output,"#chi^2 <= %e\n",chi_lim_freq);
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    //points with chisq < chisq_min + delta_chisq(dof)
    sprintf(outname,"%s_frequentistFullDrelative_scatter.txt",outname_root);
    output=fopen(outname,"w");
    fprintf(output,"#chi^2-_chimin <= %e; chimin %e\n",delta_chi_bayes,_chi_min);
    for(ix=0;ix<frequentistFullDrelative.get_rows();ix++){
        fprintf(output,"%e %e\n",frequentistFullDrelative.get_data(ix,0),frequentistFullDrelative.get_data(ix,1));
    }
    fclose(output);

    convert_to_boundary(frequentistFullDrelative,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentistFullDrelative.txt",outname_root);
    output=fopen(outname,"w");
    fprintf(output,"#chi^2-_chimin <= %e; %e\n",delta_chi_bayes,_chi_min);
    if(output==NULL){
        printf("WARNING could not open freqFullDrel\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    //point with chisq < chisq_min + 6
    convert_to_boundary(frequentist2D,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentist2D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open freq2d\n");
        exit(1);
    }
    fprintf(output,"#chi^2 <= chimin+ %e; _chi_min %e\n",delta_chi_2d, _chi_min);

    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);

    //points wiht -2 ln(marginalized_likelihood) < -2.0*ln(max_marginalized_likelihood) + 6
    convert_to_boundary(frequentist2Deff,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_frequentist2Deff.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING coudl not open freq2deff\n");
        exit(1);
    }
    fprintf(output,"#-2 ln(marginalized likelihood) - min_chi2_eff<=%e ; %e\n",delta_chi_2d,min_chi2_eff);
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);


    //sum of marginalized likelihood < _confidence_limit * total
    convert_to_boundary(bayesian2D,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_bayesian2D.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open bayes2d\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);


    //sum of individual likelihood in full D space < _confidence_limit * total
    convert_to_boundary(bayesianFullD,_dx.get_data(xdex),_dx.get_data(ydex),boundary);

    sprintf(outname,"%s_bayesianFullD.txt",outname_root);
    output=fopen(outname,"w");
    if(output==NULL){
        printf("WARNING could not open bayesFullD\n");
        exit(1);
    }
    for(ix=0;ix<boundary.get_rows();ix++){
        fprintf(output,"%e %e\n",
        boundary.get_data(ix,0),boundary.get_data(ix,1));
    }
    fclose(output);


    printf("2d ct %d\n",xct*yct);
    printf("min found max found\n");
    printf("%e %e %e %e\n",_min.get_data(xdex),found_min.get_data(0),_max.get_data(xdex),found_max.get_data(0));
    printf("%e %e %e %e\n",_min.get_data(ydex),found_min.get_data(1),_max.get_data(ydex),found_max.get_data(1));
    printf("confidence limit %e\n",confidence_limit);
    printf("delta_chi_bayes %e\n",delta_chi_bayes);
    printf("chi_lim_freq %e\n",chi_lim_freq);

}

void control_integrator::find_chi_min(array_1d<value_t> &min,
                                      array_1d<value_t> &max){
    array_1d<value_t> center;
    array_2d<value_t> seed;
    Ran dice(42);

    seed.set_cols(min.get_dim());

    index_t ix,iy;
    for(ix=0;ix<_min.get_dim()+1;ix++){
        for(iy=0;iy<_min.get_dim();iy++){
            seed.set(ix,iy,min.get_data(iy)+dice.doub()*(max.get_data(iy)-min.get_data(iy)));
        }
    }


    simplex_minimizer f_min;
    f_min.set_chisquared(_chisq);
    f_min.set_minmax(_min,_max);
    f_min.set_dice(&dice);
    f_min.use_gradient();
    f_min.find_minimum(seed, center);

    value_t trial_min=f_min.get_minimum();
    printf("simplex found %e\n",trial_min);
    if(trial_min<_chi_min){
        _chi_min=trial_min;
    }
}
