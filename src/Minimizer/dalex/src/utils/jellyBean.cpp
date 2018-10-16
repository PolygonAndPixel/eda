#include "jellyBean.h"

////////////////////////////////jellyBeanData

chiSquaredData::chiSquaredData(index_t dd, index_t cc, value_t ww, index_t nData, value_t sigma) : chisquared(dd, cc, ww){
    printf("_dim %d\n",_dim);

    if(_dim<4){
        printf("WARNING must have at least 4 params in chisquaredData\n");
        exit(1);
    }

    if(_dim%4!=0){
        printf("WARNING dimensionality must be multiple of four\n");
        exit(1);
    }

    _x_values.set_name("jellyBeanData_x");
    _y_values.set_name("jellyBeanData_y");
    _sigma.set_name("jellyBeanData_sigma");
    _wave_phase.set_name("jellyBeanData_phase");
    _wave_lambda.set_name("jellyBeanData_lambda");
    _wave_amp.set_name("jellBeanData_amp");
    _env_d_center.set_name("jellyBeanData_env_d_center");
    _env_width.set_name("jellyBeanData_env_width");
    _param_buffer.set_name("jellyBeanData_param_buffer");

    Ran local_dice(nData+dd+cc);

    _ndata=nData;
    _sig=sigma;
    _xmax=3.0;
    _dx=0.03;

    index_t ix;
    value_t ll;
    value_t ll_min,ll_max;
    ll_min=log(0.05*_xmax);
    ll_max=log(0.6*_xmax);

    value_t last_center;
    value_t d_center;
    value_t d_c;
    d_center=_xmax/value_t(_dim/4);
    last_center=0.0;
    for(ix=0;ix<_dim;ix+=4){
        _wave_phase.add(local_dice.doub()*_xmax);

        ll=local_dice.doub()*(ll_max-ll_min)+ll_min;
        _wave_lambda.add(exp(ll));

        _wave_amp.add(local_dice.doub()*5.0);

        d_c=local_dice.doub()*d_center;
        _env_d_center.add(d_c);

        _env_width.add((local_dice.doub()*0.5+0.05)*_xmax);
    }

}

void chiSquaredData::set_width(index_t ic, index_t id, value_t dd){
    _widths.set(ic, id, dd);
}

void chiSquaredData::write_data(){
    FILE *output;
    output=fopen("data_scratch.txt", "w");
    index_t ix;
    for(ix=0;ix<_x_values.get_dim();ix++){
        fprintf(output,"%e %e %e\n",_x_values.get_data(ix),_y_values.get_data(ix),_sigma.get_data(ix));
    }
    fclose(output);
}

void chiSquaredData::print_mins(){
    array_1d<value_t> trial,params;
    trial.set_name("jellyBeanData_print_mins_trial");
    params.set_name("jellyBeanData_print_mins_params");
    index_t ic,ix;
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            trial.set(ix,_centers.get_data(ic,ix));
        }
        printf("center %d val %e\n",ic,this[0](trial));
        convert_params(trial,params,ic);
        for(ix=0;ix<_dim;ix++){
            printf("    %e %e -- width %e\n",
            _centers.get_data(ic,ix),
            params.get_data(ix),_widths.get_data(0,ix));
        }
    }

    printf("\nbases\n");
    for(ix=0;ix<_dim;ix++){
        for(ic=0;ic<_dim;ic++){
            printf("%.3e ",_bases.get_data(ic,ix));
        }
        printf("\n");
    }

    printf("\nparams\n");
    for(ix=0;ix<_dim/4;ix++){
        printf("amp %e\n",_wave_amp.get_data(ix));
        printf("phase %e\n",_wave_phase.get_data(ix));
        printf("lambda %e\n",_wave_lambda.get_data(ix));
        printf("env_d_x %e\n",_env_d_center.get_data(ix));
        printf("env_width %e\n",_env_width.get_data(ix));
        printf("\n");
    }

}

void chiSquaredData::initialize_data(){

    array_1d<value_t> params;
    params.set_name("jellyBeanData_initiailize_params");
    index_t ix,ic;

    for(ix=0;ix<_dim;ix++){
        if(ix%4==2){
            params.set(ix,1.0);
        }
        else{
            params.set(ix,0.0);
        }
    }

    array_1d<value_t> samples;
    samples.set_name("jellyBeanData_initialize_samples");

    value_t xx,yy,mean,var;
    for(ix=0, xx=0.0;xx<_xmax;xx+=_dx, ix++){
        mean=data_function(params,xx);
        for(ic=0;ic<_ndata;ic++){
            samples.set(ic,normal_deviate(_dice,mean,_sig));
        }
        mean=0.0;
        for(ic=0;ic<samples.get_dim();ic++){
            mean+=samples.get_data(ic);
        }
        mean=mean/value_t(samples.get_dim());
        _x_values.set(ix,xx);
        _y_values.set(ix,mean);

        var=0.0;
        for(ic=0;ic<samples.get_dim();ic++){
            var+=power(samples.get_data(ic)-mean,2);
        }
        var=var/value_t((samples.get_dim()-1)*samples.get_dim());
        _sigma.set(ix,sqrt(var));
    }

    write_data();

}

void chiSquaredData::convert_params(const array_1d<value_t> &pt, array_1d<value_t> &out, index_t ic){
    printf("Called void convert_params");
    exit(1);
}

value_t chiSquaredData::data_function(array_1d<value_t> &params, value_t xx){

    value_t ans=0.0;

    value_t envelope;
    value_t env_x;
    value_t wave;
    value_t amp;
    value_t lambda;
    value_t phase;
    value_t last_center=0.0;
    value_t center;
    index_t ix,i_param;
    for(ix=0,i_param=0;ix<_dim;ix+=4,i_param++){
        center=last_center+fabs(params.get_data(ix))+_env_d_center.get_data(i_param);
        env_x=(xx-center)/(params.get_data(ix+1)+_env_width.get_data(i_param));

        envelope=0.1+exp(-0.5*power(env_x,2));
        last_center=center;

        amp=params.get_data(ix+2)+_wave_amp.get_data(i_param);
        lambda=params.get_data(ix+3)*0.5*_xmax+_wave_lambda.get_data(i_param);
        wave=sin(2.0*pi*(xx-_wave_phase.get_data(i_param))/lambda);
        ans+=envelope*amp*wave;
    }

    return ans;

}

value_t chiSquaredData::operator()(const array_1d<value_t> &pt){

    value_t before=value_t(time(NULL));

    if(_x_values.get_dim()==0){
        initialize_data();
    }

    value_t chisq,chisq_min;
    value_t yy;
    index_t ic,ix;
    for(ic=0;ic<_ncenters;ic++){
        chisq=0.0;
        convert_params(pt,_param_buffer, ic);
        for(ix=0;ix<_x_values.get_dim();ix++){
            yy=data_function(_param_buffer,_x_values.get_data(ix));
            chisq+=power((yy-_y_values.get_data(ix))/_sigma.get_data(ix),2);
        }

        if(ic==0 || chisq<chisq_min){
            chisq_min=chisq;
        }
    }

    _called++;
    _time_spent+=value_t(time(NULL))-before;

    if(isnan(chisq_min)){
        chisq_min=2.0*exception_value;
    }

    if(_with_logging==1){
        _log_point(pt, chisq_min);
    }

    return chisq_min;
}


jellyBeanData::jellyBeanData(index_t dd, index_t cc, value_t wc, index_t nData, value_t sigma) :
chiSquaredData(dd, cc, wc, nData, sigma){

    _parabola_centers.set_name("parabola_centers");
    _parabola_centers.set_cols(2);
    _parabola_x.set_name("parabola_x");
    _parabola_x.set_cols(2);
    _parabola_y.set_name("parabola_y");
    _parabola_y.set_cols(2);

    index_t ic,ix,iy;
    array_1d<value_t> trial;
    trial.set_name("jellyBean_data_constructor_trial");
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<2;ix++){
            _parabola_centers.set(ic,ix,-20.0+_dice->doub()*40.0);
        }

        for(ix=0;ix<2;ix++){
            trial.set(ix,normal_deviate(_dice,0.0,1.0));
        }
        trial.normalize();
        _parabola_x.set(ic,0,trial.get_data(0));
        _parabola_x.set(ic,1,trial.get_data(1));
        _parabola_y.set(ic,0,trial.get_data(1));
        _parabola_y.set(ic,1,-1.0*trial.get_data(0));

        // in projected coordinates
        for(ix=0;ix<_dim;ix++){
            _centers.set(ic,ix,-20.0+_dice->doub()*40.0);
        }

    }

    _parabola_curvature=4.0;

    /*index_t i,j;
    for(i=0;i<_dim;i++){
        for(j=0;j<_dim;j++){
            if(i==j){
                _bases.set(i,j,1.0);
            }
            else{
                _bases.set(i,j,0.0);
            }
        }
    }

    _parabola_x.set(0,0,1.0);
    _parabola_x.set(0,1,0.0);
    _parabola_y.set(0,0,0.0);
    _parabola_y.set(0,1,1.0);*/

}


void jellyBeanData::convert_params(const array_1d<value_t> &pt_in, array_1d<value_t> &out, index_t ic){

    array_1d<value_t> pt;
    pt.set_name("convert_params_pt");
    project_to_basis(pt_in,pt);

    value_t x_is,y_is;
    x_is=(pt.get_data(0)-_parabola_centers.get_data(ic,0))*_parabola_x.get_data(ic,0)
         +(pt.get_data(1)-_parabola_centers.get_data(ic,1))*_parabola_x.get_data(ic,1);

    y_is=(pt.get_data(0)-_parabola_centers.get_data(ic,0))*_parabola_y.get_data(ic,0)
         +(pt.get_data(1)-_parabola_centers.get_data(ic,1))*_parabola_y.get_data(ic,1);

    value_t rr=sqrt(x_is*x_is+y_is*y_is);

    value_t cos_theta,sin_theta;
    cos_theta=x_is/rr;
    sin_theta=y_is/rr;

    if(fabs(cos_theta)>1.0){
        if(fabs(cos_theta)>1.001){
            printf("WARNING somehow cos_theta>1.0 -- %e\n",cos_theta);
            exit(1);
        }
    }

    if(fabs(sin_theta)>1.0){
        if(fabs(sin_theta)>1.001){
            printf("WARNING somehow sin_theta>1.0 -- %e\n",sin_theta);
            exit(1);
        }
    }

    value_t r_shldbe;

    if(fabs(cos_theta)<1.0e-5 && sin_theta>0.0){
        r_shldbe=1.0e10/_parabola_curvature;
    }
    else if(fabs(cos_theta)<1.0e-10 && sin_theta<0.0){
        r_shldbe=0.25/_parabola_curvature;
    }
    else{
        r_shldbe=(1.0+sin_theta)/(2.0*_parabola_curvature*cos_theta*cos_theta);
    }

    value_t d_radius;
    d_radius=fabs(rr-r_shldbe);

    value_t y_distance;
    y_distance=(y_is+0.25/_parabola_curvature)/_widths.get_data(ic,0);
    value_t y_term=0.2*sin(y_distance)+log(0.5*y_distance+1.0);
    out.set(0,y_term);

    value_t x_shldbe,dx;
    if(y_distance<0.0){
        out.set(1,d_radius/_widths.get_data(ic,1));
    }
    else{
        x_shldbe=sqrt((y_is+0.25/_parabola_curvature)/_parabola_curvature);
        if(fabs(x_is+x_shldbe)<fabs(x_is-x_shldbe)){
            dx=fabs(x_is+x_shldbe);
        }
        else{
            dx=fabs(x_is-x_shldbe);
        }

        out.set(1,dx/_widths.get_data(ic,1));
    }

    if(out.get_data(1)<0.0){
        printf("WARNING out 1 %e\n",out.get_data(1));
        exit(1);
    }

    index_t ix;
    for(ix=2;ix<_dim;ix++){
        if(ix%4==2){
            out.set(ix,exp((pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix)));
        }
        else  if(ix%4==3){
            if(pt.get_data(ix)>_centers.get_data(ic,ix)){
                out.set(ix,log(1.0+(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix)));
            }
            else{
                out.set(ix,0.3*(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix));
            }
        }
        else{
            out.set(ix,(pt.get_data(ix)-_centers.get_data(ic,ix))/_widths.get_data(ic,ix));
        }
    }

    value_t xx;
    for(ix=0;ix<_dim;ix++){
        if(ix%4!=2){
            out.multiply_val(ix,0.01);
        }
        else{
            xx=out.get_data(ix);
            if(xx>0.1){
                out.set(ix,1.0+0.1*(xx-1.0));
            }
        }
    }


}

//////////////////////ellipse classes////////////

ellipseData::ellipseData(index_t dd, index_t cc, index_t nData, value_t sigma) :
chiSquaredData(dd, cc, 1.0, nData, sigma){

    _dir.set_name("ellipseData_dir");
    _projected.set_name("ellipseData_projected");

    index_t ix,ic;
    for(ic=0;ic<_ncenters;ic++){
        for(ix=0;ix<_dim;ix++){
            _centers.set(ic,ix,-30.0+60.0*_dice->doub());
        }
    }
}

void ellipseData::convert_params(const array_1d<value_t> &pt, array_1d<value_t> &out, index_t ic){

    index_t ix;
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_centers.get_data(ic,ix));
    }

    project_to_basis(_dir,_projected);

    for(ix=0;ix<_dim;ix++){
        out.set(ix,_projected.get_data(ix)/_widths.get_data(ic,ix));
    }
}

nonGaussianEllipseData::nonGaussianEllipseData(index_t i1, index_t i2, index_t i3, value_t d1) :
ellipseData(i1, i2, i3, d1){}

void nonGaussianEllipseData::convert_params(const array_1d<value_t> &pt, array_1d<value_t> &out, index_t ic){

    index_t ix;
    for(ix=0;ix<_dim;ix++){
        _dir.set(ix,pt.get_data(ix)-_centers.get_data(ic,ix));
    }

    project_to_basis(_dir,_projected);


    for(ix=1;ix<_dim;ix++){
        out.set(ix,_projected.get_data(ix)/_widths.get_data(ic,ix));
    }

    value_t mu,xx;
    xx=_projected.get_data(0)/_widths.get_data(ic,0);
    if(fabs(xx)<fabs(xx-0.05)){
        out.set(0,xx);
    }
    else{
        out.set(0,xx-0.05);
    }

}
