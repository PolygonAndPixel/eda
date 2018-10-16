#include "cost_fn.h"

cost_fn::cost_fn(chisq_wrapper *cc, array_1d<index_t> &aa){
    _envelope=1.0;
    build(cc,aa,1);
}

cost_fn::cost_fn(chisq_wrapper *cc, array_1d<index_t> &aa, index_t min_or_med){
    _envelope=1.0;
    build(cc,aa,min_or_med);
}

void cost_fn::build(chisq_wrapper *cc, array_1d<index_t> &aa, index_t min_or_med){

    // printf("building cost_fn with %d associates\n",aa.get_dim());
    _called=0;

    _median_associate.set_name("dchi_interior_median");
    _pt_cache.set_name("dchi_interior_pt_cache");
    _fn_cache.set_name("dchi_interior_fn_cache");
    _associates.set_name("dchi_interior_fn_associates");

    _associates.reset_preserving_room();
    _fn_cache.reset_preserving_room();
    _pt_cache.reset_preserving_room();
    _median_associate.reset_preserving_room();

    _chifn=cc;

    index_t i;

    for(i=0;i<aa.get_dim();i++){
        _associates.add(aa.get_data(i));
    }

    array_1d<value_t> norm;
    array_1d<value_t> min,max;
    index_t j;
    for(i=0;i<_associates.get_dim();i++){
        for(j=0;j<_chifn->get_dim();j++){
            if(i==0 || _chifn->get_pt(_associates.get_data(i),j)<min.get_data(j)){
                min.set(j,_chifn->get_pt(_associates.get_data(i),j));
            }
            if(i==0 || _chifn->get_pt(_associates.get_data(i),j)>max.get_data(j)){
                max.set(j,_chifn->get_pt(_associates.get_data(i),j));
            }
        }
    }

    if(min.get_dim()>0){
        for(i=0;i<_chifn->get_dim();i++){
            norm.set(i,max.get_data(i)-min.get_data(i));
        }
    }

    for(i=0;i<_chifn->get_dim();i++){
        _median_associate.set(i,0.5*(max.get_data(i)+min.get_data(i)));
    }

    array_1d<value_t> norm_sorted;
    array_1d<index_t> norm_dex;
    for(i=0;i<norm.get_dim();i++){
        norm_dex.add(i);
    }
    sort(norm, norm_sorted, norm_dex);
    if(min_or_med==1){
        _scalar_norm=norm_sorted.get_data(norm_dex.get_dim()/2);
    }
    else{
        _scalar_norm=norm_sorted.get_data(0);
    }
}


value_t cost_fn::nn_distance(const array_1d<value_t> &pt){
    value_t dd;
    index_t i,j,k;
    value_t dd_avg;
    value_t ct=0.0;
    dd_avg=0.0;

    for(i=0;i<_associates.get_dim();i++){
        dd=0.0;
        for(j=0;j<_chifn->get_dim();j++){
            dd+=power((pt.get_data(j)-_chifn->get_pt(_associates.get_data(i),j))/_scalar_norm,2);
        }
        if(dd>1.0e-20){
            dd_avg += 1.0/sqrt(dd);
            ct+=1.0;
        }
    }
    if(ct>0.0){
        dd_avg = dd_avg/ct;
        return 1.0/dd_avg;
    }
    return 0.0;
}


index_t cost_fn::get_called(){
   return _called;
}

value_t cost_fn::operator()(const array_1d<value_t> &pt){

    if(_chifn==NULL){
        printf("WARNING cannot call cost_fn operator; _chifn is NULL\n");
        exit(1);
    }

    _called++;

    value_t mu;
    index_t i_found;
    _chifn->evaluate(pt,&mu,&i_found);

    if(_associates.get_dim()==0){
        return mu;
    }

    value_t delta=_chifn->target()-_chifn->chimin();

    value_t distance;

    value_t exp_term;

    if(mu>_chifn->target()){
        exp_term=exp((_chifn->target()-mu)/_envelope);
    }
    else{
        exp_term=1.0;
    }

    if(exp_term>1.0e-7){
        distance=nn_distance(pt);
    }
    else{
        distance=0.0;
    }

    value_t val;
    val = mu-1.0*delta*distance*exp_term;
    _pt_cache.add(i_found);
    _fn_cache.add(val);
    return val;
}
