#include "chisq_wrapper.h"

chisq_wrapper::chisq_wrapper(){
    _kptr=NULL;
    _dice=NULL;
    _chifn=NULL;
    _fn.set_name("chisq_wrapper_fn");
    _characteristic_length.set_name("chisq_wrapper_characteristic_length");
    _range_min.set_name("chisq_wrapper_range_min");
    _range_max.set_name("chisq_wrapper_range_max");
    _valid_neigh.set_name("chisq_wrapper_valid_neigh");
    _valid_dd.set_name("chisq_wrapper_valid_dd");
    _ddmin=1.0e-8;
    _adaptive_target=1;
    _deltachi=-1.0;
    _chimin=2.0*exception_value;
    _target=2.0*exception_value;
    _seed=-1;
    _called=0;
    _iWhere=-1;
    _expected_min=-1.0;
    _confidence_limit=-1.0;
    _expected_delta=-1.0;
    _outname[0]=0;
    _timingname[0]=0;
    _last_written=0;
    _write_every=50000;
    _time_started=value_t(time(NULL));
    _time_batch=value_t(time(NULL));
    _last_time_spent=0.0;
    _search_type=-1;
    _search_type_log.set_name("chisq_wrapper_search_type_log");
}

chisq_wrapper::~chisq_wrapper(){
    if(_kptr!=NULL){
        delete _kptr;
    }

    if(_dice!=NULL){
        delete _dice;
    }
}

void chisq_wrapper::set_confidence_limit(value_t cc){
    _confidence_limit=cc;
}

void chisq_wrapper::set_dof(index_t dd){
    _dof=dd;
}

void chisq_wrapper::set_chisquared(chisquared *xx){
    _chifn=xx;
}

index_t chisq_wrapper::get_seed(){
    return _seed;
}

void chisq_wrapper::set_seed(index_t i){
    _seed=i;
}

void chisq_wrapper::set_target(value_t xx){
    _target=xx;
    _adaptive_target=0;
}

void chisq_wrapper::set_ddmin(value_t dd){
    _ddmin=dd;
}

void chisq_wrapper::set_min(array_1d<value_t> &vv){
    if(_kptr!=NULL){
        printf("WARNING in chisq_wrapper setting min but kptr not null\n");
        exit(1);
    }

    index_t i;
    _range_min.reset();
    for(i=0;i<vv.get_dim();i++){
        _range_min.set(i,vv.get_data(i));
    }
}

void chisq_wrapper::set_max(array_1d<value_t> &vv){
    if(_kptr!=NULL){
        printf("WARNINGin chisq_wrapper setting max but kptr not null\n");
        exit(1);
    }

    index_t i;
    _range_max.reset();
    for(i=0;i<vv.get_dim();i++){
        _range_max.set(i,vv.get_data(i));
    }
}

void chisq_wrapper::set_characteristic_length(index_t dex, value_t xx){
    if(_kptr!=NULL){
        printf("WARNING in chisq_wrapper setting characteristic_length but kptr not null\n");
        exit(1);
    }

    index_t i;
    if(dex>_characteristic_length.get_dim()){
        for(i=_characteristic_length.get_dim();i<dex;i++){
            _characteristic_length.set(i,-1.0);
        }
    }

    _characteristic_length.set(dex,xx);
}

value_t chisq_wrapper::get_characteristic_length(index_t dex){
    if(dex<0 || dex>=_characteristic_length.get_dim() ||
       _characteristic_length.get_data(dex)<0.0){

           if(_chifn->get_max(dex)-_chifn->get_min(dex)<exception_value){
               return _chifn->get_max(dex)-_chifn->get_min(dex);
           }
           else{
               return 1.0;
           }

   }

   return _characteristic_length.get_data(dex);

}

void chisq_wrapper::set_deltachi(value_t xx){
    if(_adaptive_target!=1){
        printf("WARNING chisq_wrapper trying to set detlachi, but not an adaptive target\n");
        exit(1);
    }
    _deltachi=xx;
}

void chisq_wrapper::initialize(index_t npts){
    if(_chifn==NULL){
        printf("WARNING calling chisq_wrapper_initialize with null chifn\n");
        exit(1);
    }

    if(_kptr!=NULL){
        printf("WARNING calling chisq_wrapper_initialize even though kptr not null\n");
        exit(1);
    }

    if(_range_min.get_dim()!=_chifn->get_dim()){
        printf("WARNING chisq_wrapper_initialize dim %d range_min %d\n",
        _chifn->get_dim(),_range_min.get_dim());

        exit(1);
    }

    if(_range_max.get_dim()!=_chifn->get_dim()){
        printf("WARNING chisq_wrapper_initialize dim %d range_max %d\n",
        _chifn->get_dim(),_range_max.get_dim());

        exit(1);
    }

    if(_dice==NULL){
        if(_seed<0)_seed=index_t(time(NULL));
        _dice=new Ran(_seed);
    }

    array_2d<value_t> data;
    array_1d<value_t> vv;
    vv.set_name("chisq_wrapper_initialize_vv");
    data.set_name("chisq_wrapper_data");
    index_t i,j;
    value_t mu;

    _fn.reset();
    data.set_cols(_chifn->get_dim());
    for(i=0;i<npts;i++){
        mu=2.0*exception_value;

        while(mu>exception_value){
            for(j=0;j<_chifn->get_dim();j++){
                vv.set(j,_range_min.get_data(j)+_dice->doub()*(_range_max.get_data(j)-_range_min.get_data(j)));
            }
            mu=_chifn[0](vv);
            _called++;
        }

        if(mu<_chimin){
            _chimin=mu;
            _mindex=i;
        }
        _fn.add(mu);
        data.add_row(vv);
    }

    array_1d<value_t> temp_max,temp_min;
    temp_max.set_name("chisq_wrapper_initialize_temp_max");
    temp_min.set_name("chisq_wrapper_initialize_temp_min");
    for(i=0;i<_chifn->get_dim();i++){
        temp_min.set(i,0.0);
        temp_max.set(i,get_characteristic_length(i));
    }

    _kptr=new kd_tree(data,temp_min,temp_max);
    for(i=0;i<_kptr->get_pts();i++){
        _search_type_log.set(i,_type_init);
    }

    if(_adaptive_target==1){
        if(_deltachi<0.0){
            printf("WARNING when initializing chisq_wrapper deltachi %e\n",_deltachi);
            exit(1);
        }
        _target=_chimin+_deltachi;
    }
}

void chisq_wrapper::is_it_safe(char *word){
    if(_kptr==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("kptr is null\n");
        exit(1);
    }

    if(_dice==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("dice is null\n");
        exit(1);
    }

    if(_chifn==NULL){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("chifn is null\n");
        exit(1);
    }

    if(_adaptive_target==1 && _deltachi<0.0){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("adaptive target but deltachi %e\n",_deltachi);
        exit(1);
    }

    if(_fn.get_dim()!=_kptr->get_pts()){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("fn dim %d kptr pts %d\n",_fn.get_dim(),_kptr->get_pts());
        exit(1);
    }

    if(_chifn->get_dim()!=_kptr->get_dim()){
        printf("WARNING in chisq_wrapper::%s\n",word);
        printf("chifn dim %d kptr %d\n",_chifn->get_dim(),_kptr->get_dim());
        exit(1);
    }
}

value_t chisq_wrapper::target(){
    return _target;
}

value_t chisq_wrapper::chimin(){
    return _chimin;
}

index_t chisq_wrapper::mindex(){
    return _mindex;
}

value_t chisq_wrapper::get_deltachi(){
    return _deltachi;
}

index_t chisq_wrapper::in_bounds(index_t dex, value_t val){
    is_it_safe("in_bounds");
    if(dex<0 || dex>_chifn->get_dim()){
        printf("WARNING asked for in_bounds on dex %d but %d\n",dex,_chifn->get_dim());
    }

    if(_chifn->get_max(dex)>-1.0*exception_value && val>_chifn->get_max(dex)){
        return 0;
    }

    if(_chifn->get_min(dex)<exception_value && val<_chifn->get_min(dex)){
        return 0;
    }

    return 1;
}

index_t chisq_wrapper::in_bounds(const array_1d<value_t> &pt){
    is_it_safe("in_bounds");

    index_t i;
    for(i=0;i<pt.get_dim();i++){
        if(_chifn->get_max(i)>-1.0*exception_value && pt.get_data(i)>_chifn->get_max(i)){
            return 0;
        }
        if(_chifn->get_min(i)<exception_value && pt.get_data(i)<_chifn->get_min(i)){
            return 0;
        }
    }
    return 1;

}

index_t chisq_wrapper::is_valid(const array_1d<value_t> &pt, index_t *neighdex){
    is_it_safe("is_valid");

    neighdex[0]=-1;

    if(in_bounds(pt)==0){
        return 0;
    }

    _kptr->nn_srch(pt,1,_valid_neigh,_valid_dd);
    if(_valid_dd.get_data(0)<_ddmin){
        neighdex[0]=_valid_neigh.get_data(0);
        return 0;
    }

    return 1;

}

value_t chisq_wrapper::operator()(const array_1d<value_t> &pt){
    value_t mu;
    index_t dex;
    index_t i;
    evaluate(pt,&mu,&dex);
    return mu;
}

value_t chisq_wrapper::raw_evaluate(const array_1d<value_t> &pt){
    return _chifn[0](pt);
}

void chisq_wrapper::evaluate(const array_1d<value_t> &pt, value_t *value, index_t *dex){
    is_it_safe("evaluate");

    index_t validity,neighdex;
    validity=is_valid(pt,&neighdex);

    dex[0]=-1;
    if(validity!=1){
        if(neighdex>=0){
            value[0]=_fn.get_data(neighdex);
            dex[0]=neighdex;
        }
        else{
            value[0]=2.0*exception_value;
            dex[0]=-1;
        }
        return;
    }

    value_t mu;
    mu=_chifn[0](pt);
    value[0]=mu;
    _called++;

    if(mu<exception_value){
        _kptr->add(pt);
        _fn.add(mu);
        _search_type_log.add(_search_type);
        dex[0]=_kptr->get_pts()-1;
    }

    if(mu<_chimin){
        _chimin=mu;
        _mindex=_kptr->get_pts()-1;
        if(_adaptive_target==1){
            _target=_chimin+_deltachi;
        }
    }

    if(get_pts()-_last_written>_write_every){
        write_pts();
    }

}

index_t chisq_wrapper::get_called(){
    is_it_safe("get_called");
    return _chifn->get_called();
}

value_t chisq_wrapper::get_time_spent(){
    is_it_safe("get_time_spent");
    return _chifn->get_time_spent();
}

index_t chisq_wrapper::get_pts(){
    is_it_safe("get_pts");
    return _kptr->get_pts();
}

index_t chisq_wrapper::get_dim(){
    is_it_safe("get_dim");
    return _kptr->get_dim();
}

value_t chisq_wrapper::random_double(){
    if(_dice==NULL){
        printf("WARNING chisq_wrapper random_double dice is null\n");
        exit(1);
    }

    return _dice->doub();
}

index_t chisq_wrapper::random_int(){
    if(_dice==NULL){
        printf("WARNING chisq_wrapper random_int dice is null\n");
        exit(1);
    }

    return _dice->int32();
}

value_t chisq_wrapper::get_fn(index_t dex){
    if(dex<0 || dex>=_fn.get_dim()){
        printf("WARNING asking for fn %d but %d\n",dex,_fn.get_dim());
    }

    return _fn.get_data(dex);
}


value_t chisq_wrapper::get_pt(index_t dex, index_t idim){
    is_it_safe("get_pt");

    if(dex<0 || dex>=_fn.get_dim()){
        printf("WARNING asking for pt %d but only have %d \n",
        dex,_fn.get_dim());

        exit(1);
    }

    if(idim<0 || idim>=_kptr->get_dim()){
        printf("WARNING asking for pt dim %d but only have %d\n",
        idim,_kptr->get_dim());

        exit(1);
    }

    return _kptr->get_pt(dex,idim);
}

array_1d<value_t> chisq_wrapper::get_pt(index_t dex){
    is_it_safe("get_pt(dex)");
    if(dex<0 || dex>=_fn.get_dim()){
        printf("WARNING asking wrapper for pt %d but only have %d\n",
        dex,_fn.get_dim());
        exit(1);
    }

    return _kptr->get_pt(dex);
}

void chisq_wrapper::nn_srch(array_1d<value_t> &vv, index_t kk, array_1d<index_t> &neigh, array_1d<value_t> &dd){
    is_it_safe("nn_srch");
    _kptr->nn_srch(vv,kk,neigh,dd);
}

Ran* chisq_wrapper::get_dice(){
    is_it_safe("get_dice");
    return _dice;
}

kd_tree* chisq_wrapper::get_tree(){
    is_it_safe("get_tree");
    return _kptr;
}

array_1d<value_t>* chisq_wrapper::get_fn_arr(){
    is_it_safe("get_fn_arr");
    return &_fn;
}

value_t chisq_wrapper::get_min(index_t dex){
    if(dex<0 or dex>=_range_min.get_dim()){
        printf("WARNING asked chisq_wrapper for min %d of %d\n",
        dex,_range_min.get_dim());

        exit(1);
    }

    return _range_min.get_data(dex);
}

value_t chisq_wrapper::get_max(index_t dex){
    if(dex<0 or dex>=_range_max.get_dim()){
        printf("WARNING asked chisq_wrapper for max %d of %d\n",
        dex,_range_max.get_dim());

        exit(1);
    }

    return _range_max.get_data(dex);
}

void chisq_wrapper::get_min(array_1d<value_t> &vv){
    index_t i;
    vv.set_dim(_range_min.get_dim());
    for(i=0;i<_range_min.get_dim();i++){
        vv.set(i,_range_min.get_data(i));
    }
}

void chisq_wrapper::get_max(array_1d<value_t> &vv){
    index_t i;
    vv.set_dim(_range_max.get_dim());
    for(i=0;i<_range_max.get_dim();i++){
        vv.set(i,_range_max.get_data(i));
    }
}

value_t chisq_wrapper::distance(array_1d<value_t> &p1, array_1d<value_t> &p2){
    is_it_safe("distance(arr,arr)");
    return _kptr->distance(p1,p2);
}

value_t chisq_wrapper::distance(array_1d<value_t> &p, index_t dex){
    is_it_safe("distance(arr,i)");
    return _kptr->distance(p,dex);
}

value_t chisq_wrapper::distance(index_t i1, index_t i2){
    is_it_safe("distance(i,i)");
    return _kptr->distance(i1,i2);
}

void chisq_wrapper::find_gradient(array_1d<value_t> &pt, array_1d<value_t> &grad){
    is_it_safe("find_gradient");

    value_t fCenter;
    index_t iCenter;
    evaluate(pt,&fCenter,&iCenter);

    index_t ix,i;
    value_t x2,dx;
    array_1d<value_t> trial;
    trial.set_name("chisq_wrapper_gradient_trial");

    for(i=0;i<_kptr->get_dim();i++){
        trial.set(i,pt.get_data(i));
    }

    for(ix=0;ix<_kptr->get_dim();ix++){
        if(_range_max.get_dim()==0 || _range_max.get_data(ix)-_range_min.get_data(ix)>1.0){
            dx=0.01;
        }
        else{
            dx=0.01*(_range_max.get_data(ix)-_range_min.get_data(ix));
        }

        trial.add_val(ix,dx);
        evaluate(trial,&x2,&i);
        grad.set(ix,(fCenter-x2)/(pt.get_data(ix)-trial.get_data(ix)));
        trial.set(ix,pt.get_data(ix));
    }

}

void chisq_wrapper::copy(chisq_wrapper &in){
    if(in._chifn==NULL){
        printf("WARNING cannot copy chisq_wrapper chifn is null\n");
        exit(1);
    }

    if(in._kptr==NULL){
        printf("WARNING cannot copy chisq_wrapper kptr is null\n");
        exit(1);
    }

    printf("copying chisq\n");

    if(_kptr!=NULL){
        delete _kptr;
        _kptr=NULL;
    }
    _chifn=in._chifn;
    _chimin=in._chimin;
    _target=in._target;
    _deltachi=in._deltachi;
    _ddmin=in._ddmin;
    _adaptive_target=in._adaptive_target;
    _seed=in._seed;
    _mindex=in._mindex;


    index_t i,j;

    _range_max.reset();
    _range_min.reset();
    _characteristic_length.reset();
    for(i=0;i<in._characteristic_length.get_dim();i++){
        _characteristic_length.set(i,in._characteristic_length.get_data(i));
    }
    for(i=0;i<in._range_min.get_dim();i++){
        _range_min.set(i,in._range_min.get_data(i));
    }
    for(i=0;i<in._range_max.get_dim();i++){
        _range_max.set(i,in._range_max.get_data(i));
    }

    _fn.reset();

    for(i=0;i<in.get_pts();i++){
        _fn.set(i,in.get_fn(i));
    }

    printf("making kptr\n");
    _kptr=new kd_tree(in._kptr[0]);
    if(_seed<0)_seed=index_t(time(NULL));
    _dice=new Ran(_seed);
    printf("done copying\n");
}

index_t chisq_wrapper::could_it_go_lower(value_t chimin){

    value_t tol=1.0;

    if(chimin<tol){
        return 0;
    }

    if(_adaptive_target==1){
        if(_expected_min<0.0 && _dof>0){
            _expected_min=_distribution.maximum_likelihood_chisquared_value(value_t(_dof));
        }

        if(_expected_min<0.0){
            return 1;
        }
        else{
            if(chimin<_expected_min+tol){
                return 0;
            }
            else{
                return 1;
            }
        }
    }
    else{
        if(_expected_delta<0.0 && _confidence_limit>0.0){
            _expected_delta=_distribution.confidence_limit(value_t(_kptr->get_dim()),_confidence_limit);
        }

        if(_expected_delta<0.0){
            return 1;
        }
        else{
            if(chimin<_target-_expected_delta+tol){
                return 0;
            }
            else{
                return 1;
            }

        }
    }

    return 1;

}

void chisq_wrapper::write_pts(){
    if(_outname[0]==0 || _timingname[0]==0){
        printf("CANNOT write points, timingname or outname not set\n");
        printf("outname: %s\n",_outname);
        printf("timingname: %s\n",_timingname);
        exit(1);
    }

    if(_search_type_log.get_dim() != get_pts()){
        printf("WARNING %d pts but %d log\n",
        get_pts(),_search_type_log.get_dim());
        exit(1);
    }

    FILE *output;

    index_t this_batch=get_pts()-_last_written;

    index_t i,j;
    output=fopen(_outname,"w");
    fprintf(output,"# ");
    for(i=0;i<get_dim();i++){
        fprintf(output,"p%d ",i);
    }
    fprintf(output,"chisq log\n");
    for(i=0;i<get_pts();i++){
        for(j=0;j<get_dim();j++){
            fprintf(output,"%.18e ",get_pt(i,j));
        }
        fprintf(output,"%.18e %d\n",get_fn(i),_search_type_log.get_data(i));
    }
    fclose(output);

    if(_last_written==0){
        output=fopen(_timingname,"w");
        fprintf(output,"#seed %d\n",get_seed());
    }
    else{
        output=fopen(_timingname,"a");
    }

    fprintf(output,"%d %d min %.4e target %.4e -- timing -- %.4e %.4e -- %.4e %.4e -- overhead %.4e",
        get_pts(),
        get_called(),
        chimin(),
        target(),
        value_t(time(NULL))-_time_started,
        (value_t(time(NULL))-_time_batch)/value_t(this_batch),
        get_time_spent(),
        get_time_spent()/value_t(get_pts()),
        (value_t(time(NULL))-_time_batch-get_time_spent()+_last_time_spent)/value_t(this_batch));

    fprintf(output,"\n");
    fclose(output);

    _last_written=get_pts();
    _time_batch=value_t(time(NULL));
    _last_time_spent=get_time_spent();

}
