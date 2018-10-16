#include "dalex_driver.h"

dalex_driver::dalex_driver(){
    _ct_dalex=0;
    _last_did_min=0;
    sprintf(_outname,"output/carom_output.sav");
    sprintf(_timingname,"output/carom_timing.sav");
    _good_points.set_name("dale_driver_good_points");
    _good_points.set_room(100000);
}

dalex_driver::~dalex_driver(){
    _chifn.write_pts();
}

value_t dalex_driver::evaluate(array_1d<value_t> &pt, index_t *dex){
    value_t mu;
    _chifn.evaluate(pt,&mu,dex);

    if(mu<_chifn.target() && dex[0]>=0){
        if(_good_points.contains(dex[0])==0){
            _good_points.add(dex[0]);
        }
    }

    return mu;
}

void dalex_driver::assess_good_points(){
    assess_good_points(0,_chifn.get_pts());
}

void dalex_driver::assess_good_points(index_t i_min){
    assess_good_points(i_min,_chifn.get_pts());
}

void dalex_driver::assess_good_points(index_t i_min, index_t i_max){
    index_t i;

    for(i=0;i<_good_points.get_dim();i++){
        if(_chifn.get_fn(_good_points.get_data(i))>_chifn.target()){
            _good_points.remove(i);
            i--;
        }
    }

    for(i=i_min;i<i_max+1 && i<_chifn.get_pts();i++){
        if(_chifn.get_fn(i)<_chifn.target()){
            if(_good_points.contains(i)==0){
                _good_points.add(i);
            }
        }
    }
}

void dalex_driver::set_dof(index_t dd){
    _chifn.set_dof(dd);
}

void dalex_driver::set_confidence_limit(value_t cc){
    _chifn.set_confidence_limit(cc);
}

void dalex_driver::set_seed(index_t ii){
    _chifn.set_seed(ii);
}

void dalex_driver::set_min(array_1d<value_t> &vv){
    _chifn.set_min(vv);
}

void dalex_driver::set_max(array_1d<value_t> &vv){
    _chifn.set_max(vv);
}

void dalex_driver::set_characteristic_length(index_t dex, value_t vv){
    _chifn.set_characteristic_length(dex,vv);
}

void dalex_driver::set_chisquared(chisquared *xx){
    _chifn.set_chisquared(xx);
}

void dalex_driver::set_deltachi(value_t xx){
    _chifn.set_deltachi(xx);
}

void dalex_driver::set_target(value_t tt){
    _chifn.set_target(tt);
}

void dalex_driver::set_write_every(index_t ww){
    _chifn.set_write_every(ww);
}

index_t dalex_driver::get_dim(){
    return _chifn.get_dim();
}

index_t dalex_driver::get_called(){
    return _chifn.get_called();
}

value_t dalex_driver::get_chimin(){
    return _chifn.chimin();
}

void dalex_driver::set_outname(char *nn){
    index_t i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _outname[i]=nn[i];
    }
    _outname[i]=0;
}

void dalex_driver::set_timingname(char *nn){
    index_t i;
    for(i=0;i<letters-1 && nn[i]!=0;i++){
        _timingname[i]=nn[i];
    }
    _timingname[i]=0;
}

void dalex_driver::initialize(index_t npts){
    _chifn.initialize(npts);
    _cloud.build(&_chifn);
    assess_good_points(0);
}

void dalex_driver::mcmc_init(){
    dalex_initializer initializer;
    initializer.set_chifn(&_chifn);
    initializer.search();
}

void dalex_driver::search(index_t limit){
    index_t pt_start;

    value_t min0=_chifn.chimin();
    _chifn.set_outname(_outname);
    _chifn.set_timingname(_timingname);
    printf("before init min %e\n",_chifn.chimin());
    _chifn.set_search_type(_type_init);
    mcmc_init();
    printf("min now %e -> %e\n",min0,_chifn.chimin());
    printf("called %d\n",_chifn.get_pts());

    _cloud.set_limit(limit);

    while(_chifn.get_pts()<limit){

        pt_start=_chifn.get_pts();
        _cloud.search();
        _ct_dalex+=_chifn.get_pts()-pt_start;
    }

    index_t i;
    printf("minpt -- %e\n",_chifn.get_fn(_chifn.mindex()));
    for(i=0;i<_chifn.get_dim();i++){
        printf("    %.3e\n",_chifn.get_pt(_chifn.mindex(),i));
    }
    printf("\n\n");
    _chifn.write_pts();
}
