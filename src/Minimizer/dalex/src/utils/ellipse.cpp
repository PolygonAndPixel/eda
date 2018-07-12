#include "ellipse.h"

int ellipse::contains(const array_1d<double> &pt){
    return contains(pt, 0);
}

int ellipse::contains(const array_1d<double> &pt, const int use_extreme){
    int i,j;
    double sum=0.0;
    double component,threshold;
    for(i=0;i<_bases.get_rows();i++){
        component=0.0;
        if(i<_max.get_dim() && i<_min.get_dim()){
            threshold=0.025*(_max.get_data(i)-_min.get_data(i));
        }
        for(j=0;j<_bases.get_rows();j++){
            component+=(pt.get_data(j)-_center.get_data(j))*_bases.get_data(i,j);
        }
        if(use_extreme==1 && i<_max.get_dim() && component-_max.get_data(i)>threshold){
            return 0;
        }
        if(use_extreme==1 && i<_min.get_dim() && _min.get_data(i)-component>threshold){
            return 0;
        }
        sum+=power(component/_radii.get_data(i),2);
        if(sum>1.0){
            return 0;
        }
    }
    return 1;
}

void ellipse::build(const array_2d<double> &pts_in){

    array_1d<double> dummy_center;
    build(dummy_center, pts_in);

}

void ellipse::build(const array_1d<double> &center_in,
                    const array_2d<double> &pts_in){

    Ran dice(99);

    _min.reset_preserving_room();
    _max.reset_preserving_room();

    int n_pts = pts_in.get_rows();
    int dim = pts_in.get_cols();

    if(n_pts<dim){
        printf("WARNING cannot construct ellipse from %d pts in %d dimensions\n",
        n_pts,dim);
        exit(1);
    }

    _bases.reset_preserving_room();
    _radii.reset_preserving_room();

    int i;
    if(center_in.get_dim()!=pts_in.get_cols()){
        _find_center(pts_in);
    }
    else{
        for(i=0;i<center_in.get_dim();i++){
            _center.set(i,center_in.get_data(i));
        }
    }

    array_1d<double> dir,dir_max;
    dir.set_name("ellipse_build_dir");
    dir_max.set_name("ellipse_build_dir_max");

    int j,k;
    double component,norm,norm_max;

    while(_bases.get_rows()!=dim){
        for(i=0;i<n_pts;i++){
            for(j=0;j<dim;j++){
                dir.set(j,pts_in.get_data(i,j)-_center.get_data(j));
            }
            for(j=0;j<_bases.get_rows();j++){
                component=0.0;
                for(k=0;k<dim;k++){
                    component+=dir.get_data(k)*_bases.get_data(j,k);
                }
                for(k=0;k<dim;k++){
                    dir.subtract_val(k,component*_bases.get_data(j,k));
                }
            }
            norm=dir.normalize();
            if(i==0 || norm>norm_max){
                norm_max=norm;
                for(j=0;j<dim;j++){
                    dir_max.set(j,dir.get_data(j));
                }
            }
        }
        if(norm_max<1.0e-20){
            printf("WARNING ellipse basis had norm %e; setting randomly\n",
            norm_max);
            norm=-1.0;
            while(norm<1.0e-20){
                for(i=0;i<dim;i++){
                    dir.set(i,normal_deviate(&dice,0.0,1.0));
                }
                for(j=0;j<_bases.get_rows();j++){
                    component=0.0;
                    for(k=0;k<dim;k++){
                        component+=dir.get_data(k)*_bases.get_data(j,k);
                    }
                    for(k=0;k<dim;k++){
                        dir.subtract_val(k,component*_bases.get_data(j,k));
                    }
                }
                norm=dir.normalize();
                if(norm>1.0e-20){
                    for(k=0;k<dim;k++){
                        dir_max.set(k,dir.get_data(k));
                    }
                    norm_max=1.0e-10;
                }
            }
        }
        _bases.add_row(dir_max);
        _radii.add(norm_max);
    }

    int is_valid=0;

    while(is_valid==0){
        is_valid=1;
        for(i=0;i<n_pts;i++){
            if(contains(pts_in(i))==0){
                is_valid=0;
                break;
            }
        }

        if(is_valid==0){
            _set_radii(pts_in);
        }
    }

    //set extremities
    for(i=0;i<n_pts;i++){
        for(j=0;j<dim;j++){
            dir.set(j,pts_in.get_data(i,j)-_center.get_data(j));
        }

        for(j=0;j<dim;j++){
            component=0.0;
            for(k=0;k<dim;k++){
                component+=dir.get_data(k)*_bases.get_data(j,k);
            }
            if(j>=_min.get_dim() || component<_min.get_data(j)){
                _min.set(j,component);
            }
            if(j>=_max.get_dim() || component>_max.get_data(j)){
                _max.set(j,component);
            }
        }
    }

}


void ellipse::_find_center(const array_2d<double> &pts_in){
    int n_pts=pts_in.get_rows();
    int dim=pts_in.get_cols();

    array_1d<double> min,max;
    min.set_name("ellipse_build_min");
    max.set_name("ellipse_build_max");
    int i,j,k;
    _center.set_dim(dim);
    _center.zero();
    for(i=0;i<n_pts;i++){
        for(j=0;j<dim;j++){
            if(j>=min.get_dim() || pts_in.get_data(i,j)<min.get_data(j)){
                min.set(j,pts_in.get_data(i,j));
            }
            if(j>=max.get_dim() || pts_in.get_data(i,j)>max.get_data(j)){
                max.set(j,pts_in.get_data(i,j));
            }
        }
    }

    array_1d<double> geo_center;
    geo_center.set_name("ellipse_geo_center");
    for(i=0;i<dim;i++){
        geo_center.set(i,0.5*(min.get_data(i)+max.get_data(i)));
    }

    if(_use_geo_center==1){
        for(i=0;i<dim;i++){
            _center.set(i,geo_center.get_data(i));
        }
        return;
    }

    double dd;
    double ddmin;
    for(i=0;i<pts_in.get_rows();i++){
        dd=0.0;
        for(j=0;j<dim;j++){
           dd+=power((pts_in.get_data(i,j)-geo_center.get_data(j))/
                     (max.get_data(j)-min.get_data(j)),2);
        }
        if(i==0 || dd<ddmin){
            ddmin=dd;
            for(j=0;j<dim;j++){
                _center.set(j,pts_in.get_data(i,j));
            }
        }
    }
}

void ellipse::_set_radii(const array_2d<double> &pts_in){

    array_1d<int> bad_dexes;
    int i;
    for(i=0;i<pts_in.get_rows();i++){
        if(contains(pts_in(i))==0){
            bad_dexes.add(i);
        }
    }

    int ipt;
    array_1d<double> norm,norm_sorted;
    array_1d<int> norm_dex,norm_ct;
    for(i=0;i<pts_in.get_cols();i++){
        norm_ct.set(i,0);
    }

    int j,k;
    double component;
    for(i=0;i<bad_dexes.get_dim();i++){
        ipt=bad_dexes.get_data(i);
        for(j=0;j<pts_in.get_cols();j++){
            component=0.0;
            for(k=0;k<pts_in.get_cols();k++){
                component+=(pts_in.get_data(ipt,k)-_center.get_data(k))*_bases.get_data(j,k);
            }
            norm.set(j,fabs(component)/_radii.get_data(j));
            norm_dex.set(j,j);
        }
        sort(norm,norm_sorted,norm_dex);
        j=norm_dex.get_data(pts_in.get_cols()-1);
        norm_ct.add_val(j,1);
    }

    norm_dex.reset_preserving_room();
    for(i=0;i<pts_in.get_cols();i++){
        norm_dex.set(i,i);
    }
    array_1d<int> norm_ct_sorted;
    sort(norm_ct,norm_ct_sorted,norm_dex);
    i=norm_dex.get_data(pts_in.get_cols()-1);

    _radii.multiply_val(i,1.1);

}

void ellipse::copy(ellipse &other){
    _bases.reset_preserving_room();
    _center.reset_preserving_room();
    _radii.reset_preserving_room();
    _min.reset_preserving_room();
    _max.reset_preserving_room();
    int i;
    for(i=0;i<other.dim();i++){
        _center.set(i,other.center(i));
        _radii.set(i,other.radii(i));
    }
    _bases.set_dim(other.dim(),other.dim());
    int j;
    for(i=0;i<other.dim();i++){
        for(j=0;j<other.dim();j++){
            _bases.set(i,j,other.bases(i,j));
        }
    }

    for(i=0;i<other._min.get_dim();i++){
        _min.set(i,other._min.get_data(i));
    }
    for(i=0;i<other._max.get_dim();i++){
        _max.set(i,other._max.get_data(i));
    }
}

void ellipse_list::add(ellipse &ell_in){
    ellipse *buffer;
    int i;
    if(_ct>0){
        buffer=new ellipse[_ct];
        for(i=0;i<_ct;i++){
            buffer[i].copy(_ellipse_list[i]);
        }
        delete [] _ellipse_list;
        _ellipse_list=new ellipse[_ct+1];
        for(i=0;i<_ct;i++){
            _ellipse_list[i].copy(buffer[i]);
        }
        delete [] buffer;
    }
    else{
        _ellipse_list=new ellipse[1];
    }
    _ellipse_list[_ct].copy(ell_in);
    _ct++;
}
