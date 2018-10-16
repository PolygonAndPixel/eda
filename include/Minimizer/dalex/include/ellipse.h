#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "containers.h"
#include "goto_tools.h"

class ellipse{

    public:

        ellipse(){
            _bases.set_name("ellipse_bases");
            _radii.set_name("ellipse_radii");
            _center.set_name("ellipse_center");
            _min.set_name("ellipse_min");
            _max.set_name("ellipse_max");
            _use_geo_center=0;
        }

        ~ellipse(){}
        void build(const array_1d<value_t>&, const array_2d<value_t>&);
        void build(const array_2d<value_t>&);

        void use_geo_center(){
            _use_geo_center=1;
        }

        void do_not_use_geo_center(){
            _use_geo_center=0;
        }
        index_t get_dim(){return _bases.get_rows();}
        index_t contains(const array_1d<value_t>&, const index_t);
        index_t contains(const array_1d<value_t>&);
        value_t bases(index_t i,index_t j){return _bases.get_data(i,j);}
        value_t center(index_t i){return _center.get_data(i);}
        value_t radii(index_t i){return _radii.get_data(i);}
        index_t dim(){return _center.get_dim();}
        void copy(ellipse&);

    private:
        array_2d<value_t> _bases;
        array_1d<value_t> _radii,_center,_min,_max;
        index_t _use_geo_center;

        void _set_radii(const array_2d<value_t>&);
        void _find_center(const array_2d<value_t>&);
};


class ellipse_list{
    public:
        ellipse_list(){
            _ellipse_list=NULL;
            _ct=0;
        }

        ~ellipse_list(){
            if(_ellipse_list!=NULL){
                delete [] _ellipse_list;
            }
        }

        void add(ellipse&);
        index_t ct(){return _ct;}
        void reset(){
            _ct=0;
            delete [] _ellipse_list;
            _ellipse_list=NULL;
        }

        ellipse* operator()(index_t ii){
            if(ii<0 || ii>=_ct){
                printf("No such ellipse: %d\n",ii);
                exit(1);
            }
            return &_ellipse_list[ii];
        }

    private:
        index_t _ct;
        ellipse *_ellipse_list;
};

#endif
