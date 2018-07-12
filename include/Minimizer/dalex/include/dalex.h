#ifndef DALEX_H
#define DALEX_H

#include "containers.h"
#include "goto_tools.h"
#include "chisq_wrapper.h"
#include "simplex.h"
#include "eigen_wrapper.h"
#include "explorers.h"
#include "ellipse.h"

class dalex{

    public:
        ~dalex(){
            FILE *ellipse_file;
            char out_name[500];
            sprintf(out_name,"ellipse_corners_%d.txt",_chifn->get_seed());
            ellipse_file=fopen(out_name, "w");
            int i,j,zone;
            double sign;
            double component;
            for(zone=0;zone<_exclusion_zones.ct();zone++){
                for(i=0;i<_chifn->get_dim();i++){
                    for(sign=-1.0;sign<1.5;sign+=2.0){
                        for(j=0;j<_chifn->get_dim();j++){
                            component=_exclusion_zones(zone)->center(j) +
                                      sign*_exclusion_zones(zone)->radii(i)*_exclusion_zones(zone)->bases(i,j);
                            fprintf(ellipse_file,"%e ",component);
                        }
                        fprintf(ellipse_file,"\n");
                    }
                }
            }
            fclose(ellipse_file);
        };

        dalex(){
            _chifn=NULL;
            _reset_threshold=0.5;
            _reset_chimin=2.0*exception_value;
            _strikes=0;
            _strikeouts=0;
            _has_struck=0;
            _limit=-1;
            _good_points.set_name("dalex_good_points");
            _good_points.set_room(100000);
            _tendril_path.set_name("dalex_tendril_path");
            _tendril_path.set_cols(2);
            _tendril_path.set_row_room(20000);
            _target_factor=1.0;
            _last_checked_good=0;
            _simplex_mindex=-1;

            _basis_chimin=2.0*exception_value;
            _basis_associates.set_name("dalex_basis_associates");
            _basis_mm.set_name("dalex_basis_mm");
            _basis_vv.set_name("dalex_basis_vv");
            _basis_norm.set_name("dalex_basis_norm");
            _basis_associate_norm.set_name("dalex_basis_associate_norm");
            _basis_vectors.set_name("dalex_basis_vectors");
            _basis_ddsq.set_name("dalex_basis_ddsq");

            _minimizers.set_name("dalex_minimizers");
            _tendril_init=0;
        };

        void build(chisq_wrapper*);

        void set_limit(int ii){
            _limit=ii;
        }

        void search();
        void simplex_search();
        void simplex_search(int);
        void simplex_search(array_1d<int>&);
        int simplex_boundary_search(const int, const int, ellipse_list&, int*);
        int _exploration_simplex(int,int,array_1d<int>&);
        void tendril_search();
        void init_fill();
        void find_tendril_candidates(double);
        void get_new_tendril(int*,int*);
        void min_explore(int, int);
        void initialize_min_exploration();

        int bisection(int, int, double, double);
        int bisection(int, const array_1d<double>&, double, double);
        int bisection(const array_1d<double>&, const array_1d<double>&, double, double);

        double get_basis(int i, int j){
            return _basis_vectors.get_data(i,j);
        }

        void set_target_factor(double tt){
            _target_factor=tt;
        }

        double get_norm(int);

    private:

       double distance(int i1, int i2){
           int i,j;
           double component,dd;
           dd=0.0;
           for(i=0;i<_chifn->get_dim();i++){
               component=0.0;
               for(j=0;j<_chifn->get_dim();j++){
                   component+=(_chifn->get_pt(i1,j)-_chifn->get_pt(i2,j))*_basis_vectors.get_data(i,j);
               }
               dd+=power(component/_basis_norm.get_data(i),2);
           }
           return sqrt(dd);
       }

      double cardinal_distance(int i1, int i2){
          int i;
          double dd=0.0;
          for(i=0;i<_chifn->get_dim();i++){
              dd+=power((_chifn->get_pt(i1,i)-_chifn->get_pt(i2,i))/
                         _chifn->get_characteristic_length(i),2);
          }
          return sqrt(dd);
      }

       void add_good_point(int ii){
           if(_chifn->get_fn(ii)<target() && _good_points.contains(ii)==0){
               _good_points.add(ii);
           }
       }

        void evaluate(const array_1d<double> &pt, double *mu_out, int *i_out){
            _chifn->evaluate(pt,mu_out,i_out);
            if(mu_out[0]<target() && _good_points.contains(i_out[0])==0){
                add_good_point(i_out[0]);
            }
        }

        void _update_good_points(){
            _update_good_points(_last_checked_good);
        }

        void _update_good_points(int i_start){
            safety_check("_update_good_points");

            int i;
            for(i=i_start;i<_chifn->get_pts();i++){
                if(_chifn->get_fn(i)<target() && _good_points.contains(i)==0){
                    add_good_point(i);
                }
            }
            _last_checked_good=_chifn->get_pts();
        }

        void assess_good_points(){
            int i;
            for(i=0;i<_good_points.get_dim();i++){
                if(_chifn->get_fn(_good_points.get_data(i))>target()){
                    _good_points.remove(i);
                    i--;
                }
            }
        }

        double chimin(){
            safety_check("chimin");
            return _chifn->chimin();
        }

        int mindex(){
            safety_check("mindex");
            return _chifn->mindex();
        }

        void safety_check(char *word){
            if(_chifn==NULL){
                printf("ERROR: dalex called %s but _chifn is NULL\n", word);
                exit(1);
            }

        }

        double target(){
            safety_check("target");
            return _target_factor*_chifn->target();
        }

        chisq_wrapper *_chifn;
        double _target_factor;
        double _reset_threshold,_reset_chimin;
        int _simplex_mindex;
        array_1d<int> _good_points;

        ////////code related to finding basis vectors
        array_1d<int> _basis_associates;
        array_1d<double> _basis_mm,_basis_bb,_basis_model,_basis_vv;
        array_1d<double> _basis_norm,_basis_associate_norm;
        array_2d<double> _basis_vectors,_basis_ddsq;
        double _basis_chimin;

        double basis_error(array_2d<double>&, array_1d<double>&);
        void find_trial_bases(int, array_1d<double>&, array_2d<double> &out_bases);
        void validate_bases(array_2d<double>&, char*);
        void guess_bases(array_2d<double>&);
        void find_covariance_matrix(int, array_2d<double>&);
        void find_bases();


        ///////code related to explorers
        explorers _min_explorers;
        int _last_checked_good;

        /////code related to minimizers
        array_1d<int> _minimizers;
        void refine_minimum();
        void iterate_on_minimum();

        //////code related to tendrils
        void get_negative_gradient(int, cost_fn&, ellipse&, array_1d<double>&);
        array_2d<int> _tendril_path;

        array_1d<int> _particles,_origins;
        array_1d<int> _particle_candidates,_origin_candidates;
        array_1d<int> _strikes_arr;

        void compass_search(ellipse&);

        int _limit;
        int _strikes,_strikeouts;
        int _has_struck;
        int _tendril_init;

        ellipse_list _exclusion_zones;
        ellipse_sampler _ellipse_sampler;
};

#endif
