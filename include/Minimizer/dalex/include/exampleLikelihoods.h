#ifndef EXAMPLE_LIKELIHOODS_H
#define EXAMPLE_LIKELIHOODS_H

#include "jellyBean.h"

class integrableJellyBean : public jellyBeanData{

    public:
        integrableJellyBean() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,20.0);
            _widths.set(0,1,2.0);
            _widths.set(0,2,10.0);
            _widths.set(0,3,10.0);

            _parabola_curvature=0.5;

            _bases.set(0,0,0.3);
            _bases.set(0,1,sqrt(1.0-0.3*0.3));
            _bases.set(0,2,0.0);
            _bases.set(0,3,0.0);

            _bases.set(1,0,-1.0*_bases.get_data(0,1));
            _bases.set(1,1,_bases.get_data(0,0));
            _bases.set(1,2,0.0);
            _bases.set(1,3,0.0);

            _bases.set(2,0,0.0);
            _bases.set(2,1,0.0);
            _bases.set(2,2,0.4);
            _bases.set(2,3,sqrt(1.0-0.4*0.4));

            _bases.set(3,0,0.0);
            _bases.set(3,1,0.0);
            _bases.set(3,2,_bases.get_data(2,3));
            _bases.set(3,3,-1.0*_bases.get_data(2,2));
        }
};

class integrableJellyBeanXX : public chisquared{

    public:
        integrableJellyBeanXX(int ii_dim): chisquared(ii_dim){

            int dxx=ii_dim-4;

            int i;
            for(i=0;i<dxx;i++){
                _dxx_center.set(i,_dice->doub()*(10.0)-5.0);
                _dxx_widths.set(i,_dice->doub()*0.9+0.1);
            }
            array_1d<double> trial;
            int j,k;
            double component;
            while(_dxx_bases.get_rows()<dxx){
                for(i=0;i<dxx;i++){
                    trial.set(i,normal_deviate(_dice,0.0,1.0));
                }
                for(j=0;j<_dxx_bases.get_rows();j++){
                    component=0.0;
                    for(k=0;k<dxx;k++){
                        component+=trial.get_data(k)*_dxx_bases.get_data(j,k);
                    }
                    for(k=0;k<dxx;k++){
                        trial.subtract_val(k,component*_dxx_bases.get_data(j,k));
                    }
                }
                component=trial.normalize();
                if(component>1.0e-10){
                   _dxx_bases.add_row(trial);
                }
            }

        }

        virtual double operator()(const array_1d<double> &pt){
            _called++;
            array_1d<double> four_d_pt;
            int i;
            for(i=0;i<4;i++){
                four_d_pt.set(i,pt.get_data(i));
            }
            double d4_val = _d4_chisq(four_d_pt);
            double dxx_val = 0.0;
            double component;
            int ix;
            int dxx=_dxx_center.get_dim();
            double ratio;
            for(ix=0;ix<dxx;ix++){
                component=0.0;
                for(i=0;i<dxx;i++){
                    component+=(pt.get_data(4+i)-_dxx_center.get_data(i))*_dxx_bases.get_data(ix,i);
                }
                ratio=fabs(component/_dxx_widths.get_data(ix));
                dxx_val+=ratio*log(ratio+1.0);
            }

            return d4_val + dxx_val;
        }

    private:
        integrableJellyBean _d4_chisq;
        array_1d<double> _dxx_center;
        array_2d<double> _dxx_bases;
        array_1d<double> _dxx_widths;

};

class gaussianJellyBean4 : public jellyBeanData{

    public:
        gaussianJellyBean4() :
        jellyBeanData(4,1,1.0,100,0.4){
            _widths.set(0,0,250.0);
            _widths.set(0,1,2.0);
            _widths.set(0,2,1.0);
            _widths.set(0,3,3.0);

            _parabola_curvature=8.0;
        }

};

class gaussianJellyBean12 : public jellyBeanData{

    public:
        gaussianJellyBean12() :
        jellyBeanData(12,1,1.0,100,0.4){
            Ran constructor_dice(44);
            int i;

            _widths.set(0,0,100.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=0.5;

            for(i=2;i<12;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }
};


class gaussianJellyBean16 : public jellyBeanData{

    public:
        gaussianJellyBean16() :
        jellyBeanData(16,1,1.0,100,0.4){
            Ran constructor_dice(44);
            int i;

            _widths.set(0,0,20.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=0.1;

            for(i=2;i<16;i++){
                _widths.set(0,i,3.5+3.0*constructor_dice.doub());
            }
        }
};


class gaussianJellyBean8 : public jellyBeanData{

    public:
        gaussianJellyBean8() :
        jellyBeanData(8,1,1.0,100,0.4){
            Ran constructor_dice(44);
            int i;

            _widths.set(0,0,100.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=0.5;

            for(i=2;i<8;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }
};

class nonGaussianLump8 : public gaussianJellyBean8{

    public:
        nonGaussianLump8(){}

    virtual void convert_params(const array_1d<double> &pt_in, array_1d<double> &out, int ic){
        array_1d<double> pt;
        pt.set_name("convert_params_pt");
        project_to_basis(pt_in,pt);

        int ix;
        for(ix=0;ix<_dim;ix++){
            if(ix%4==3){
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

        double xx;
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

};

class nonGaussianLump12 : public gaussianJellyBean12{

    public:
        nonGaussianLump12(){}

    virtual void convert_params(const array_1d<double> &pt_in, array_1d<double> &out, int ic){
        array_1d<double> pt;
        pt.set_name("convert_params_pt");
        project_to_basis(pt_in,pt);

        int ix;
        for(ix=0;ix<_dim;ix++){
            if(ix%4==3){
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

        double xx;
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

};



class gaussianJellyBean24 : public jellyBeanData{

    public:
        gaussianJellyBean24() :
        jellyBeanData(24,1,1.0,100,0.4){
            Ran constructor_dice(99);
            int i;

            _widths.set(0,0,250.0);
            _widths.set(0,1,2.0);
            _parabola_curvature=8.0;

            for(i=2;i<24;i++){
                _widths.set(0,i,0.5+3.0*constructor_dice.doub());
            }
        }

};


class nonGaussianLump24 : public gaussianJellyBean24{

    public:
        nonGaussianLump24(){}

    virtual void convert_params(const array_1d<double> &pt_in, array_1d<double> &out, int ic){
        array_1d<double> pt;
        pt.set_name("convert_params_pt");
        project_to_basis(pt_in,pt);

        int ix;
        for(ix=0;ix<_dim;ix++){
            if(ix%4==3){
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

        double xx;
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

};

#endif
