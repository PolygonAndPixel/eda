#include "containers.h"
#include "wrappers.h"

#ifndef GOTO_H
#define GOTO_H
#define pi 3.141592654

#define exception_value 1.0e30

void kill(char*);

value_t raiseup(value_t,value_t);

inline value_t power(value_t arg,index_t raised){

  //return arg raised to the integer power 'raised'

  index_t n;
  value_t ans;

  if(raised==0)return 1.0;
  else{ans=1.0;
  for(n=0;n<raised;n++){
    ans=ans*arg;
  }

  return ans;
  }

}

void transcribe(char*,char*);

struct Ran{

//this structure will be based on the Xorshift random number generator
// discovered by George Marsaglia and published in
//Journal of Statistical Software, volume 8, no. 14 pp 1-6

//parameters are drawn from the table on page 347 of
//Numerical Recipes (3rd edition)
//William H. press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
//Cambridge University Press, 2007

unsigned long long x;
Ran(unsigned long long seed){

x=seed^88172645463325252LL;
x^=(x<<21);
x^=(x>>35);
x^=(x<<4);
//printf("starting rand with %ld from seed %d\n",x,seed);
}

void thework(){
  x^=(x<<21);
  x^=(x>>35);
  x^=(x<<4);
}

value_t doub(){
  thework();
  return x*5.42101086242752217e-20;
}

index_t int32(){
  thework();
  index_t ans=index_t(x);
  if(ans<0)ans=-1*ans;
  return ans;
}

};

void polint(value_t*,value_t*,index_t,value_t,value_t*,value_t*);

value_t interpolate(value_t*,value_t*,value_t,index_t);

void sort(value_t*,index_t*,index_t);

void check_sort(value_t*,index_t*,index_t);

void sort_and_check(value_t*,value_t*,index_t*,index_t);

value_t normal_deviate(Ran*,value_t,value_t);

void naive_gaussian_solver(array_1d<value_t>&,array_1d<value_t>&,
array_1d<value_t>&,index_t);

value_t compare_arr(array_1d<value_t>&,array_1d<value_t>&);

index_t compare_int_arr(array_1d<index_t>&, array_1d<index_t>&);

value_t bisection(function_wrapper&,array_1d<value_t>&,array_1d<value_t>&,value_t,value_t,array_1d<value_t>&);

value_t integrate_cos_n(value_t, value_t, index_t);

inline void expand_grid(index_t ii ,array_1d<index_t> &grid_ct, array_1d<index_t> &out){
    index_t ix,iy,denom,subtract;

    for(ix=0;ix<grid_ct.get_dim();ix++){
        denom=1;
        for(iy=ix+1;iy<grid_ct.get_dim();iy++){
            denom*=grid_ct.get_data(iy);
        }
        iy=ii/denom;
        out.set(ix,iy);
        subtract=ii/denom;
        ii-=subtract*denom;
    }

}


inline void expand_grid(long index_t ii ,array_1d<index_t> &grid_ct, array_1d<index_t> &out){
    index_t ix,iy;

    long index_t denom,subtract,quotient;

    for(ix=0;ix<grid_ct.get_dim();ix++){
        denom=1;
        for(iy=ix+1;iy<grid_ct.get_dim();iy++){
            denom*=grid_ct.get_data(iy);
        }
        quotient=ii/denom;
        iy=index_t(quotient);
        out.set(ix,iy);
        subtract=ii/denom;
        ii-=subtract*denom;
    }

}

struct chisquared_distribution{

    value_t _fix;

    value_t _pdf_fn(value_t x, value_t dof,value_t *lnpdf){
        value_t logans;
        logans=-0.5*x+(0.5*dof-1.0)*log(x)-_fix;

        lnpdf[0]=logans;
        return exp(logans);
    }

    value_t maximum_likelihood_chisquared_value(value_t dof){
        _fix=0.0;
        value_t x,lnpdf,pdf,maxpdf,x_max;
        maxpdf=-1.0;
        for(x=0.0;x<dof+1000.0;x+=1.0){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(lnpdf>maxpdf){
                maxpdf=lnpdf;
                x_max=x;
            }
        }

        return x_max;
    }

    value_t confidence_limit(value_t dof, value_t pct){
        _fix=0.0;

        value_t x,pdf,pdfold,xold,total,lim,ans,lnpdf;
        value_t maxf,maxx,lmax=-1.0e10;
        value_t start,step,stop,target;
        value_t x1,x2;

        maxf=-1.0e20;

        target=-1.0;

        start=1.0e-10;
        stop=dof+1000.0;
        step=0.01;

        xold=0.0;
        pdfold=0.0;
        total=0.0;

        _fix=0.0;
        for(x=start;x<stop;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(lnpdf>lmax)lmax=lnpdf;
        }
        _fix=lmax;

        lmax=-1.0e10;

        for(x=start;x<stop;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            if(lnpdf>lmax){
                maxf=pdf;
                maxx=x;
                lmax=lnpdf;
            }
            total+=0.5*(pdf+pdfold)*(x-xold);

            xold=x;
            pdfold=pdf;
        }


        lim=pct*total;
        ans=0.0;
        xold=start;

        pdfold=0.0;
        for(x=start;ans<lim;x+=step){
            pdf=_pdf_fn(x,dof,&lnpdf);
            ans+=0.5*(pdfold+pdf)*(x-xold);
            xold=x;
            pdfold=pdf;

        }

        return xold;

    }
};


class ellipse_sampler{

    public:
        ellipse_sampler();
        ~ellipse_sampler(){
            if(_dice!=NULL){
                delete _dice;
            }
        }

        index_t is_initialized(){return _initialized;}

        void initialize(index_t, index_t);
        void get_pt(array_1d<value_t>&);

    private:
        index_t _initialized;
        index_t _dim;
        index_t _steps;
        value_t _dx;
        array_1d<index_t> _dim_record;
        Ran *_dice;
        array_2d<value_t> _cos_n_grid;
};

#endif
