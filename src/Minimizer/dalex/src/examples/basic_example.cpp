#include "dalex_driver.h"
#include "chisq.h"

class my_chisq_fn : public chisquared{

    public:

    // the constructor should call the constructor for
    // chisquared, passing in one int for the dimensionality
    // of the likelihood function
    my_chisq_fn() : chisquared(6){
       _called=0;
       _dim=6;
       _center.set(0,1.3);
       _center.set(1,6.7);
       _center.set(2,4.5);
       _center.set(3,-0.5);
       _center.set(4,-2.7);
       _center.set(5,11.4);

    }

    // you must define an operator that accepts an array_1d<double>
    // representing the point in parameter space and returns
    // a double representing the chi^2 value at that point
    virtual double operator()(const array_1d<double> &in){
        int i;
        double ans=0;
        for(i=0;i<_dim;i++){
            ans+=power((in.get_data(i)-_center.get_data(i))/2.0,2);
        }
        return ans;
    }

    array_1d<double> _center;
};

int main(){
    my_chisq_fn chifn;
    array_1d<double> min,max;
    int i;
    for(i=0;i<chifn.get_dim();i++){
        min.set(i,-20.0);
        max.set(i,20.0);
    }

    dalex_driver dalex_test;
    dalex_test.set_deltachi(12.03); // set delta chi^2 defining chi^2_lim
    dalex_test.set_seed(883); // seed the random number generator

    dalex_test.set_min(min); // set the minimum range in parameter space
    dalex_test.set_max(max); // set the maximum range in parameter space
    // the above limits just bound the initial random samples and the
    // seeds of the first optimization search; if Dalex finds interesting
    // regions of parameter space outside of these bounds, it will go there

    dalex_test.set_chisquared(&chifn); // pass in the chi^2 function
    dalex_test.set_timingname("output/workspace/basic_timing.txt");
    dalex_test.set_outname("output/workspace/basic_output.txt");
    dalex_test.initialize(12); // randomly sample 2D points
    dalex_test.set_write_every(10000); // how often to write output
    dalex_test.search(40000); // sample 40,000 points
}
