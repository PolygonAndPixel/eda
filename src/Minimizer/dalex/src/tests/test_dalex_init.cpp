#include "dalex_initializer.h"
#include "exampleLikelihoods.h"

int main(){

    int iterations=100;

    array_1d<double> min,max;
    int ii, dim, jj;
    chisq_wrapper *chifn_wrapped;
    integrableJellyBean *chifn_integrable;
    dalex_initializer *initializer;

    FILE *output;
    output=fopen("output/scratch/init_test_data.txt", "w");

    int seed,d_seed;
    d_seed=17;
    dim=4;
    fprintf(output,"dim %d\n",dim);
    for(jj=0;jj<dim;jj++){
        min.set(jj,-40.0);
        max.set(jj,40.0);
    }
    for(seed=99,ii=0;ii<iterations;ii++){
        seed+=d_seed;

        chifn_integrable=new integrableJellyBean();
        chifn_wrapped = new chisq_wrapper();
        chifn_wrapped->set_min(min);
        chifn_wrapped->set_max(max);
        chifn_wrapped->set_chisquared(chifn_integrable);
        chifn_wrapped->set_deltachi(12.0);
        chifn_wrapped->set_seed(seed);
        chifn_wrapped->initialize(dim*2);

        initializer = new dalex_initializer();
        initializer->set_chifn(chifn_wrapped);
        initializer->search();

        fprintf(output,"%d %e %d\n",
        chifn_wrapped->get_called(),chifn_wrapped->chimin(),seed);

        delete chifn_integrable;
        delete chifn_wrapped;
        delete initializer;
    }
    fprintf(output,"\n");

    gaussianJellyBean12 *d12_chisq;

    dim=12;
    fprintf(output,"dim %d\n",dim);

    for(jj=0;jj<dim;jj++){
        min.set(jj,-40.0);
        max.set(jj,40.0);
    }
    for(seed=99,ii=0;ii<iterations;ii++){
        seed+=d_seed;

        d12_chisq=new gaussianJellyBean12();
        chifn_wrapped = new chisq_wrapper();
        chifn_wrapped->set_min(min);
        chifn_wrapped->set_max(max);
        chifn_wrapped->set_chisquared(d12_chisq);
        chifn_wrapped->set_deltachi(12.0);
        chifn_wrapped->set_seed(seed);
        chifn_wrapped->initialize(dim*2);

        initializer = new dalex_initializer();
        initializer->set_chifn(chifn_wrapped);
        initializer->search();

        fprintf(output,"%d %e %d\n",
        chifn_wrapped->get_called(),chifn_wrapped->chimin(),seed);

        delete d12_chisq;
        delete chifn_wrapped;
        delete initializer;
    }

    /*fprintf(output,"\n");

    gaussianJellyBean24 *d24_chisq;

    dim=24;
    fprintf(output,"dim %d\n",dim);

    for(jj=0;jj<dim;jj++){
        min.set(jj,-40.0);
        max.set(jj,40.0);
    }
    for(seed=99,ii=0;ii<iterations;ii++){
        seed+=d_seed;

        d24_chisq=new gaussianJellyBean24();
        chifn_wrapped = new chisq_wrapper();
        chifn_wrapped->set_min(min);
        chifn_wrapped->set_max(max);
        chifn_wrapped->set_chisquared(d24_chisq);
        chifn_wrapped->set_deltachi(12.0);
        chifn_wrapped->set_seed(seed);
        chifn_wrapped->initialize(dim*2);

        initializer = new dalex_initializer();
        initializer->set_chifn(chifn_wrapped);
        initializer->search();

        fprintf(output,"%d %e %d\n",
        chifn_wrapped->get_called(),chifn_wrapped->chimin(),seed);

        delete d24_chisq;
        delete chifn_wrapped;
        delete initializer;
    }*/

    fclose(output);

}
