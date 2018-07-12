#include "exampleLikelihoods.h"
#include "dalex_driver.h"

int main(int iargc, char *argv[]){


    char out_name[letters];
    int i;
    for(i=0;argv[1][i]!=0 && i<letters-1;i++){
        out_name[i]=argv[1][i];
    }
    printf("outname %s\n",out_name);

    array_1d<double> min,max;
    for(i=0;i<12;i++){
        min.set(i,-40.0);
        max.set(i,40.0);
    }

    array_1d<int> seed_list;
    for(i=1;i<1000000;i+=1523){
        seed_list.add(i);
    }

    jellyBeanData *chisq;
    dalex_driver *dalex_driver_test;

    FILE *output_file;
    output_file=fopen(out_name,"w");
    fclose(output_file);

    int i_s;
    int overage=0;
    for(i_s=0;i_s<seed_list.get_dim();i_s++){
        chisq=new gaussianJellyBean12;
        dalex_driver_test=new dalex_driver;

        dalex_driver_test->set_seed(seed_list.get_data(i_s));
        dalex_driver_test->set_confidence_limit(0.95);
        dalex_driver_test->set_deltachi(21.0);
        dalex_driver_test->set_chisquared(chisq);
        dalex_driver_test->set_min(min);
        dalex_driver_test->set_max(max);
        dalex_driver_test->initialize(100);
        dalex_driver_test->mcmc_init();
        if(dalex_driver_test->get_chimin()>999.0){
            overage++;
        }

        output_file=fopen(out_name,"a");
        fprintf(output_file,"%d %d %e\n",
        seed_list.get_data(i_s),dalex_driver_test->get_called(),dalex_driver_test->get_chimin());
        fclose(output_file);

        printf("\nseed %d min %e overage %d out of %d %d\n\n",
        seed_list.get_data(i_s),dalex_driver_test->get_chimin(),overage,i_s+1,
        dalex_driver_test->get_called());
        delete dalex_driver_test;
        delete chisq;
    }


}
