#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "jellyBean.h"
#include "exampleLikelihoods.h"
#include "dalex_driver.h"

int main(int iargc, char *argv[]){

int seed=-1;
double delta_chi=21.03;
int nsamples=-1;
char output_name[letters];
char timing_name[letters];

int i,j;

timing_name[0]=0;
output_name[0]=0;
for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 's':
                i++;
                seed=atoi(argv[i]);
                break;
            case 'n':
                i++;
                nsamples=atoi(argv[i]);
                break;
            case 'o':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    output_name[j]=argv[i][j];
                }
                output_name[j]=0;
                break;
            case 't':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    timing_name[j]=argv[i][j];
                }
                timing_name[j]=0;
                break;
        }
    }
}

if(iargc==1 || seed<0 || nsamples<0 ||
   timing_name[0]==0 || output_name[0]==0){
    printf("call signature\n");
    printf("./bin/curved_d_example -s seed -n nsamples -o output_name -t timing_name\n");
    exit(1);
}

gaussianJellyBean12 chifn;

dalex_driver dalex_test;

// define the delta_chi^2 in chi^2_lim = chi^2_min + delta_chi^2
dalex_test.set_deltachi(delta_chi);

// seed the random number generator
dalex_test.set_seed(seed);

// pass in the chi^2 function
dalex_test.set_chisquared(&chifn);

// Set the original parameter space bounds over which to search.
// Dalex will stray outside of these bounds, if needs be.
// Mostly, these bounds will help Dalex place the original,
// random chi^2 samples.
array_1d<double> max,min;
for(i=0;i<chifn.get_dim();i++){
    min.set(i,-40.0);
    max.set(i,40.0);
}

dalex_test.set_min(min);
dalex_test.set_max(max);

// specify the name of the file where timing information will be written
dalex_test.set_timingname(timing_name);

// specify the name of the file where the output data will be written
dalex_test.set_outname(output_name);

// randomly sample at least 2D points in parameter space
// (this helps establish the KD Tree that Dalex uses to
// keep track of its history).
dalex_test.initialize(24);

// actually perform the search.
dalex_test.search(nsamples);

}
