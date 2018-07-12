#include "aps_extractor.h"

int main(int iargc, char *argv[]){

int i,j;
array_1d<double> lnchi_hist,xmin,xmax,chi_hist;
array_1d<int> totalHist;
array_1d<int> totalDistLnChi,totalDistChi;

for(i=0;i<120;i++){
    chi_hist.set(i,i*1.0+0.5);
    lnchi_hist.set(i,-2.0+i*0.1);
    totalHist.set(i,0);

    totalDistChi.set(i,0);
    totalDistLnChi.set(i,0);
}

char inputName[letters],word[letters],outputRoot[letters];

array_1d<int> xdexes,ydexes;
xdexes.set_name("xdexes");
ydexes.set_name("ydexes");

int dim,ncenters,limit=-1;
double target_chi,delta_chi;

dim=22;

target_chi=-1.0;
delta_chi=-1.0;

for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'h':
                printf("d = dim\nc = chisq\nl = limit\n");
                printf("o = output\ni = input\nx = xdexes\n");
                exit(1);
            case 'd':
                i++;
                dim=atoi(argv[i]);
                break;
            case 't':
                i++;
                target_chi=atof(argv[i]);
                break;
            case 'c':
                i++;
                delta_chi=atof(argv[i]);
                break;
            case 'l':
                i++;
                limit=atoi(argv[i]);
                break;
            case 'o':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    outputRoot[j]=argv[i][j];
                }
                outputRoot[j]=0;
                break;
            case 'i':
                i++;
                for(j=0;j<letters-1 && argv[i][j]!=0;j++){
                    inputName[j]=argv[i][j];
                }
                inputName[j]=0;
                break;
            case 'x':
                i++;
                xdexes.add(atoi(argv[i]));
                i++;
                ydexes.add(atoi(argv[i]));
                break;
        }
    }
}

if(delta_chi<0.0 && target_chi<0.0){
    printf("WARNING must specify either target_chi or delta_chi\n");
    exit(1);
}

if(delta_chi>0.0 && target_chi>0.0){
    printf("WARNING cannot specify both target_chi and delta_chi\n");
    exit(1);
}

array_2d<double> data;
array_1d<double> chisq,vv,vvprojected,mu,sig;
array_1d<int> ling;

for(i=0;i<dim;i++){
    xmin.set(i,2.0*exception_value);
    xmax.set(i,-2.0*exception_value);
}

int hdex;
double nn,chival;
array_1d<int> *hptr;

data.set_name("main_data");
chisq.set_name("main_chisq");
ling.set_name("main_ling");
vv.set_name("main_vv");
mu.set_name("main_mu");
sig.set_name("main_sig");
vvprojected.set_name("main_vvprojected");

FILE *input,*output;

double maxerr=-1.0;

double chi_min=2.0*exception_value;
int ct=0;

input=fopen(inputName,"r");
if(input==NULL){
    printf("\nWARNING could not open %s\n\n",inputName);
    exit(1);
}
for(i=0;i<dim+5;i++){
    fscanf(input,"%s",word);
}
while(fscanf(input,"%le",&nn)>0 && (limit<0 || ct<limit)){
    ct++;
    vv.set(0,nn);
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&nn);
        vv.set(i,nn);
    }

    fscanf(input,"%le",&chival);
    if(chival<chi_min)chi_min=chival;


    fscanf(input,"%le",&nn);

    fscanf(input,"%le",&nn);

    fscanf(input,"%d",&j);

}
fclose(input);
printf("chimin %e\n",chi_min);

if(target_chi<0.0){
    target_chi=chi_min+delta_chi;
}

printf("target %e delta %e\n",target_chi,delta_chi);

input=fopen(inputName,"r");
for(i=0;i<dim+5;i++){
    fscanf(input,"%s",word);
}

ct=0;
printf("word %s dim %d\n",word,dim);
while(fscanf(input,"%le",&nn)>0 && (limit<0 || ct<limit)){
    ct++;
    vv.set(0,nn);
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&nn);
        vv.set(i,nn);
    }

    fscanf(input,"%le",&chival);
    chisq.add(chival);

    fscanf(input,"%le",&nn);
    mu.add(nn);
    fscanf(input,"%le",&nn);
    sig.add(nn);
    fscanf(input,"%d",&j);
    ling.add(j);

    hdex=get_dex(lnchi_hist,log(chival-chi_min)+1.0e-6);

    hptr=NULL;

    if(hptr!=NULL){
        for(i=hdex;i<hptr->get_dim();i++){
            hptr->add_val(i,1);
        }
    }

    for(i=hdex;i<totalHist.get_dim();i++){
        totalHist.add_val(i,1);
    }

    if(log(chival-chi_min)+1.0e-6<lnchi_hist.get_data(lnchi_hist.get_dim()-1)+0.1){
        totalDistLnChi.add_val(hdex,1);
    }
    if(chival-chi_min<chi_hist.get_data(chi_hist.get_dim()-1)+0.5){
        hdex=get_dex(chi_hist,chival-chi_min);
        totalDistChi.add_val(hdex,1);
    }

    if(chival<=target_chi){
        for(i=0;i<dim;i++){
            if(vv.get_data(i)<xmin.get_data(i)){
                xmin.set(i,vv.get_data(i));
            }
            if(vv.get_data(i)>xmax.get_data(i)){
                xmax.set(i,vv.get_data(i));
            }
        }
    }

}

fclose(input);

printf("read in everything\n");

aps_extractor apsExtractor;
apsExtractor.set_filename(inputName);
apsExtractor.set_target(target_chi);
if(limit>0)apsExtractor.set_cutoff(limit);

char outname[letters];
int ix,iy;

if(xdexes.get_dim()==0){
    for(ix=0;ix<dim;ix++){
        for(iy=ix+1;iy<dim;iy++){
            xdexes.add(ix);
            ydexes.add(iy);
        }
    }
}

for(i=0;i<xdexes.get_dim();i++){
    ix=xdexes.get_data(i);
    iy=ydexes.get_data(i);

    sprintf(outname,"%s_%d_%d_frequentist.sav",outputRoot,ix,iy);
    apsExtractor.write_good_points(outname,ix,iy);

    //sprintf(outname,"%s_%d_%d_bayesian.sav",outputRoot,ix.get_data(i),iy.get_data(i));
    //apsExtractor.draw_bayesian_bounds(outname,ix.get_data(i),iy.get_data(i),0.95);

}

sprintf(outname,"%s_good.sav",outputRoot);
apsExtractor.write_good_points(outname);


sprintf(outname,"%s_histograms.sav",outputRoot);
output=fopen(outname,"w");
fprintf(output,"#lnchi simplex compass ricochet total\n");
for(i=0;i<lnchi_hist.get_dim();i++){
    fprintf(output,"%le %d\n",
    lnchi_hist.get_data(i),
    totalHist.get_data(i));
}
fclose(output);

sprintf(outname,"%s_distributions.sav",outputRoot);
output=fopen(outname,"w");
fprintf(output,"#lnDeltachi dN/dlnchi chi dN/dchi\n");
for(i=0;i<lnchi_hist.get_dim();i++){
    fprintf(output,"%le %d %le %d\n",
    lnchi_hist.get_data(i),totalDistLnChi.get_data(i),
    chi_hist.get_data(i)+chi_min,totalDistChi.get_data(i));
}
fclose(output);

double volume=1.0;
for(i=0;i<dim;i++){
    volume*=xmax.get_data(i)-xmin.get_data(i);
    printf("%d %e %e\n",i,xmin.get_data(i),xmax.get_data(i));
}

printf("volume %e \n",volume);

}
