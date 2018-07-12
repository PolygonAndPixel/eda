#include "kd.h"

int get_dex(double value, double min, double dx){
    int ii=0;
    double test=min;
    while(test<value){
        ii++;
        test+=dx;
    }

    if(ii==0) return ii;

    if(test-value>0.5*dx){
        ii--;
    }
    return ii;

}

int main(int iargc, char *argv[]){

    int resolution=50;

    int i,j;
    int dim;
    char input_name[letters],output_root[letters];
    double chisquared_target;

    array_1d<int> ix, iy;
    ix.set_name("ix");
    iy.set_name("iy");


    for(i=1; i<iargc; i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'i':
                    i++;
                    transcribe(argv[i], input_name);
                break;
                case 'o':
                    i++;
                    transcribe(argv[i], output_root);
                break;
                case 'c':
                    i++;
                    chisquared_target=atof(argv[i]);
                break;
                case 'x':
                    i++;
                    ix.add(atoi(argv[i]));
                    i++;
                    iy.add(atoi(argv[i]));
                break;
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                break;
                case 'h':
                    printf("d -- dim\nc -- chisq\n");
                    printf("i -- input\no -- output\n");
                    printf("x -- dexes\n");
                    exit(1);
                break;
            }
        }
    }


    if(ix.get_dim()==0){
        for(i=0;i<dim;i++){
            for(j=i+1;j<dim;j++){
                ix.add(i);
                iy.add(j);
            }
        }
    }


    array_1d<double> vv;
    vv.set_name("vv");

    char word[letters];
    double mu;
    array_1d<double> min,max;
    array_1d<double> dx;
    min.set_name("xmin");
    max.set_name("xmax");
    dx.set_name("dx");

    FILE *input;

    array_2d<double> raw_points;
    raw_points.set_name("raw_points");
    raw_points.set_cols(dim);

    for(i=0;i<dim;i++){
        min.set(i,2.0*exception_value);
        max.set(i,-2.0*exception_value);
    }

    input=fopen(input_name,"r");
    for(i=0;i<dim+5;i++){
        fscanf(input,"%s",word);
    }
    while(fscanf(input,"%le",&mu)>0){
        vv.set(0,mu);
        for(i=1;i<dim;i++){
            fscanf(input, "%le",&mu);
            vv.set(i,mu);
        }

        fscanf(input,"%le",&mu);
        if(mu<=chisquared_target){
            raw_points.add_row(vv);
            for(i=0;i<dim;i++){
                if(vv.get_data(i)<min.get_data(i)){
                    min.set(i,vv.get_data(i));
                }
                if(vv.get_data(i)>max.get_data(i)){
                    max.set(i,vv.get_data(i));
                }
            }
        }

        for(i=0;i<3;i++)fscanf(input,"%d",&j);
    }
    fclose(input);

    printf("raw_points %d %e\n",raw_points.get_rows(),chisquared_target);

    for(i=0;i<dim;i++){
        dx.set(i,(max.get_data(i)-min.get_data(i))/double(resolution));
    }

    FILE *output;

    char output_name[2*letters];

    int ii;
    int xdex,ydex;
    int _ix, _iy;
    array_2d<double> points,boundary;
    array_2d<int> has_been_plotted;
    has_been_plotted.set_name("has_been_plotted");
    points.set_name("points");
    boundary.set_name("boundary");
    points.set_cols(2);
    boundary.set_cols(2);
    has_been_plotted.set_dim(resolution+1,resolution+1);

    for(ii=0;ii<ix.get_dim();ii++){
        points.reset_preserving_room();
        boundary.reset_preserving_room();

        xdex=ix.get_data(ii);
        ydex=iy.get_data(ii);

        has_been_plotted.zero();

        for(i=0;i<raw_points.get_rows();i++){
            _ix=get_dex(raw_points.get_data(i,xdex),min.get_data(xdex),dx.get_data(xdex));
            _iy=get_dex(raw_points.get_data(i,ydex),min.get_data(ydex),dx.get_data(ydex));

            if(has_been_plotted.get_data(_ix,_iy)==0){
                has_been_plotted.set(_ix,_iy,1);
                j=points.get_rows();
                points.set(j,0,min.get_data(xdex)+_ix*dx.get_data(xdex));
                points.set(j,1,min.get_data(ydex)+_iy*dx.get_data(ydex));
            }
        }

        convert_to_boundary(points,dx.get_data(xdex),dx.get_data(ydex),boundary);
        sprintf(output_name,"%s_cc%.2f_%d_%d.txt",output_root,chisquared_target,xdex,ydex);
        output=fopen(output_name,"w");
        for(i=0;i<boundary.get_rows();i++){
            fprintf(output,"%e %e\n",boundary.get_data(i,0),boundary.get_data(i,1));
        }
        fclose(output);

    }

}
