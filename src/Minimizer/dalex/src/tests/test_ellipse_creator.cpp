#include "containers.h"
#include "goto_tools.h"
#include "ellipse.h"

int main(int iargc, char **argv){

    int dim=12;
    int n_pts=1000;
    int seed1=99;
    int seed2=342;

    int i,j,k;

    for(i=1;i<iargc;i++){
        if(argv[i][0]=='-'){
            switch(argv[i][1]){
                case 'd':
                    i++;
                    dim=atoi(argv[i]);
                    break;
            }
        }
    }

    array_1d<double> center;
    array_2d<double> bases;
    array_2d<double> pts;
    Ran dice(seed1);


    array_1d<double> dir;
    double component,norm;
    while(bases.get_rows()!=dim){
        for(i=0;i<dim;i++){
            dir.set(i,normal_deviate(&dice, 0.0, 1.0));
        }
        for(i=0;i<bases.get_rows();i++){
            component=0.0;
            for(j=0;j<dim;j++){
                component+=dir.get_data(j)*bases.get_data(i,j);
            }
            for(j=0;j<dim;j++){
                dir.subtract_val(j,component*bases.get_data(i,j));
            }
        }
        norm=dir.normalize();
        if(norm>1.0e-20){
            bases.add_row(dir);
        }
    }

    array_1d<double> radii;
    for(i=0;i<dim;i++){
        center.set(i,normal_deviate(&dice, 0.0, 10.0));
        radii.set(i,0.1+dice.doub()*2.0);
    }

    ellipse_sampler sampler;
    sampler.initialize(dim, seed2);
    array_1d<double> base_pt,true_pt;
    for(i=0;i<n_pts;i++){
        sampler.get_pt(base_pt);
        for(j=0;j<dim;j++){
            true_pt.set(j,center.get_data(j));
        }
        for(j=0;j<dim;j++){
            for(k=0;k<dim;k++){
                true_pt.add_val(k,base_pt.get_data(j)*radii.get_data(j)*bases.get_data(j,k));
            }
        }
        pts.add_row(true_pt);
    }

    ellipse ell;
    ell.build(pts);

    double dd=0.0;
    for(i=0;i<dim;i++){
        dd+=power(ell.center(i)-center.get_data(i),2);
    }
    printf("center distance %e\n",sqrt(dd));

    for(i=0;i<n_pts;i++){
        if(ell.contains(pts(i))==0){
            printf("WARNING does not contain a pt\n");
            exit(1);
        }
    }


    double dot,best_dot,local_best_dot;
    int best_dex;
    array_1d<int> chosen,chosen_control;
    array_1d<int> best_pair, local_best_pair;
    int pairs=0;
    while(chosen_control.get_dim()!=dim){
        best_dot=-1.0e30;
        for(i=0;i<dim;i++){
            if(chosen_control.contains(i)==0){
                local_best_dot=-1.0e30;
                for(j=0;j<dim;j++){
                    if(chosen.contains(j)==0){
                        dot=0.0;
                        for(k=0;k<dim;k++){
                            dot+=bases.get_data(i,k)*ell.bases(j,k);
                        }
                        dot=fabs(dot);
                        if(dot>local_best_dot){
                            local_best_dot=dot;
                            local_best_pair.set(0,i);
                            local_best_pair.set(1,j);
                        }
                    }
                }
                if(local_best_dot>best_dot){
                    best_dot=local_best_dot;
                    best_pair.set(0,local_best_pair.get_data(0));
                    best_pair.set(1,local_best_pair.get_data(1));
                }
            }
        }
        chosen_control.add(best_pair.get_data(0));
        chosen.add(best_pair.get_data(1));
        dot=0.0;
        for(i=0;i<dim;i++){
            dot+=bases.get_data(best_pair.get_data(0),i)*ell.bases(best_pair.get_data(1),i);
        }
        printf("%d %d -- %e -- %e %e\n",
               best_pair.get_data(0),
               best_pair.get_data(1),
               dot,
               radii.get_data(best_pair.get_data(0)),
               ell.radii(best_pair.get_data(1)));
    }

    printf("\n");
    if(dim<=6){
        printf("brute matching\n");
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
                dot=0.0;
                for(k=0;k<dim;k++){
                    dot+=bases.get_data(i,k)*ell.bases(j,k);
                }
                printf("%d %d -- %e -- %e %e\n",
                i,j,dot,radii.get_data(i),ell.radii(j));
            }
            printf("\n");
        }

        printf("\nbases\n");
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
                printf("%e ",bases.get_data(i,j));
            }
            printf("\n");
        }
    }

    FILE *output;
    output=fopen("test_ellipse_pts.txt","w");
    for(i=0;i<n_pts;i++){
        for(j=0;j<dim;j++){
            fprintf(output,"%e ",pts.get_data(i,j));
        }
        fprintf(output,"\n");
    }
    fclose(output);

    ////////test ellipse list
    array_2d<double> bases1,bases2,bases3;
    array_2d<double> pts1,pts2,pts3;
    array_1d<double> radii1,radii2,radii3;
    array_1d<double> center1,center2,center3;
    ellipse ellipse1,ellipse2,ellipse3;

    array_2d<double> *bs_ptr;
    array_2d<double> *pt_ptr;
    array_1d<double> *rad_ptr;
    array_1d<double> *c_ptr;
    ellipse *ell_ptr;

    ellipse_list ell_list;

    int ip;
    for(ip=0;ip<3;ip++){
        if(ip==0){
            bs_ptr=&bases1;
            pt_ptr=&pts1;
            rad_ptr=&radii1;
            c_ptr=&center1;
            ell_ptr=&ellipse1;
        }
        else if(ip==1){
            bs_ptr=&bases2;
            pt_ptr=&pts2;
            rad_ptr=&radii2;
            c_ptr=&center2;
            ell_ptr=&ellipse2;
        }
        else if(ip==2){
            bs_ptr=&bases3;
            pt_ptr=&pts3;
            rad_ptr=&radii3;
            c_ptr=&center3;
            ell_ptr=&ellipse3;
        }

        while(bs_ptr->get_rows()!=dim){
            for(i=0;i<dim;i++){
                dir.set(i,normal_deviate(&dice,0.0,1.0));
            }
            for(j=0;j<bs_ptr->get_rows();j++){
                component=0.0;
                for(k=0;k<dim;k++){
                    component+=dir.get_data(k)*bs_ptr->get_data(j,k);
                }
                for(k=0;k<dim;k++){
                    dir.subtract_val(k,component*bs_ptr->get_data(j,k));
                }
            }
            component=dir.normalize();
            if(component>1.0e-12){
                bs_ptr->add_row(dir);
            }
        }

        for(i=0;i<dim;i++){
            rad_ptr->set(i,dice.doub()*2.1);
            c_ptr->set(i,normal_deviate(&dice,0.0,10.0));
        }

        for(i=0;i<n_pts;i++){
            sampler.get_pt(base_pt);
            for(j=0;j<dim;j++){
                true_pt.set(j,c_ptr->get_data(j));
            }
            for(j=0;j<dim;j++){
                for(k=0;k<dim;k++){
                    true_pt.add_val(k,base_pt.get_data(j)*rad_ptr->get_data(j)*bs_ptr->get_data(j,k));
                }
            }
            pt_ptr->add_row(true_pt);
        }
        ell_ptr->build(pt_ptr[0]);
        ell_list.add(ell_ptr[0]);
    }

    double diff;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            diff=fabs(ell_list(0)->bases(i,j)-ellipse1.bases(i,j));
            if(diff>1.0e-10){
                printf("WARNING bases %d %d %e\n",i,j,diff);
                exit(1);
            }
        }
        diff=fabs(ell_list(0)->center(i)-ellipse1.center(i));
        if(diff>1.0e-10){
            printf("WARNING center %d %e\n",i,diff);
            exit(1);
        }
        diff=fabs(ell_list(0)->radii(i)-ellipse1.radii(i));
        if(diff>1.0e-10){
            printf("WARNING radii %d %e\n",i,diff);
            exit(1);
        }
    }

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            diff=fabs(ell_list(1)->bases(i,j)-ellipse2.bases(i,j));
            if(diff>1.0e-10){
                printf("WARNING bases %d %d %e\n",i,j,diff);
                exit(1);
            }
        }
        diff=fabs(ell_list(1)->center(i)-ellipse2.center(i));
        if(diff>1.0e-10){
            printf("WARNING center %d %e\n",i,diff);
            exit(1);
        }
        diff=fabs(ell_list(1)->radii(i)-ellipse2.radii(i));
        if(diff>1.0e-10){
            printf("WARNING radii %d %e\n",i,diff);
            exit(1);
        }
    }

    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            diff=fabs(ell_list(2)->bases(i,j)-ellipse3.bases(i,j));
            if(diff>1.0e-10){
                printf("WARNING bases %d %d %e\n",i,j,diff);
                exit(1);
            }
        }
        diff=fabs(ell_list(2)->center(i)-ellipse3.center(i));
        if(diff>1.0e-10){
            printf("WARNING center %d %e\n",i,diff);
            exit(1);
        }
        diff=fabs(ell_list(2)->radii(i)-ellipse3.radii(i));
        if(diff>1.0e-10){
            printf("WARNING radii %d %e\n",i,diff);
            exit(1);
        }
    }
}
