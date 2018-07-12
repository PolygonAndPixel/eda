#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "goto_tools.h"


void transcribe(char *w1, char *w2){
    int i;
    for(i=0;i<letters-1 && w1[i]!=0;i++){
        w2[i]=w1[i];
    }
    w2[i]=0;
}

double normal_deviate(Ran *chaos, double mu, double sig){

 int i;
 double x1,x2,y1,y2;

 x1=chaos->doub();
 x2=chaos->doub();

 y1=sqrt(-2.0*log(x1))*cos(2.0*pi*x2);
 y2=sqrt(-2.0*log(x1))*sin(2.0*pi*x2);


 return mu+y1*sig;

}


void kill(char *words){
 double junk;
 FILE *output;

 //write the character string words[] to the terminal
 //then hang, waiting for an input double

 printf("%s\n",words);
 scanf("%lf",&junk);
}

/*just going to lift this out of Numerical Recipes p 109*/
void polint(double *xa, double *ya, int n, double x, double *y, double *dy){
/*recall: because this is from NR, it has arrays start at element unity*/
	int i,m,ns=1,isfd;
	double den,dif,dift,ho,hp,w;
	double c[n+1],d[n+1];
	char scream[100];
	dif=fabs(x-xa[1]);
	for(i=1;i<=n;i++){
		if((dift=fabs(x-xa[i]))<dif){
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for(m=1;m<n;m++){
		for(i=1;i<=n-m;i++){
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if((den=ho-hp)==0.0){printf("Error in routine polint ");
				for(isfd=1;isfd<=n;isfd++)printf(" (%e, %e) ",xa[isfd],ya[isfd]);
				printf(" want %e \n",x);
				sprintf(scream,"stop");
				kill(scream);
			
			
			}
			/*This error can occur only if two input xas re (to within roundoff) identical*/
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y+=(*dy=(2*ns<(n-m)?c[ns+1]:d[ns--]));
	}
}


double interpolate(double *x, double *y, double target, int el){

  //this is a wrapper for the Numerical Recipes polynomial interpolation
  //routine polint()

  //x[] is the array of x values

  //y[] is the array of y values

  //'target' is the x value for which you want to interpolate a y value

  //'el' is the number of elements in x[] and y[]

  //x[] must be monotonic, but it can be either ascending or descending

 double *xt,*yt;
 double xint[5],yint[5],err,ans;
 int i,n,min;


 if(x[0]>x[1]){
  xt=new double[el];
  yt=new double[el];
  for(i=0;i<el;i++){
   xt[i]=x[el-1-i];
   yt[i]=y[el-1-i];

  }
 }
 else{
  xt=x;
  yt=y;
 }


 for(i=0;xt[i]<target && i<el;i++);
 if(xt[i]==target){ans=yt[i];}
 else{
  if(i<2)min=0;
  else if(i>el-2)min=el-4;
  else min=i-2;
  for(n=1;n<5;n++){
   xint[n]=xt[min+n-1];
   yint[n]=yt[min+n-1];
  }
  //printf("min was %d %e i %d %e target %e lim %e %e\n",min,xt[min],i,xt[i],target,xt[0],xt[el-1]);
  polint(xint,yint,4,target,&ans,&err);
 }

 if(x[0]>x[1]){
  delete [] xt;
  delete [] yt;
 }

 //if(!(ans>0.0))scanf("%lf",&junk);
 return ans;

}


//the routines below are just the merge sort routine from Numerical Recipes;

int scanner(double *m, int *indices, int index, int el){
/*this will take the matrix m and put everything in it with value
greater than element m[index] to the right of that element
and everything less than m[index] to the left; it then returns
the new index of that anchored element (which you now *know* is in
the right spot*/

	double toobig,toosmall;
	int i,j,n,k,pp;
	int itoobig,itoosmall;
        FILE *output;

	itoobig=0;
	itoosmall=0;
	for(i=0;i<el;i++){
	  if(m[i]<=m[index])itoosmall++;
	
	}
	
	toobig=m[index];
	itoobig=indices[index];
	
	indices[index]=indices[itoosmall-1];
	m[index]=m[itoosmall-1];
	
	index=itoosmall-1;
	m[index]=toobig;
	indices[index]=itoobig;

	i=0;
	j=el-1;
	n=0;
	pp=1;
	while(index<el-1 && n<el-1){
	 for(;m[i]<=m[index] && i<index;i++,n++);
	 itoobig=indices[i];
	 toobig=m[i];
	
	 for(;m[j]>m[index] && j>index;j--,n++){
	
	 }
	 itoosmall=indices[j];
	 toosmall=m[j];
	
	
	
	 if(toosmall<toobig){
	
	 //in case it ran to the end and i found m[index]
	  m[i]=toosmall;
	  indices[i]=itoosmall;
	
	  m[j]=toobig;
	  indices[j]=itoobig;
	 }
	
	}

	return index;
}

void sort(double *m, int *indices, int el){
	double *newm,junk;
	int k,i,*newi,ii,diff;
	
	if(el==2){
	 if(m[0]>m[1]){
	
	   junk=m[1];
	   m[1]=m[0];
	   m[0]=junk;
	
	   k=indices[1];
	   indices[1]=indices[0];
	   indices[0]=k;
	
	 }
	}
	else if(el>2){
	
	diff=0; //just in case all elements are identical
	for(ii=1;ii<el;ii++)if(m[ii]!=m[0])diff=1;
	
	if(diff==1){

	   i=scanner(m,indices,el/2,el);

	   if(i+1<el){
	     newm=&m[i+1];
	    sort(m,indices,i+1);
	     newi=&indices[i+1];
	    sort(newm,newi,el-i-1);
	  }
	  else{
	    sort(m,indices,el-1);
	  }
	
	 }
	
	}
}

void check_sort(double *sorted, double *unsorted, int *inn, int *inn_0, int n){
    int i,j;
    double diff;

    for(i=0;i<n;i++){
        for(j=0;j<n && inn_0[j]!=inn[i];j++);
	
        diff=fabs(sorted[i]-unsorted[j]);
	if(unsorted[inn[i]]!=0.0){
	    diff=diff/fabs(unsorted[inn[i]]);
	}
	if(diff>1.0e-10 || inn_0[j]!=inn[i]){
	    printf("WARNING sort failed to associate\n");
	    exit(1);
	}
    }

    for(i=1;i<n;i++){
        if(sorted[i]<sorted[i-1]){
	    printf("WARNING sort failed to sort\n");
	    exit(1);
	}
    }


}


void sort_and_check(double *list, double *sorted, int *inn, int n){
    int i,*inn_0;

    inn_0=new int[n];
    for(i=0;i<n;i++){
        sorted[i]=list[i];
	inn_0[i]=inn[i];
    }
    sort(sorted,inn,n);
   // printf("sorted\n");
    check_sort(sorted,list,inn,inn_0,n);
   // printf("checked\n");

    delete [] inn_0;
}



void naive_gaussian_solver(
array_1d<double> &aa_in, array_1d<double> &bb_in,
array_1d<double> &xx, int params){


    array_1d<double> buffer,aa,bb;
    buffer.set_dim(params);
    aa.set_dim(params*params);
    bb.set_dim(params);

    buffer.set_name("naive_buffer");
    aa.set_name("naive_aa");
    bb.set_name("naive_bb");

    array_1d<int> dexes;
    dexes.set_dim(params);

    dexes.set_name("naive_dexes");

    int i,k;
    for(i=0;i<params*params;i++){
        aa.set(i,aa_in.get_data(i));

    }
    for(i=0;i<params;i++){
        bb.set(i,bb_in.get_data(i));
        dexes.set(i,i);
    }

    double amax,nn;
    int imax,ii,j;
    int row,col,rowmax,colmax;

    for(ii=0;ii<params;ii++){
        for(row=ii;row<params;row++){
	    for(col=ii;col<params;col++){
	        nn=fabs(aa.get_data(row*params+col));
		if((row==ii && col==ii) || nn>amax){
		
		    amax=nn;
		    rowmax=row;
		    colmax=col;
		
		}
	    }
	}
	
	if(rowmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(ii*params+i));
	    for(i=0;i<params;i++)aa.set(ii*params+i,aa.get_data(rowmax*params+i));
	    for(i=0;i<params;i++)aa.set(rowmax*params+i,buffer.get_data(i));
	
	    nn=bb.get_data(ii);
	    bb.set(ii,bb.get_data(rowmax));
	    bb.set(rowmax,nn);
	
	}
	
	if(colmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(i*params+ii));
	    for(i=0;i<params;i++)aa.set(i*params+ii,aa.get_data(i*params+colmax));
	    for(i=0;i<params;i++)aa.set(i*params+colmax,buffer.get_data(i));
	
	    j=dexes.get_data(ii);
	    dexes.set(ii,dexes.get_data(colmax));
	    dexes.set(colmax,j);
	
	
	}
	
	for(row=ii+1;row<params;row++){
	    nn=aa.get_data(row*params+ii)/aa.get_data(ii*params+ii);

	    for(col=0;col<params;col++){
	        aa.subtract_val(row*params+col,aa.get_data(ii*params+col)*nn);
		
	    }
	
	    bb.subtract_val(row,bb.get_data(ii)*nn);	
	}
	
	/*printf("\n");
	for(i=0;i<params;i++){
	    for(j=0;j<params;j++){
	        printf("%.4e ",aa[i*params+j]);
	    }
	    printf("\n");
	}*/
	
    }

    double err,maxerr,mindiag=-1.0,minfail,mm;
    int ifail;


    maxerr=-1.0;
    for(row=0;row<params;row++){
        for(col=0;col<params;col++){
	    if(row>col){
	        err=fabs(aa.get_data(row*params+col));
		if(err>maxerr)maxerr=err;
	    }
	    else if(row==col){
	        err=fabs(aa.get_data(row*params+col));
		if(mindiag<0.0 || err<mindiag)mindiag=err;
	    }
	}
    }


    /*if(maxerr>1.0e-6 || isnan(maxerr)){
        //printf("tridiagonalization: maxerr %e mindiag %e\n",maxerr,mindiag);
	//exit(1);
    }*/

    for(ii=params-1;ii>=0;ii--){
        buffer.set(ii,bb.get_data(ii));
	for(row=params-1;row>ii;row--){
	    buffer.subtract_val(ii,buffer.get_data(row)*aa.get_data(ii*params+row));
	}
	mm=buffer.get_data(ii)/aa.get_data(ii*params+ii);
	buffer.set(ii,mm);

    }

    for(i=0;i<params;i++){
        xx.set(dexes.get_data(i),buffer.get_data(i));
    }


    for(ii=0;ii<params;ii++){
        nn=0.0;
	for(col=0;col<params;col++){
	    nn+=xx.get_data(col)*aa_in.get_data(ii*params+col);
	}
	
	err=fabs(nn-bb_in.get_data(ii));
	if(bb_in.get_data(ii)!=0.0)err=err/fabs(bb_in.get_data(ii));
	if(err>maxerr || ii==0){
	    maxerr=err;
	    //if(maxerr>1.0e-6)printf("maxerr %e -- %e %e\n",maxerr,nn,bb_in.get_data(ii));
	}
    }


    if(maxerr>1.0e-5 || isnan(maxerr) || isinf(maxerr)){

	
	nn=0.0;
	minfail=-10.0;
	for(i=0;i<params;i++){
	    for(j=i+1;j<params;j++){
	       nn=0.0;
	       for(k=0;k<params;k++){
	           nn+=power(aa_in.get_data(i*params+k)+aa_in.get_data(j*params+k),2);
	       }
	       if(minfail<0.0 || nn<minfail){
	           minfail=nn;
		   ifail=j;
	       }
	    }
	}
	
	throw ifail;
    }


   //printf("naive gaussian solver maxerr %e\n",maxerr);
}


double compare_arr(array_1d<double> &v1, array_1d<double> &v2){
    double err=0.0,maxerr=-1.0;
    int i,dim;

    if(v1.get_dim()!=v2.get_dim()){
        //printf("WARNING in compare_arr the two arrays do not have same dim\n");
	//printf("%d %d\n",v1.get_dim(),v2.get_dim());
        maxerr=exception_value;
    }

    if(v1.get_dim()<v2.get_dim()){
        dim=v1.get_dim();
    }
    else{
        dim=v2.get_dim();
    }

    for(i=0;i<dim;i++){
        err=fabs(v1.get_data(i)-v2.get_data(i));
	if(v1.get_data(i)!=0.0)err=err/fabs(v1.get_data(i));
	if(err>maxerr)maxerr=err;
    }

    return maxerr;

}

int compare_int_arr(array_1d<int> &p1, array_1d<int> &p2){
    /*
    returns 0 if the arrays have different contents (order does not matter)

    returns 1 if they have the same contents
    */
    if(p1.get_dim()!=p2.get_dim()){
        return 0;
    }

    int i,j,ans=1,found_it;

    for(i=0;i<p1.get_dim() && ans==1;i++){
        found_it=0;
        for(j=0;j<p2.get_dim() && found_it==0;j++){
            if(p1.get_data(i)==p2.get_data(j))found_it=1;
        }

        if(found_it==0)ans=0;

    }

    return ans;

}

double bisection(function_wrapper &fn,
                 array_1d<double> &_lowball,
                 array_1d<double> &_highball,
                 double target, double tol,
                 array_1d<double> &found){



   double mu;
   int i;
   array_1d<double> trial,lowball,highball;
   trial.set_name("global_bisection_trial");
   lowball.set_name("global_bisection_lowball");
   highball.set_name("global_bisection_highball");

   for(i=0;i<_lowball.get_dim();i++){
       lowball.set(i,_lowball.get_data(i));
       highball.set(i,_highball.get_data(i));
   }

   double flow=fn(lowball);
   double fhigh=fn(highball);


   if(flow>fhigh){
       mu=flow;
       flow=fhigh;
       for(i=0;i<lowball.get_dim();i++){
           trial.set(i,lowball.get_data(i));
           lowball.set(i,highball.get_data(i));
       }
       fhigh=mu;
       for(i=0;i<lowball.get_dim();i++){
           highball.set(i,trial.get_data(i));
       }

   }

   if(flow>target){
       printf("WARNING in bisection lowball %e target %e\n",
       flow,target);
       exit(1);
   }

   double fout;
   if(target-flow>fhigh-target){
       for(i=0;i<lowball.get_dim();i++){
           found.set(i,highball.get_data(i));
       }
       fout=fhigh;
   }
   else{
       for(i=0;i<lowball.get_dim();i++){
           found.set(i,lowball.get_data(i));
       }
       fout=flow;
   }

   double err=fabs(fout-target);
   int ct;
   double ftrial;
   for(ct=0;ct<200 && err>tol;ct++){
       for(i=0;i<lowball.get_dim();i++){
           trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
       }

       ftrial=fn(trial);
       if(ftrial<target){
           for(i=0;i<lowball.get_dim();i++){
               lowball.set(i,trial.get_data(i));
           }

           if(target-ftrial<err){
               for(i=0;i<lowball.get_dim();i++){
                   found.set(i,trial.get_data(i));
               }
               fout=ftrial;
               err=fabs(target-ftrial);
           }
       }
       else{
           for(i=0;i<lowball.get_dim();i++){
               highball.set(i,trial.get_data(i));
           }

           if(ftrial-target<err){
               for(i=0;i<lowball.get_dim();i++){
                   found.set(i,trial.get_data(i));
               }
               err=fabs(ftrial-target);
           }

       }
   }

   return fout;

}

double integrate_cos_n(double b1, double b2, int n){

    double upper,lower;

    if(n<0){
        printf("integrate_cos_n cannot handle n<0\n");
        exit(1);
    }

    if(n==0){
        return b2-b1;
    }
    else if(n==1){
        return sin(b2)-sin(b1);
    }
    else{
        upper=power(cos(b2),n-1)*sin(b2);
        lower=power(cos(b1),n-1)*sin(b1);
        return (upper-lower)/double(n)+\
                double(n-1)*integrate_cos_n(b1, b2, n-2)/double(n);
    }
}

ellipse_sampler::ellipse_sampler(){
    _dice=NULL;
    _steps=1000;
    _dx=1.0/double(_steps);
    _initialized=0;

    _cos_n_grid.set_name("ellipse_n_cos_grid");
    _dim_record.set_name("ellipse_dim_record");
}

void ellipse_sampler::initialize(int dd, int seed){
    _dim=dd;
    if(_dice!=NULL){
        delete _dice;
    }
    _dice=new Ran(seed);
    _cos_n_grid.set_dim(_dim+1,_steps);
    int i,j;
    double theta;
    for(i=0;i<_cos_n_grid.get_rows();i++){
        for(j=0;j<_steps;j++){
            theta=asin(j*_dx);
            _cos_n_grid.set(i,j,integrate_cos_n(0.0,theta,i));
        }
    }

    _initialized=1;
}

void ellipse_sampler::get_pt(array_1d<double> &out){
    _dim_record.reset_preserving_room();
    double remainder=1.0;
    double sqrt_remainder=1.0;
    int dim_used,dim_left;
    int i_dim,i;
    int i_max_vol;
    double max_vol;
    double roll,coord;
    for(dim_used=0;dim_used<_dim;dim_used++){
        dim_left=_dim-dim_used-1;
        i_dim=-1;
        while(i_dim<0 || _dim_record.contains(i_dim)==1){
            i_dim=_dice->int32()%_dim;
        }
        _dim_record.add(i_dim);
        i_max_vol=int(sqrt_remainder/_dx);
        if(i_max_vol==_cos_n_grid.get_cols()){
            i_max_vol--;
        }
        max_vol=_cos_n_grid.get_data(dim_left,i_max_vol);
        roll=2.0*(_dice->doub()-0.5)*max_vol;
        for(i=0;i<_cos_n_grid.get_cols() && _cos_n_grid.get_data(dim_left,i)<fabs(roll);i++);
        if(i!=0){
            i--;
        }

        coord=10.0;
        while(coord*coord>remainder){
            coord=(i+_dice->doub())*_dx;
        }
        out.set(i_dim,coord);
        if(roll<0.0){
            out.multiply_val(i_dim,-1.0);
        }
        remainder-=coord*coord;
        sqrt_remainder=sqrt(remainder);
    }
}
