/*
 Azim-Berdy Besya 25/06/2014 
 Gillespie Algorithm for birth-death process
  
 species
        X0 is number of mRNA M
	X1 is number of protein P
        X2 protein dimer
	

 reactions
    production mRNA:         0  ---- po  ----------------> X0
    decay mRNA:              X0 ---- p1 = gamma_m M -----> 0 
    production protein:      0  ---- p2 = k_T M ---------> X1 
    decay protein monomer:   X1 ---- p3 = gamma_p1 X1 -----> 0 
    decay protein monomer:   X2 ---- p4 = gamma_p2 X2 -----> 0 

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include "mt.h"


/////////////////////////////////////////////////////////////
// user definition part
//
#define N        5		// number of reaction
#define M        3		// number of chemical species
//#define V        7.99 //6.18 //24.69  //49.39  //24.69     // cell volume
#define cc       5.8321//437 //47 //189 //378 //189        // number of moelcules at stable point c
#define ca       0.5268 //20 //5 //20 //40 //20  // number of moelcules at stable point a
#define nn       2              // Hill function order
#define gm       20.0           // gamma_m
#define gp       10.0       // gamma_P_1 / gamma_P_2
#define K        1.
#define K1       0.2   
#define S0       0.11          // S_0
#define Sa       3.3          // S_a
#define WARMUP    20000.0      	// warm up time
#define n_pd      10000
#define Nr      100
#define Nrp      1
 
int x[M];			// population of chemical species
double c[N];	                // coefficients
double p[N];			// propencities of reaction
double er[9];                      // error
double m[9];                      // mean escape time
double Tc[9];
double pd[n_pd];
int pt,p0, xa, xc;
int pti;
int dpt;
double mu, V;
double sigma;
bool deviateAvailable;
double rsquared;
double var1;
double var2;
double storedDeviate;
double polar;

///  Functions /////////////////////////////////////////////////////////////////
double randn_notrig(double mu, double sigma) {
        
    //        If no deviate has been stored, the polar Box-Muller transformation is
    //        performed, producing two independent normally-distributed random
    //        deviates.  One is stored for the next round, and one is returned.
    if (!deviateAvailable) {
        
        //        choose pairs of uniformly distributed deviates, discarding those
        //        that don't fall within the unit circle
        do {
            var1=2.0*(1.0 - genrand_real3()) - 1.0;
            var2=2.0*(1.0 - genrand_real3()) - 1.0;
            rsquared=var1*var1+var2*var2;
        } while ( rsquared>=1.0 || rsquared == 0.0);
        
        //        calculate polar tranformation for each deviate
        polar=sqrt(-2.0*log(rsquared)/rsquared);
        
        //        store first deviate and set flag
        storedDeviate=var1*polar;
        deviateAvailable=true;
        
        //        return second deviate
        return var2*polar*sigma + mu;
    }
    
    //        If a deviate is available from a previous call to this function, it is
    //        returned, and the flag is set to false.
    else {
        deviateAvailable=false;
        return storedDeviate*sigma + mu;
    }
}
void init( int x[], double c[], double T, int p0){
  // population of chemical species
  pt = p0  ;
  x[0] = 0;       // mRNA
  x[2] = pt/2.+V*K1/8.-sqrt((pt/2.+V*K1/8.)*(pt/2.+V*K1/8.)-pt*pt/4.) ; 
  x[1] = pt - 2*x[2];      // protein monomer
  

  // parameters
  c[0] = gp*log(2.0)/T;            //  gamma_p
  c[1] = 20./T;                 //   k_p
  c[2] = S0*K*c[0]*gm/c[1];    //   k_0
  c[3] = Sa*K*c[0]*gm/c[1];      //   k_a
}

int select_reaction(double p[], int pn, double sp, double r){
	int rc = -1;
	double s = 0.0;
	int i;
	r = r * sp;
	for(i=0; i<pn; i++){
		s += p[i];
		if(r < s){
			rc = i;
			break;
		}
	}
	return rc;
}

double sum(double a[], int n){
	int i;
	double s=0.0;
	for(i=0; i<n; i++) 
		s += a[i];
	return(s);
}

/*/void update_p( double p[], double c[], int x[]){
  //  propensities
  p[0] = V*c[2]+V*c[3]*pow((x[1]/V),nn)/(pow(K,nn)+pow((x[1]/V),nn));	
  p[1] = gm*x[0];	
  p[2] = c[1]*x[0];	
  p[3] = c[0]*x[1];
}
void update_x(int x[], int rc){
  switch(rc){
		case 0: x[0]++; break;
		case 1: x[0]--; break;
		case 2: x[1]++; break;
		case 3: x[1]--; break;
		}
} /*/

void update_p( double p[], double c[], int x[]){
  //  propensities
  p[0] = V*c[2]+V*c[3]*pow((x[1]/V),nn)/(pow(K,nn)+pow((x[1]/V),nn));	
  p[1] = gm*x[0];	
  p[2] = c[1]*x[0];	
  p[3] = c[0]*x[1];
  p[4] = c[0]*x[2]/gp;
} 
void update_x(int x[], int rc){
  pti = x[1] + 2*x[2] ;
  switch(rc){
		case 0: x[0]++; break;
		case 1: x[0]--; break;
		case 2: x[1]++; break;
		case 3: x[1]--; break;
		case 4: x[2]--; break;
		}
   pt = x[1] + 2*x[2] ;
   dpt=pt-pti ;
   if(dpt != 0){
	   mu = pt/2.+V*K1/8.-sqrt((pt/2.+V*K1/8.)*(pt/2.+V*K1/8.)-pt*pt/4.) ;
	   sigma = K1*mu/(4.*pt-8.*mu+8.+K1);
	   x[2] = abs(floor( randn_notrig(mu, sigma) + 0.5 ))  ;
	   x[1] = pt-2*x[2] ;
    }
} // 

	
///   Main Program ////////////////////////////////////////////////////////
int main(void){  
	double sp = 0.0;	// sum of propencities
	double tau=0.0;			// step of time
	double r;			// random number
	int rc;			// reaction number selected
	int i;
	int j;
	double T=1. ;
	//for(i=0; i < n_pd ; i++){pd[i]= 0.0 ; }

	// main loop
	
	//for(i=0; i < 8 ; i++){Tc[i]= 0.25 + i * 0.25   ; }
	//Tc[8]= 4. ;
	for(j=0; j < Nrp; j++){
	  xa = 300 ;//+ 20*(j+1) ;
	V = xa/ca ;
	xc = floor(cc*xa/ca) ;
	  //T = Tc[j]   ;
	  double mm=0.0;
	  double tt=0.0;
	  for(i=0; i < Nr; i++){  //
	    double ts=0.0;			// time
	    p0 = xa ;
	      init(x, c, T, p0);
	    pt = xa ;
	    while( pt < xc){ 
	    // while(  ts < ENDTIME ){

	      // update propencity
	      update_p(p, c, x);
	      sp = sum(p, N);	

	      // sample tau
	      if(sp > 0){
		tau = -log(1.0 - genrand_real3()) / sp;
	      }else{
		break;
	      }
	
	      // select reaction
	      r = 1.0 -  genrand_real3();
	      rc = select_reaction(p, N, sp, r);

	      /*/ probablity distribution update
	      pt = x[1]+ 2*x[2] ;
		if(ts > WARMUP){
		  pd[pt] += tau; //pd[x[1]] += tau;
		  } /*////prob_dist

	      // update chemical species
	      update_x(x, rc);
	      
	      // time
	      ts += tau;
	      pt = x[1]+2*x[2] ;
	    }
	      // output
	    mm += ts/Nr ;
	    tt += ts*ts/Nr  ;
	  }
	  er[j] = sqrt(tt-mm*mm)/sqrt(Nr)   ;
	  m[j] = mm ;
	}
	FILE * fp;
        fp = fopen ("ca1t202.m", "w+");
	fprintf(fp, "function xi = ca1t202\n");
	fprintf(fp, "xi = [\n");
	for(i=0; i < Nrp; i++){
	  fprintf(fp, "%f  %f", m[i],er[i]);
	  fprintf(fp, "\n"); 
	}
	fprintf(fp, "];\n"); 
	fclose(fp); //

	    /*////// saving probability distribution in pd.m
	FILE * fp;
        fp = fopen ("pd.m", "w+");
	fprintf(fp, "function xi = pd\n");
	fprintf(fp, "xi = [\n");
	for(i=0; i < 5000; i++){
	  pd[i] /= (ENDTIME-WARMUP);
	  fprintf(fp, "%f", pd[i]);
	  fprintf(fp, "\n");  
	}
	fprintf(fp, "];\n");  
        fclose(fp);  /*//////save_prob      

	printf("run_F\n");
      
	return(0);
}






