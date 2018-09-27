/*
Gillespie Algorithm for switching process
*/

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include "MersenneTwister.h"
#include "mt.h"

using namespace std;

inline void init();
inline double sum();
inline int select_reaction();
inline void update_p();
inline void update_x();
inline double randn_notrig();

/////////////////////////////////////////////////////////////
// user definition part
//
MTRand *R;
#define N        5		// number of reaction
#define M        3		// number of chemical species
#define cc       5.8321//7.6548//51.8834//437 //47 //189 //378 //189        // number of moelcules at stable point c
#define ca       0.5268//0.8071//3.6042 //20 //5 //20 //40 //20  // number of moelcules at stable point a
#define nn       2              // Hill function order
#define gm       20.0           // gamma_m
#define gpr      10.//10//100.0       // gamma_P_1 / gamma_P_2
#define al       1.
#define K        1.
#define K1       2.*K/(al*gpr)
#define S0       0.11//0.1 //0.11          // S_0
#define Sa       3.3//4.6 //3.3          // S_a
#define WARMUP    20000.0      	// warm up time
#define n_pd      10000
#define Nr       100
#define Nrp      5
 
int x[M];			// population of chemical species
double c[N];	                // coefficients
double p[N];			// propencities of reaction
double er[9];                      // error
double m[9];                      // mean escape time
double Tc[9];
double pd[n_pd];
int pt,pti, dpt, xa, xc, rc, i, j;
double mu, V, sp, tau, r, T;
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
	  //var1=2.0*(1.0 - genrand_real3()) - 1.0;
	  //var2=2.0*(1.0 - genrand_real3()) - 1.0;
	    var1=2.0*R->rand() - 1.0;
	    var2=2.0*R->rand() - 1.0;
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
void init( int x[], double c[], double T, int pt,double V){ //init(x, c, T, pt);
  // population of chemical species
  x[0] = 0;       // mRNA
  x[2] = pt/2.+V*K1/8.-sqrt((pt/2.+V*K1/8.)*(pt/2.+V*K1/8.)-pt*pt/4.) ; 
  x[1] = pt - 2*x[2];      // protein monomer
  // parameters
  c[0] = log(2.0)/T;            //  gamma_p_2
  c[1] = 20.;                 //   k_p
  c[2] = V*S0*K*c[0]*gm/c[1];    //   k_0
  c[3] = V*Sa*K*c[0]*gm/c[1];      //   k_a
}

int select_reaction(double p[], int n, double sp, double r){ //rc = select_reaction(p, N, sp, r);
	int rc = -1;
	double s = 0.0;
	int i;
	r = r * sp;
	for(i=0; i<n; i++){
		s += p[i];
		if(r < s){
			rc = i;
			break;
		}
	}
	return rc;
}

double sum(double p[], int n){  //sp = sum(p, N);
	int i;
	double sp=0.0;
	for(i=0; i<n; i++) 
		sp += p[i];
	return(sp);
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

void update_p( double p[], double c[], int x[], double V){ //update_p(p, c, x);
  //  propensities
  double tmp1 = x[1]/V ;
  tmp1 =tmp1*tmp1 ;
  double tmp2 = K*K ;
  p[0] = c[2]+c[3]*tmp1/(tmp2+tmp1);	
  p[1] = gm*x[0];	
  p[2] = c[1]*x[0];	
  p[3] = gpr*c[0]*x[1];
  p[4] = c[0]*x[2];
} 
void update_x(int x[], int rc, double V){ //update_x(x, rc);
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
        R = new MTRand(); 
	sp = 0.0;	// sum of propencities
	tau=0.0;			// step of time
	//T; //=1. ;
	//for(i=0; i < n_pd ; i++){pd[i]= 0.0 ; }

	// main loop
	
	
	//for(i=0; i < Nrp-1 ; i++){Tc[i]= 0.25 + i * 0.25   ; }
	//Tc[ Nrp-1]= 4. ;
	/*	Tc[0]= 1. ;
	Tc[1]= sqrt(2.) ;
	Tc[2]= 2. ;
	Tc[3]= 2.*sqrt(2.) ;
	Tc[4]= 4. ; */


	for(j=0; j < Nrp; j++){
	  T = 1.; //Tc[j]   ;
	  xa = 20*(j+1); //40 ; //*(j+1) ;
	    xc = floor(cc*xa/ca) ;
	    V = xa/ca ;
	  double mm=0.0;
	  double tt=0.0;
	  for(i=0; i < Nr; i++){  //
	    double ts=0.0;			// time
	    pt = xa ;
	    init(x, c, T, pt,V);
	    while( pt < xc){ 
	    // while(  ts < ENDTIME ){

	      // update propencity
	      update_p(p, c, x,V);
	      sp = sum(p, N);	

	      // sample tau
	      if(sp > 0){
		//tau = -log(1.0 - genrand_real3()) / sp;
		tau = -log(1.0-R->randExc()) / sp;
	      }else{
		break;
	      }
	
	      // select reaction
	      //r = 1.0 -  genrand_real3();
	      r = 1.0-R->randExc() ;
	      rc = select_reaction(p, N, sp, r);

	      /*/ probablity distribution update
	      pt = x[1]+ 2*x[2] ;
		if(ts > WARMUP){
		  pd[pt] += tau; //pd[x[1]] += tau;
		  } /*////prob_dist

	      // update chemical species
	      update_x(x, rc,V);
	      
	      // time
	      ts += tau;
	      pt = x[1]+2*x[2] ;
	    }
	      // output
	    mm += ts ;
	    tt += ts*ts  ;
	  }
	  er[j] = sqrt(tt/Nr-mm*mm/(Nr*Nr))/sqrt(Nr)   ;
	  m[j] = mm/Nr ;
	}
	FILE * fp;
        fp = fopen ("transc.m", "w+");
	fprintf(fp, "function xi = cpp50\n");
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

