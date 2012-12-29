// classes defines the growth factor
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_integration.h>

#define MAXSTP 10000
#define TINY 1.0e-30
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))
static double maxarg1,maxarg2, minarg1, minarg2;

//void derivs (double a, double b[], double c[]);

using namespace std;

class growth {

 protected:
  float a, H0, Omega_M, Omega_k, Omega_L, wt;
  double g, f;

 public:
   
 growth(float inp1, float inp3, float inp4, cosmo & cosm_model) {
    // cosmological params at z =0
    Omega_M = inp1;
    Omega_k = inp3;
    Omega_L = 1.0- Omega_M- Omega_k;
    wt= inp4;
 }

/*growth factor */
  

void growthfactor(float redshift, double *Da, double *Dadz){
   double x1, x2, dplus, ddot;
   const float redshift1=200.0;
   double *ystart;
   
   x1 = 1.0/(1.0+100000.0);
   x2 = 1.0/(1.0+redshift);
   ystart = (double *)malloc(2*sizeof(double));
   ystart[0] = x1;
   ystart[1] = 0.0;
      
   odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, false);
      //printf("Dplus = %f;  Ddot = %f \n", ystart[0], ystart[1]);	
      
   dplus = ystart[0];
   ddot  = ystart[1];
   x1 = 1.0/(1.0+100000.0);
   x2 = 1.0;
   ystart[0] = x1;
   ystart[1] = 0.0;
      
   odesolve(ystart, 2, x1, x2, 1.0e-6, 1.0e-6, false);	
      //printf("Dplus = %f;  Ddot = %f \n", ystart[0], ystart[1]);
      
   *Da    = dplus/ystart[0];
   *Dadz = ddot/ystart[0];
   //printf("\nredshift= %lf Growth factor = %f;  Derivative = %f \n", redshift, *Da, *Dadz);
   free(ystart);
   
   return;
}


void odesolve(double ystart[], int nvar, double x1, double x2, double eps, double h1, bool print_stat){	
   int i, nstp, nok, nbad, feval;
   double x,hnext,hdid,h;
   double *yscal,*y,*dydx;
   const double hmin=0.0;
      
   feval = 0;
   yscal= (double *)malloc(nvar*sizeof(double));
   y= (double *)malloc(nvar*sizeof(double));
   dydx= (double *)malloc(nvar*sizeof(double));
   x=x1;
   h=SIGN(h1,x2-x1);
   nok = nbad = 0;
   for (i=0; i<nvar; ++i) {y[i]=ystart[i];}
      
   for (nstp=0; nstp<MAXSTP; ++nstp) {
      derivs(x, y, dydx);
      ++feval;
      for (i=0; i<nvar; ++i)
      {yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;}
      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,&feval);
      if (hdid == h) ++nok; else ++nbad;
      if ((x-x2)*(x2-x1) >= 0.0) {
         for (i=0; i<nvar; ++i) {ystart[i]=y[i];}
         free(dydx);
         free(y);
         free(yscal);
         if (print_stat){
            printf("ODEsolve:\n");
            printf(" Evolved from x = %f to x = %f\n", x1, x2);
            printf(" successful steps: %d\n", nok);
            printf(" bad steps: %d\n", nbad);
            printf(" function evaluations: %d\n", feval);
         }				
         return;
      }
      if (fabs(hnext) <= hmin) {
         printf("Step size too small in ODEsolve");
         exit(1);
      }
      h=hnext;
   }
   printf("Too many steps in ODEsolve");
   exit(1);
}



void derivs(double a, double y[], double dydx[]){
   double H;
   // float Omega_M, Omega_L, wt;
   H = sqrt(Omega_M/(a*a*a) + (Omega_L)*pow(a, -3.0*(1.0+wt)));
   dydx[0] = y[1]/(a*H);
   dydx[1] = -2.0*y[1]/a + 1.5*Omega_M*y[0]/(H*pow(a, 4.0));
   return;
}


void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps, 
                      double yscal[], double *hdid, double *hnext, int *feval)
{
   int i;
   double errmax,h,htemp,xnew,*yerr,*ytemp;
      
   yerr= (double *)malloc(n*sizeof(double));
   ytemp= (double *)malloc(n*sizeof(double));
   h=htry;
      
   for (;;) {
      rkck(y,dydx,n,*x,h,ytemp,yerr);
      *feval += 5;
      errmax=0.0;
      for (i=0; i<n; ++i) {errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));}
      errmax /= eps;
      if (errmax <= 1.0) break;
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
      xnew=(*x)+h;
      if (xnew == *x) {
         printf("Stepsize underflow in ODEsolve rkqs");
         exit(1);
      }
   }
   if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
   else *hnext=5.0*h;
   *x += (*hdid=h);
   for (i=0; i<n; ++i) {y[i]=ytemp[i];}
   free(ytemp);
   free(yerr);
}


/* Cash-Karp Runge-Kutta Step, based on the 
work of Fehlberg, who ﬁgured out that six function evaluations could 
be used to determine two Runge-Kutta steps, one fourth-order and one 
ﬁfth-order. The diﬀerence between these estimates can be used as a 
truncation error for adjusting the stepsize. */
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
                      double yerr[])
{
   int i;
   static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
   b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
   b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
   b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
   b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
   c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
   dc5 = -277.00/14336.0;
   double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
   dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
      
   ak2= (double *)malloc(n*sizeof(double));
   ak3= (double *)malloc(n*sizeof(double));
   ak4= (double *)malloc(n*sizeof(double));
   ak5= (double *)malloc(n*sizeof(double));
   ak6= (double *)malloc(n*sizeof(double));
   ytemp= (double *)malloc(n*sizeof(double));
      
   for (i=0; i<n; ++i)
      ytemp[i]=y[i]+b21*h*dydx[i];
   derivs(x+a2*h,ytemp,ak2);
   for (i=0; i<n; ++i)
      ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
   derivs(x+a3*h,ytemp,ak3);
   for (i=0; i<n; ++i)
      ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
   derivs(x+a4*h,ytemp,ak4);
   for (i=0; i<n; ++i)
      ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
   derivs(x+a5*h,ytemp,ak5);
   for (i=0; i<n; ++i)
      ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
   derivs(x+a6*h,ytemp,ak6);
   for (i=0; i<n; ++i)
      yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
   for (i=0; i<n; ++i)
      yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
   
   free(ytemp);
   free(ak6);
   free(ak5);
   free(ak4);
   free(ak3);
   free(ak2);
}



};
  
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef FMAX
#undef FMIN
#undef MAXSTP
#undef TINY
#undef SIGN


 
 
