//
// Data preparation for geoid fit
//
// 14-08-2018 - v. 3
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>

#ifdef OMP
#include <omp.h>
#endif

#define STR_SIZE	40

using namespace std;

	long double const zer=0.0000000000000000000;
	long double const one=1.0000000000000000000;
//	long double const a=6378137.0;
//	long double const finv=298.257223563;
//	long double const a=10000.0;
//	long double const finv=100.0;
	long double const pi=4.*atan(1.);
	long double const ras=648000./pi;	// Conversion from rad to arcsec
	long double const d2r=pi/180.;	    // Conversion from deg to rad
	long double const nf=360./sqrt(3.);	    // Longitude partition factor
    long double me2,e2,bax,a,finv,dN;
    
//void transf(long double *phi,long double *lam,long double *x,long double *y,long double *z,long double *dell,long int m);
void transf(long double *phi,long double *lam,long double *x,long double *y,long double *z,long int m);

int main()
{
	long int i,j,k,m,sf,n,top;
    int cps=CLOCKS_PER_SEC;
//	long double rmax=one*RAND_MAX;
    string outfi;
	long double f,sp,rp,Ne,h,R;
	long double latstep,lonstep,rl,sp2;
    long double nlat,nlon,cf;
    long double  *phi,*lam,*x,*y,*z;  //,*dell;   // Big arrays

	clock_t t1, t2;

//  Determine parameters
    cout<<"\n Enter semimajor axis a (in m) : ";
    cin>>a;
 	cout<<"\n Enter inverse eccentricity : ";
    cin>>finv;
    cout<<"\n Enter nominal latitude step (in degrees) : ";
    cin>>latstep;
    cout<<"\n Enter random error of cartesian coordinates  (in m) : ";
    cin>>dN;
 
    f=one/finv;
    bax=a*(one-f);
	e2=f*(2.-f);
	me2=(one-e2);

	srand(time(0));

 /*   
//  --  Parameters
	dN=0.0;				// (max-min) range of geoid undulation (m)
	latstep=5.0;		// increment step for latitude, in degrees   (index i) - limit = 0.05
	lonstep=5.0;		// increment step for longitude, in degrees  (index j) - limit = 0.05
	nlat=180/latstep;
	
//  --
*/
   
    sf=(int)(90/latstep);
    sp2=latstep/2.;
    top=12*sf*sf; 
    cout<<"\n\n max total points = "<<top<<endl;
    
   phi   = (long double *)calloc(top, sizeof(long double));
   lam   = (long double *)calloc(top, sizeof(long double));
   x   = (long double *)calloc(top, sizeof(long double));
   y   = (long double *)calloc(top, sizeof(long double));
   z   = (long double *)calloc(top, sizeof(long double));
// dell = (long double *)calloc(top, sizeof(long double));
 
   t1=clock();
   
    m=0;
    
 //  Start from Equator
    n=2*(int)nf;
    nlon=720./n;
    for (i = 0 ; i < n ; i++)   {
        phi[m]=sp2*(i%2);
        lam[m]=0.5*nlon*i;
        m++;
        if (phi[m-1]!=0.0)  {
            phi[m]=-phi[m-1];
            lam[m]=lam[m-1];
            m++;
        } 
    }

    //  Proceed towards North - South
    for (k = 1; k < sf ; k++)   {
        nlat=k*latstep;
        cf=cos(d2r*nlat);
        n=2*(int)(nf*cf);
        lonstep=720./n;
        for (i = 0 ; i < n ; i++)   {
            phi[m]=nlat + sp2*(i%2);
            lam[m]=0.5*lonstep*i;
            m++;
            phi[m]=-phi[m-1];
            lam[m]=lam[m-1];
            m++;
        }
    }
 
 //  Add polar sections
    phi[m]=90.;
    lam[m]=0.;
    m++;
    phi[m]=-90.;
    lam[m]=0.;
    m++;
    
      //  Compute cartesians
//     transf(phi,lam,x,y,z,dell,m);
     transf(phi,lam,x,y,z,m);
     
    //  Write output file
    
     cout<<"\n Enter output file name (with extension) : ";
     cin>>outfi;
     ofstream res(outfi);
     res<<fixed<<"  "<<m;
	 res<<scientific<<setprecision(16)<<endl;
     for (i = 0 ; i < m ; i++) 
         res<<endl<<setw(26)<<phi[i]<<setw(26)<<lam[i]<<setw(26)<<x[i]<<setw(26)<<y[i]<<setw(26)<<z[i];  //<<setw(26)<<dell[i];
  
    res<<endl<<"\n     a = "<<a;
    res<<endl<<"\n finv = "<<finv;
    res<<endl<<"\n Step = "<<latstep;
    res<<endl<<"\n SDev = "<<dN<<" (Random error of cartesian coordinates, in m) \n";
    res<<"\n Max total points = "<<top<<endl;
    
    t2=clock()-t1;
    
    res<<"\n Computing time : "<<float(t2)/cps<<" seconds \n";
    cout<<"\n Computing time : "<<float(t2)/cps<<" seconds \n";
    
    res.close();
    free(phi); free(lam); free(x); free(y); free(z);    

     return 0;
}

//void transf(long double *phi,long double *lam,long double *x,long double *y,long double *z,long double *dell,long int m)
void transf(long double *phi,long double *lam,long double *x,long double *y,long double *z,long int m)
{
    long double d1,d2,cf,dc;
    long int i;
    long double rmax=one*RAND_MAX;

    for (i=0;i<m;i++)  {
        d1=sin(d2r*phi[i]);
        cf=sqrt(one-e2*d1*d1);
        d2=cos(d2r*phi[i]);
        x[i]=(a/cf)*d2*cos(d2r*lam[i]);
        dc=(-0.50+rand()/rmax)*dN;
        x[i]=x[i]+dc;
        y[i]=(a/cf)*d2*sin(d2r*lam[i]);
        dc=(-0.50+rand()/rmax)*dN;
        y[i]=y[i]+dc;
        z[i]=(a/cf)*me2*d1;
        dc=(-0.50+rand()/rmax)*dN;
        z[i]=z[i]+dc;
 /*        
        d1=x[i]/a;
        d1=d1*d1;
        d2=y[i]/a;
        d2=d2*d2;
        d2=d1+d2;
        d1=z[i]/bax;
        d1=d1*d1;
        dell[i]=d2+d1-one;
*/
    }    
}
