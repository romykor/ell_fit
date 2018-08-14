#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <quadmath.h>
#ifdef OMP
#include <omp.h>
#endif

#define MATR_DIM	9
#define STR_SIZE	40
typedef __float128 quadfloat;

// ---------------------------------------
//  Cholesky Inversion and other functions
// ---------------------------------------

using namespace std;

ofstream dot("outmat9.txt");

void display(long double x[9][9])
{
	int i,j,n=9;
    for (i = 0;i < n;i++) {
    	for (j = 0;j < n;j++)
    		dot<<setw(22)<<x[i][j];
    	dot<<endl;
   	}
   	cout<<endl;
}

void cholesky(long double a[MATR_DIM][MATR_DIM],long double b[MATR_DIM][MATR_DIM])
{
    int i,j,n=MATR_DIM,k;
    long double s;
    long double c[MATR_DIM][MATR_DIM],d[MATR_DIM][MATR_DIM],e[MATR_DIM][MATR_DIM];

	dot<<scientific<<setprecision(12);
	dot<<"\n Matrix A \n";
	display(a);

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)  {
			c[i][j] = 0.;
			d[i][j] = 0.;
		}

	for (i = 0;i < n;i++)  {
		s=0.;
		for (k = 0;k < i;k++)   // upper limit reduced
			s=s+c[k][i]*c[k][i];
		c[i][i]=sqrt(a[i][i]-s);
		for (j = i+1;j < n;j++)  {
			s=0.;
			for (k = 0;k < i;k++)
				s=s+c[k][i]*c[k][j];
			c[i][j]=(a[i][j]-s)/c[i][i];
			c[j][i]=0.;
		}
	}

	for (i = 0;i < n;i++)
		d[i][i]=1./c[i][i];
	for (i = 0;i < n-1;i++)  {
		for (j = i+1;j < n;j++)  {
			s=0.;
			for (k = i;k < j;k++)
				s=s+d[i][k]*c[k][j];
			d[i][j]=-d[j][j]*s;
			d[j][i]=0.;
		}
	}

	for (i = 0;i<n;i++)
		for (j = i;j < n;j++)  {
			s=0.;
			for (k = j;k <n ;k++)
				s=s+d[i][k]*d[j][k];
			b[i][j]=s;
			b[j][i]=b[i][j];
		}

	dot<<"\n Matrix B = A^(-1) \n";
	display(b);

//  Check
	for (i = 0;i < n;i++)
		for (j = 0;j < n;j++)  {
			s=0.;
			for (k = 0;k < n;k++)
				s=s+a[i][k]*b[k][j];
			e[i][j]=s;
	}
	dot<<"\n Check!  A*B  (should be Unit matrix) \n";
	display(e);
}

void scalma(long double x,long double **a, long double **b)
{
	int i,j,n=MATR_DIM;
    for (i = 0;i < n;i++)
    	for (j = 0;j < n;j++)
    		b[i][j]=x*a[i][j];
}

void display_quad(quadfloat **x)
{
	int i, j, n = 9;
	char tmp_buff[STR_SIZE];
	memset(tmp_buff, 0, sizeof(tmp_buff));

	for (i = 0;i < n;i++) {
	 	for (j = 0;j < n;j++) {
			quadmath_snprintf(tmp_buff, sizeof(tmp_buff), "%.20Qe", x[i][j]);
			dot <<"   " << tmp_buff;
		}
		dot<<endl;
	}
}

void cholesky_quad(quadfloat **a, quadfloat **b)
{
    int i, j, n = MATR_DIM, k;
 	const quadfloat one = 1.0000000000000000000000000000;
   
    quadfloat s;
    //quadfloat c[MATR_DIM][MATR_DIM], d[MATR_DIM][MATR_DIM], e[MATR_DIM][MATR_DIM];
    quadfloat **c, **d, **e;
	quadfloat *tmp_qf[3];
	c 	= (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	d	= (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	e 	= (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	tmp_qf[0] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	tmp_qf[1] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	tmp_qf[2] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	for (k = 0; k < MATR_DIM; k++) {
		c[k] = &tmp_qf[0][k * MATR_DIM];
		d[k] = &tmp_qf[1][k * MATR_DIM];
		e[k] = &tmp_qf[2][k * MATR_DIM];
	}

	dot << scientific << setprecision(12);
	dot << "\n Matrix A \n";
	display_quad(a);

	for (i = 0; i < n; i++)  {
		s = 0.;
		for (k = 0; k < i; k++)   // upper limit reduced
			s += c[k][i] * c[k][i];
		c[i][i] = sqrtq(a[i][i] - s);
		for (j = i + 1; j < n; j++)  {
			s = 0.;
			for (k = 0; k < i; k++)
				s += c[k][i] * c[k][j];
			c[i][j] = (a[i][j] - s) / c[i][i];
			c[j][i] = 0.;
		}
	}

	for (i = 0; i < n; i++)
		d[i][i] = one / c[i][i]; 
	for (i = 0; i < n - 1; i++)  {
		for (j = i + 1; j < n; j++)  {
			s = 0.;
			for (k = i; k < j; k++)
				s += d[i][k] * c[k][j];
			d[i][j] = -d[j][j] * s;
			d[j][i] = 0.;
		}
	}

	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)  {
			s = 0.;
			for (k = j; k < n ; k++)
				s += d[i][k] * d[j][k];
			b[i][j] = s;
			b[j][i] = b[i][j];
		}

	dot<<"\n Matrix B = A^(-1) \n";
	display_quad(b);

//  Check
#ifdef OMP
#pragma omp parallel for collapse(2) shared(a, b, e) private(s)
#endif
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)  {
			s = 0.;
			for (k = 0; k < n; k++)
				s += a[i][k] * b[k][j];
			e[i][j] = s;
	}
	dot << "\n Check!  A*B  (should be Unit matrix) \n";
	display_quad(e);
    
    dot.close();

    free(c); free(d); free(e);
    free(tmp_qf[0]); free(tmp_qf[1]); free(tmp_qf[2]); 
    
}

void scalma_quad(quadfloat x,quadfloat **a, quadfloat **b)
{
	int i, j, n = MATR_DIM;
#ifdef OMP
#pragma omp parallel for collapse(2) shared(a, b, x)
#endif
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			b[i][j] = x * a[i][j];
}

