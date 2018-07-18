//
//  Fit an ellipsoid to geoid undulation data
//  using a grid of EGM2008 values
//
//  V. 1  -  17 July 2018
//

#include "matfun.h"

using namespace std;

const unsigned long n = 64442;			// actual value may vary
//quadfloat N[MATR_DIM][MATR_DIM],Niv[MATR_DIM][MATR_DIM],Vx[MATR_DIM][MATR_DIM];

int main()
{
	int k,m, cps = 0;
	unsigned long i;
	const quadfloat one=1.000000000000000000000000000000000000000000;
	//const quadfloat one = expq(0); /* Maybe this is closer, who knows :P */
	quadfloat x,y,z,sR,s0;
	quadfloat U[MATR_DIM], C[MATR_DIM];
	quadfloat **N, **Niv, **Vx;
	quadfloat *tmp_qf[3];
	N 	= (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	Niv = (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	Vx 	= (quadfloat **)calloc(MATR_DIM, sizeof(quadfloat *));
	tmp_qf[0] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	tmp_qf[1] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	tmp_qf[2] = (quadfloat *)calloc(MATR_DIM * MATR_DIM, sizeof(quadfloat));
	for (k = 0; k < MATR_DIM; k++) {
		N[k] 	= &tmp_qf[0][k * MATR_DIM];
		Niv[k] 	= &tmp_qf[1][k * MATR_DIM];
		Vx[k] 	= &tmp_qf[2][k * MATR_DIM];
	}

	clock_t t1, t2;

	double cof[9];
	cof[0] = pow(10., 0);
	cof[1] = pow(10., 0);
	cof[2] = pow(10., 0);
	cof[3] = pow(10., 0);
	cof[4] = pow(10., 0);
	cof[5] = pow(10., 0);
	cof[6] = pow(10., 0);
	cof[7] = pow(10., 0);
	cof[8] = pow(10., 0);
	//quadfloat d[n+1],r[n+1],a[n+1][MATR_DIM];
	// Define matrices with "malloc"
	quadfloat *d, *r, **a;
	quadfloat *big_unified; /* Workaround to get a contiguous space for all the 'a' matrix, in order to keep spatial locality */
	d = (quadfloat *)calloc(n + 1, sizeof(quadfloat));
	r = (quadfloat *)calloc(n + 1, sizeof(quadfloat));
	big_unified = (quadfloat *)calloc((n + 1) * MATR_DIM, sizeof(quadfloat));
	a = (quadfloat **)calloc(n + 1, sizeof(quadfloat *));

#pragma omp parallel for shared(a, big_unified)
	for (i = 0; i <= n; i++)
		a[i] = &big_unified[i * MATR_DIM];

//	ifstream dat("geocart.txt");
	ifstream dat("geomod-t.txt");
	if (!dat.is_open())
		return 1;

	t1 = clock();

	// Phase 1  -  Data input  -  Big simple loop
	char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
	for (i = 1; i <= n; i++)   {
		dat >> x_s >> y_s >> z_s;
		/* strtoflt128 converts string to __float128, see documentation of libquadmath */
		x = strtoflt128(x_s, NULL);
		y = strtoflt128(y_s, NULL);
		z = strtoflt128(z_s, NULL);

		d[i] = x*x + y*y + z*z;
		a[i][0] = (d[i]-3. * y*y) * cof[0];
		a[i][1] = (d[i]-3. * z*z) * cof[1];
		a[i][2] = x * y * cof[2];
		a[i][3] = x * z * cof[3];
		a[i][4] = y * z * cof[4];
		a[i][5] = x * cof[5];
		a[i][6] = y * cof[6];
		a[i][7] = z * cof[7];
		a[i][8] = one * cof[8];
	}
	//
	dat.close();

	// Phase 2  -  Evaluation of normal symmetric matrix N
	for (k = 0; k < MATR_DIM; k++) {
		U[k] = 0.;
		for (m = k; m < MATR_DIM; m++)  {
			N[k][m] = 0.;
			for (i = 1; i <= n; i++) {		//  again big loop
				N[k][m] = N[k][m] + a[i][k]*a[i][m];
				if (m == k)
					U[k] = U[k] + a[i][k]*d[i];
			}
			N[m][k] = N[k][m];
		}
	}

	// Phase 3  -  Compute inverse matrix  Niv  (Cholesky method)
	cholesky_quad(N, Niv);

	// Phase 4  -  Evaluation of coefficient matrix C
	for (k = 0; k < MATR_DIM; k++)  {
		C[k] = 0.;
		for (m = 0; m < MATR_DIM; m++)
			C[k] = C[k] + Niv[k][m]*U[m];
	}

	ofstream res("elfit-res.txt");
	if (!res.is_open())
		return 2;

	res<<scientific<<setprecision(10);
	// Phase 5  -  Evaluation of residuals
	sR = 0.;
	for (i = 1;i <= n;i++)  {		//  again big loop
		r[i] = -d[i];
		for (m = 0; m < MATR_DIM; m++)
			r[i] = r[i] + a[i][m]*C[m];
		sR = sR + r[i]*r[i];
	}

	// Phase 6  -  Evaluation of variances
	s0 = sqrtq(sR/(n-MATR_DIM));
	scalma_quad(s0,Niv,Vx);				// Scalar-matrix multiplication

	// Phase 7  -  Transform solution due to scaling
	for (m = 0; m < 9; m++)
		C[m] = C[m]*cof[m];

	t2 = clock() - t1;

	// Phase 8  -  Output of results

	res<<fixed<<setprecision(4)<<"\n Computing time : "<<float(t2)/cps<<" seconds \n";
	res<<scientific<<setprecision(14);
	res<<"\n Polynomial coefficients of ellipsoid"<<endl;
	for (m = 0; m < MATR_DIM; m++)
		res<<"\n C["<<m<<"] = "<<setw(25)<<C[m];

    res<<"\n\n Variance-covariance matrix Vx(upper triangular)"<<endl;
    for (k = 0;k < MATR_DIM; k++) {
    	for (m = k; m < MATR_DIM; m++)
    		res<<setw(26)<<Vx[k][m];
    	res<<endl;
   	}

	res<<"\n\n A-posteriori error:  s0 = "<<setw(25)<<s0<<endl;

	return 0;
}

