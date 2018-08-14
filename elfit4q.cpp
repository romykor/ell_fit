/*
 *  Fit an ellipsoid to geoid undulation data
 *  using a grid of EGM2008 values
 *
 *  V. 4  -  14  August  2018
 */

#include "matfun-q.h"

using namespace std;

int main()
{
	int k,m, cps =CLOCKS_PER_SEC;
	unsigned long i,n;
    double phi,lam, semi_a, finv, gristep,dC;
	const quadfloat one = 1.0000000000000000000000000000;
	//const quadfloat one = expq(0); /* Maybe this is closer, who knows :P */
	quadfloat x, y, z,sR, s0;   // ,dell
	quadfloat U[MATR_DIM], C[MATR_DIM];
	quadfloat **N, **Niv, **Vx;
	quadfloat n_km, u_k; /* helper variables to be used (for safety) in OMP implementation instead of N[m][k] and U[k] respectively */
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
    string nam_in, nam_out; 
    char eq;

	memset(C, 0, MATR_DIM * sizeof(quadfloat));
	memset(U, 0, MATR_DIM * sizeof(quadfloat));

	char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE];  // ,w_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
//	memset(w_s, 0, STR_SIZE * sizeof(char));

    double cof[9];
	cof[0] = pow(10., -11);
	cof[1] = pow(10., -11);
	cof[2] = pow(10., -11);
	cof[3] = pow(10., -11);
	cof[4] = pow(10., -11);
	cof[5] = pow(10., -6);
	cof[6] = pow(10., -6);
	cof[7] = pow(10., -6);
	cof[8] = pow(10., 0);

    //  Define file names
    
    cout<<"\n Enter input data file name (with extension) : ";
    cin>>x_s;
	nam_in  = x_s;
	ifstream dat(nam_in);
	if (!dat.is_open())
		return 1;
    
    cout<<"\n Enter output file name (with extension) : ";
    cin>>y_s;
	nam_out  = y_s;
 
	// Phase 1A  -  Data input  

    dat>>n;         // number of data points
    cout<<"\n Data points : "<<n<<endl;

    quadfloat *d, *r, **a;
	quadfloat *big_unified; /* Workaround to get a contiguous space for all the 'a' matrix, in order to keep spatial locality */
	d = (quadfloat *)calloc(n + 1, sizeof(quadfloat));
	r = (quadfloat *)calloc(n + 1, sizeof(quadfloat));
	big_unified = (quadfloat *)calloc((n + 1) * MATR_DIM, sizeof(quadfloat));
	a = (quadfloat **)calloc(n + 1, sizeof(quadfloat *));

#pragma omp parallel for shared(a, big_unified)
	for (i = 0; i <= n; i++)
		a[i] = &big_unified[i * MATR_DIM];

	// Phase 1B  -  Data input  -  Big simple loop
    
	for (i = 1; i <= n; i++)   {
		dat >>phi>>lam>> x_s >> y_s >> z_s;  //>> w_s;
		/* strtoflt128 converts string to __float128, see documentation of libquadmath */
		x = strtoflt128(x_s, NULL);
		y = strtoflt128(y_s, NULL);
		z = strtoflt128(z_s, NULL);
//		dell = strtoflt128(w_s, NULL);

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
	
	dat >> x_s >>eq>>semi_a;
	dat >> x_s >>eq>>finv;
	dat >> x_s >>eq>>gristep;
	dat >> x_s >>eq>>dC;
    
	dat.close();

	t1 = clock();

	// Phase 2  -  Evaluation of normal symmetric matrix N
	for (k = 0; k < MATR_DIM; k++) {
		memset(&u_k, 0, sizeof(quadfloat));
		for (m = k; m < MATR_DIM; m++)  {
			memset(&n_km, 0, sizeof(quadfloat));
#ifdef OMP
#pragma omp parallel for shared(k, m) reduction(+:n_km, u_k)
#endif
			for (i = 1; i <= n; i++) {		//  again big loop
				n_km += a[i][k]*a[i][m];
				if (m == k)
					u_k += a[i][k]*d[i];
			}
			N[k][m] = n_km;
			N[m][k] = n_km;
			U[k] = u_k;
		}
	}

	// Phase 3  -  Compute inverse matrix  Niv  (Cholesky method)
	cholesky_quad(N, Niv);

	// Phase 4  -  Evaluation of coefficient matrix C
	for (k = 0; k < MATR_DIM; k++)  {
		for (m = 0; m < MATR_DIM; m++)
			C[k] += Niv[k][m]*U[m];
	}

	// Phase 5  -  Evaluation of residuals
    
	memset(&sR, 0, sizeof(quadfloat));
#ifdef OMP
#pragma omp parallel for reduction(+:sR)
#endif
	for (i = 1; i <= n; i++)  {		//  again big loop
		r[i] = -d[i];
		for (m = 0; m < MATR_DIM; m++)
			r[i] += a[i][m]*C[m];
		sR += r[i]*r[i];
	}

	// Phase 6  -  Evaluation of variances
	s0 = sqrtq(sR/(n-MATR_DIM));
	scalma_quad(s0,Niv,Vx);				// Scalar-matrix multiplication

	// Phase 7  -  Transform solution due to scaling
	for (m = 0; m < 9; m++)
		C[m] = C[m]*cof[m];

	t2 = clock() - t1;

	// Phase 8  -  Output of results

	ofstream res(nam_out);
	if (!res.is_open())
		return 2;

	res<<fixed<<setprecision(4)<<"\n Computing time : "<<float(t2)/cps<<" seconds \n";
    res<<" using "<<n<<"  data points  "<<endl;
	res<<fixed<<setprecision(10);
    
    res<<"\n Parameters of model ellipsoid"<<endl;
    res<<"\n Semi-major axis = "<<semi_a<<" m ";
    res<<"\n Inverse eccentricity = "<<finv;
    res<<"\n Nominal data grid step = "<<gristep<<" degrees \n";
    res<<"\n Random error of coordinates = "<<dC<<" m \n";
    
	res<<"\n Polynomial coefficients of ellipsoid \n"<<scientific;
	for (m = 0; m < MATR_DIM; m++) {
		quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", C[m]);
		res<<"\n C["<<m<<"] = "<<setw(25)<<x_s;
	}

    res<<"\n\n Variance-covariance matrix Vx(upper triangular)"<<endl;
	for (k = 0;k < MATR_DIM; k++) {
		for (m = k; m < MATR_DIM; m++) {
			quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", Vx[k][m]);
			res<<setw(26)<<x_s;
		}
		res<<endl;
	}

	quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", s0);
	res<<"\n\n A-posteriori error:  s0 = "<<setw(25)<<x_s<<endl;
    
    res.close();
    free(N); free(Niv); free(Vx);
    free(tmp_qf[0]); free(tmp_qf[1]); free(tmp_qf[2]); 
    free(d); free(r); free(a); free(big_unified);

	return 0;
}

