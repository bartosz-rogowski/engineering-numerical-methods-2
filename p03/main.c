#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#define PI 4.*atan(1.0)

//przedzial do calkowania
double a = -1.0;
double b = 1.0;

double x_max = 6.0;
double d_alpha = 0.05;
double d_ksi = 0.001;

double sign(double x)
{
	double temp = (x > 0.0) ? 1.0 : ( (x < 0.0) ? -1.0 : 0.0 );
	return temp;
}

double phi(int i, double ksi)
{
	if(i == 1)
		return ksi*(ksi - 1.0)/2.0;
	if(i == 2)
		return -(ksi + 1.0)*(ksi - 1.0);
	if(i == 3)
		return ksi*(ksi + 1.0)/2.0;
	return 0.0;
}

double x_ksi(double x_2m_minus_1, double x_2m_plus_1, double ksi)
{
	return 0.5*(x_2m_minus_1 + x_2m_plus_1) + 0.5*(x_2m_plus_1 - x_2m_minus_1)*ksi;
}

double second_derivative(double x, double dx, int i)
{
	return ( phi(i, x+dx) - 2.0* phi(i, x) +  phi(i, x-dx))/pow(dx, 2);
}

double first_derivative(double x, double dx, int i)
{
	return ( phi(i, x+dx) - phi(i, x-dx) )/(2.0*dx); 
}

void eigen_solver(double *S, double *O, double *c, double *E, int N)
{
	gsl_matrix_view A = gsl_matrix_view_array (S, N, N);
	gsl_matrix_view B = gsl_matrix_view_array (O, N, N);
	gsl_vector *eval =  gsl_vector_alloc (N);
	gsl_matrix *evec =  gsl_matrix_alloc (N, N);

	gsl_eigen_gensymmv_workspace *work = gsl_eigen_gensymmv_alloc (N);
	gsl_eigen_gensymmv (&A.matrix, &B.matrix, eval, evec, work);
	gsl_eigen_gensymmv_free (work);
	gsl_eigen_gensymmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

	for(int i = 0; i<N; i++) 
	{
		E[i] = gsl_vector_get(eval,i);
	}

	for(int j = 0; j<N; ++j) 
	{
		for(int i = 0; i<N; ++i)
		{
			c[i*N+j] = gsl_matrix_get(evec, i, j);
		}
	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
}


double int_S(int i, int j, double J_m, double x_2m_minus_1, double x_2m_plus_1)
{	
	double xi, wi;
	size_t n = 5;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int k = 0; k < n; ++k)
	{
		gsl_integration_glfixed_point(a,b,k,&xi,&wi,tab);
		calka += wi*0.5/J_m*first_derivative(xi, d_ksi, i)*first_derivative(xi, d_ksi, j);
		calka += wi*0.5*J_m*pow(x_ksi(x_2m_minus_1, x_2m_plus_1, xi), 2)*phi(i, xi)*phi(j, xi);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}

double int_O(int i, int j, double J_m)
{	
	double xi, wi;
	size_t n = 5;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int k = 0; k < n; ++k)
	{
		gsl_integration_glfixed_point(a,b,k,&xi,&wi,tab);
		calka += wi*J_m*phi(i, xi)*phi(j, xi);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}


void solve(double** A_data, double* b_data, double* x_data, int N)
{
	double* A_data2= calloc(N*N, sizeof(double));
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < N; ++j)
		{
			A_data2[i*N+j] = A_data[i][j];
		}
	}
	gsl_matrix_view A = gsl_matrix_view_array(A_data2, N, N);
	gsl_vector_view b = gsl_vector_view_array(b_data, N);
	gsl_vector* x = gsl_vector_alloc(N);
 
	int s;

	gsl_permutation* p = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(&A.matrix, p, &s);
	gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);

	for(int i = 0; i<N; ++i) 
	{
		x_data[i] = gsl_vector_get(x,i);
	}
	free(A_data2);
	gsl_permutation_free(p);
	gsl_vector_free(x);
}


void MES(int M, FILE** fp)
{
	int N = 2*M+1;
	double* c_kp = calloc(N*N, sizeof(double));
	double* x_wezly = calloc(N, sizeof(double));
	for(double alpha = 0.4; alpha - 2.0 < 0; alpha += d_alpha)
	{
		double* x = calloc(N, sizeof(double));
		double* S = calloc(N*N, sizeof(double));
		double* O = calloc(N*N, sizeof(double));
		for(int k = 1; k <= N; ++k)
		{
			x[k-1] = x_max * pow(fabs(2.0*k-N-1.0)/(N*1.0), alpha) * sign((2.0*k-N-1.0)/(N*1.0));
		}

		for(int m = 1; m <= M; ++m)
		{
			for(int i = 1; i <= 3; ++i)
			{
				for(int j = 1; j <= 3; ++j)
				{
					int p = 2*m + i-2;
					int q = 2*m + j-2;
					// tablica x rozpoczyna sie od zera, stad reindeksacja o -1
					double xa = x[2*m - 1 - 1];
					double xb = x[2*m + 1 - 1];

					double Jm = 0.5*(xb - xa);
					S[(p-1)*N+q-1] += int_S(i, j, Jm, xa, xb);
					O[(p-1)*N+q-1] += int_O(i, j, Jm);
				}
			}
		}

		if(fabs(alpha - 0.4) < 0.01 && M == 5)
		{
			printf("S\n");
			for(int i = 0; i < N; ++i)
			{
				for(int j = 0; j < N; ++j)
					printf("%9.5g ", S[i*N+j]);
			printf("\n");
			}
			printf("\n\nO\n");
			for(int i = 0; i < N; ++i)
			{
				for(int j = 0; j < N; ++j)
					printf("%9.3g ", O[i*N+j]);
			printf("\n");
			}
		}

		double* E = calloc(N, sizeof(double));
		double* c = calloc(N*N, sizeof(double));
		eigen_solver(S, O, c, E, N);

		if(fabs(alpha - 1.4) < 0.01)
		{ //if alpha == 1.4
			for(int p = 0; p < N*N; ++p)
			{
				c_kp[p] = c[p];
				if(p < N)
					x_wezly[p] = x[p];
			}
		}

		if(fabs(alpha - 0.4) < 0.01 && M == 5)
		{
			printf("\n\nE\n");
			for(int i = 0; i < N; ++i)
				printf("%9.5g ", E[i]);
			printf("\n");
		}

		
		fprintf((*fp), "%g %g %g %g %g %g\n", alpha, E[0], E[1], E[2], E[3], E[4]);
		free(O);
		free(S);
		free(x);
		free(E);
		free(c);
	}

	fprintf((*fp), "\n\n");

	//alpha = 1.4
	for(int p = 0; p < 5; ++p)
	{
		for(double x = x_wezly[0]; x - x_wezly[N-1] <= 0.0; x += 0.01)
		{
			for(int m = 1; m <= M; ++m)
			{
				double xa = x_wezly[2*m - 1 - 1];
				double xb = x_wezly[2*m + 1 - 1];
				double Jm = 0.5*(xb - xa);
				if(xa <= x && x < xb)
				{
					double u = 0.0;
					for(int i = 1; i <= 3; ++i)
					{
						int k = 2*m+i-2;
						double xi = (x - 0.5*(xa + xb))/Jm;
						u += c_kp[(k-1)*N+p]*phi(i, xi);
					}
					fprintf((*fp), "%g %g\n", x, u);
				}
			}
		}
		fprintf((*fp), "\n\n");
	}

	free(c_kp);
	free(x_wezly);
}

int fact(int n)
{
	int fact = 1;
	for(int i = 1; i <= n; ++i)
		fact *= i; 
	return fact;
}


double u_dokl(double x, int mu)
{
	double H = 1.0;
	if(mu == 1)
		H = 2.0*x;
	if(mu == 2)
		H = 4.0*pow(x, 2) - 2.0;
	if(mu == 3)
		H = 8.0*pow(x, 3) - 12.0*x;
	if(mu == 4)
		H = 16.0*pow(x, 4) - 48.0*pow(x, 2) + 12.0;
	return H*exp(-0.5*pow(x, 2)) / sqrt(pow(2, mu)*fact(mu)*sqrt(PI));// * sqrt(sqrt(0.25));

}

int main()
{
	FILE* fp5= fopen("M5.dat", "w+");
	FILE* fp10= fopen("M10.dat", "w+");
	FILE* fp30= fopen("M30.dat", "w+");
	FILE* fp_u_dokl = fopen("u_dokl.dat", "w+");

	for(double x = -6.0; x - 6.0 < 0.0; x += 0.05)
	{
		fprintf(fp_u_dokl, "%g %g %g %g %g %g\n", x, u_dokl(x, 0), u_dokl(x, 1), u_dokl(x, 2), u_dokl(x, 3), u_dokl(x, 4));
	}

	MES(5, &fp5);
	MES(10, &fp10);
	MES(30, &fp30);

	fclose(fp5);
	fclose(fp10);
	fclose(fp30);
	fclose(fp_u_dokl);

	return 0;
}