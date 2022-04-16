#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#define PI 4.*atan(1.0)

//krance przedzialu, w ktorym rozpatrujemy problem
double a = -PI;
double b = PI;

//rozw. analityczne
double u_dokl(double x)
{
	return (x - PI)*(x + PI)*exp(-x);
}

double f(double x)
{
	return 2.0*(1.0 - 3.0*x + x*x - PI*PI)*exp(-x);
}

double L(double x, double dx, int i, double(*baza)(int, double))
{
	double temp = ((*baza)(i, x+dx) - 2.0*(*baza)(i, x) + (*baza)(i, x-dx))/pow(dx, 2);
	temp -= ((*baza)(i, x+dx) - (*baza)(i, x-dx) )/(2.0*dx); 
	return temp;
}

double baza_1(int i, double x)
{
	return cos((i - 0.5)*x)*exp(-x);
}

double baza_2(int i, double x)
{
	return (x - PI)*(x + PI)*pow(x, i-1);
}

//zwraca element a[k][i]
double najmniejsze_kwadraty_a(int k, int i, double(*baza)(int, double))
{	
	double xi, wi;
	double dx = 0.0001;
	size_t n = 40;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int j = 0; j < n; ++j)
	{
		gsl_integration_glfixed_point(a,b,j,&xi,&wi,tab);
		calka += wi*L(xi, dx, k, baza)*L(xi, dx, i, baza);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}

//zwraca element b[k]
double najmniejsze_kwadraty_b(int k, double(*baza)(int, double))
{	
	double xi, wi;
	double dx = 0.0001;
	size_t n = 40;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int j = 0; j < n; ++j)
	{
		gsl_integration_glfixed_point(a,b,j,&xi,&wi,tab);
		calka += wi*L(xi, dx, k, baza)*f(xi);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}

//zwraca element a[k][i]
double galerkin_a(int k, int i, double(*baza)(int, double))
{	
	double xi, wi;
	double dx = 0.0001;
	size_t n = 40;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int j = 0; j < n; ++j)
	{
		gsl_integration_glfixed_point(a,b,j,&xi,&wi,tab);
		calka += wi*(*baza)(k, xi)*L(xi, dx, i, baza);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}

//zwraca element b[k]
double galerkin_b(int k, double(*baza)(int, double))
{	
	double xi, wi;
	double dx = 0.0001;
	size_t n = 40;
	gsl_integration_glfixed_table* tab = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int j = 0; j < n; ++j)
	{
		gsl_integration_glfixed_point(a,b,j,&xi,&wi,tab);
		calka += wi*(*baza)(k, xi)*f(xi);
	}
	gsl_integration_glfixed_table_free(tab);
	return calka;
}


//znajduje rozwiazania Ac = b i zapisuje wyniki do x_data
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


int main()
{
	FILE* f_kolokacja = fopen("kolokacja.dat", "w+");
	FILE* f_kwadraty = fopen("kwadraty.dat", "w+");
	FILE* f_galerkin = fopen("galerkin.dat", "w+");
	double d_x = 0.001;

	double(*bazy[2])(int, double) = {&baza_1, &baza_2};

	//petla po bazach
	for(int baza = 0; baza < 2; ++baza)
	{
		//kolokacja
		for(int N = 6; N < 11; ++N)
		{
			double dx = 2.0*PI/(N+1.0);
			double** a = calloc(N, sizeof(double*));
			double* b = calloc(N, sizeof(double));
			for(int k = 1; k < N+1; ++k)
			{
				a[k-1] = calloc(N, sizeof(double));
				double x_k = k*dx - PI;
				for(int i = 1; i < N+1; ++i)
				{
					a[k-1][i-1] = L(x_k, d_x, i, bazy[baza]);
				}
				b[k-1] = f(x_k);
			}

			double* c = calloc(N, sizeof(double));
			solve(a, b, c, N); //wyniki zapisywane do c

			for(double x = -PI; x - PI < 0.0; x+=0.01)
			{
				double u_n = 0.0;
				for(int k = 0; k < N; ++k)
				{
					u_n += c[k] * bazy[baza](k+1, x);
				}
				fprintf(f_kolokacja, "%g %g\n", x, u_dokl(x) - u_n);
			}

			for(int k = 0; k < N; ++k)
			{
				free(a[k]);
			}
			free(a);
			free(b);
			free(c);
			fprintf(f_kolokacja, "\n\n");
		}



		//najmniejsze kwadraty
		for(int N = 6; N < 11; ++N)
		{
			double dx = 2.0*PI/(N+1.0);
			double** a = calloc(N, sizeof(double*));
			double* b = calloc(N, sizeof(double));
			for(int k = 1; k < N+1; ++k)
			{
				a[k-1] = calloc(N, sizeof(double));
				double x_k = k*dx;
				for(int i = 1; i < N+1; ++i)
				{
					a[k-1][i-1] = najmniejsze_kwadraty_a(k, i, bazy[baza]);
				}
				b[k-1] = najmniejsze_kwadraty_b(k, bazy[baza]);

			}
			
			double* c = calloc(N, sizeof(double));
			solve(a, b, c, N); //wyniki zapisywane do c

			for(double x = -PI; x - PI < 0.0; x+=0.01)
			{
				double u_n = 0.0;
				for(int k = 0; k < N; ++k)
				{
					u_n += c[k]*bazy[baza](k+1, x);
				}
				fprintf(f_kwadraty, "%g %g\n", x, u_dokl(x) - u_n);
			}

			for(int k = 0; k < N; ++k)
			{
				free(a[k]);
			}
			free(a);
			free(b);
			free(c);
			fprintf(f_kwadraty, "\n\n");
		}



		//Galerkin
		for(int N = 6; N < 11; ++N)
		{
			double dx = 2.0*PI/(N+1.0);
			double** a = calloc(N, sizeof(double*));
			double* b = calloc(N, sizeof(double));
			for(int k = 1; k < N+1; ++k)
			{
				a[k-1] = calloc(N, sizeof(double));
				double x_k = k*dx;
				for(int i = 1; i < N+1; ++i)
				{
					a[k-1][i-1] = galerkin_a(k, i, bazy[baza]);
				}
				b[k-1] = galerkin_b(k, bazy[baza]);

			}
					
			double* c = calloc(N, sizeof(double));
			solve(a, b, c, N); //wyniki zapisywane do c

			for(double x = -PI; x - PI < 0.0; x+=0.01)
			{
				double u_n = 0.0;
				for(int k = 0; k < N; ++k)
				{
					u_n += c[k]*bazy[baza](k+1, x);
				}
				fprintf(f_galerkin, "%g %g\n", x, u_dokl(x) - u_n);
			}

			for(int k = 0; k < N; ++k)
			{
				free(a[k]);
			}
			free(a);
			free(b);
			free(c);
			fprintf(f_galerkin, "\n\n");
		}
	}
	fclose(f_kolokacja);
	fclose(f_kwadraty);
	fclose(f_galerkin);
	return 0;
}