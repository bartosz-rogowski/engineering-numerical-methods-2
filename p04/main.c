#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#define PI 4.*atan(1.0)

//przedzial do calkowania
double a = -1.0;
double b = 1.0;


double phi(int i, int j, double ksi)
{
	if(i == 0 && j == 0)
		return 0.5 - 0.75*ksi + 0.25*pow(ksi, 3);
	if(i == 1 && j == 0)
		return 0.5 + 0.75*ksi - 0.25*pow(ksi, 3);
	if(i == 0 && j == 1)
		return 0.25*(1. - ksi - pow(ksi, 2) + pow(ksi, 3));
	if(i == 1 && j == 1)
		return 0.25*(-1. - ksi + pow(ksi, 2) + pow(ksi, 3));
	return 0.0;
}

double w(int i, double ksi_1, double ksi_2)
{
	if(i == 1)
		return 0.25*(1. - ksi_1)*(1. - ksi_2);
	if(i == 2)
		return 0.25*(1. + ksi_1)*(1. - ksi_2);
	if(i == 3)
		return 0.25*(1. + ksi_1)*(1. + ksi_2);
	if(i == 4)
		return 0.25*(1. - ksi_1)*(1. + ksi_2);
	return 0.0;
}

double f(double x, double y)
{
	return sin(2.*y)*pow(sin(x), 2);
}


double second_derivative(double x, int i, int j)
{
	double dx = 0.001;
	return ( phi(i, j, x+dx) - 2.0* phi(i, j, x) + phi(i, j, x-dx) )/pow(dx, 2);
}

// double first_derivative(double x, double dx, int i)
// {
// 	return ( phi(i, x+dx) - phi(i, x-dx) )/(2.0*dx); 
// }


double ksis_to_x_or_y(double* x_or_y, int m, double ksi_1, double ksi_2)
{
	double temp = 0.;
	for(int i = 0; i <= 3; ++i)
	{
		temp += x_or_y[4*(m-1)+i] * w(i+1, ksi_1, ksi_2);
	}
	return temp;
}

int* map_l_to_pairs(int l)
{
	//default for l = 1
	int* temp = calloc(2, sizeof(int));

	if(l == 2)
	{
		temp[0] = 1;
	}

	if(l == 3)
	{
		temp[0] = 1;
		temp[1] = 1;
	}

	if(l == 4)
	{
		temp[1] = 1;
	}
	return temp;
}


double int_F(double* x, double* y, int m, int l1, int i1, int i2, double Jm)
{	
	int* q = map_l_to_pairs(l1); //contains pairs (q[0] = alpha, q[1] = beta)
	double xi_1, wi_1, xi_2, wi_2;
	size_t n = 20;
	gsl_integration_glfixed_table* tab_1 = gsl_integration_glfixed_table_alloc(n);
	gsl_integration_glfixed_table* tab_2 = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int k1 = 0; k1 < n; ++k1)
	{
		gsl_integration_glfixed_point(a,b,k1,&xi_1,&wi_1,tab_1);
		for(int k2 = 0; k2 < n; ++k2)
		{
			gsl_integration_glfixed_point(a,b,k2,&xi_2,&wi_2,tab_2);
			double temp = phi(q[0], i1, xi_1) * phi(q[1], i2, xi_2);
			temp *= f(ksis_to_x_or_y(x, m, xi_1, xi_2), ksis_to_x_or_y(y, m, xi_1, xi_2));
			calka += wi_1 * wi_2 * temp * Jm * (-1.0);
		}
		// calka *= wi_1;
	}
	gsl_integration_glfixed_table_free(tab_1);
	gsl_integration_glfixed_table_free(tab_2);
	free(q);
	return calka;
}

double int_S(double* x, double* y, int m, int l1, int i1, int i2, int l2, int j1, int j2, double Jm)
{	
	int* q_1 = map_l_to_pairs(l1); //contains pairs (q_1[0] = alpha_1, q_1[1] = beta_1)
	int* q_2 = map_l_to_pairs(l2); //contains pairs (q_2[0] = alpha_2, q_2[1] = beta_2)
	double xi_1, wi_1, xi_2, wi_2;
	size_t n = 20;
	gsl_integration_glfixed_table* tab_1 = gsl_integration_glfixed_table_alloc(n);
	gsl_integration_glfixed_table* tab_2 = gsl_integration_glfixed_table_alloc(n);
	double calka = 0.;
	for(int k1 = 0; k1 < n; ++k1)
	{
		gsl_integration_glfixed_point(a,b,k1,&xi_1,&wi_1,tab_1);
		for(int k2 = 0; k2 < n; ++k2)
		{
			gsl_integration_glfixed_point(a,b,k2,&xi_2,&wi_2,tab_2);
			double temp = phi(q_1[0], i1, xi_1) * phi(q_1[1], i2, xi_2);
			temp *= (second_derivative(xi_1, q_2[0], j1) * phi(q_2[1], j2, xi_2) + 
				second_derivative(xi_2, q_2[1], j2) * phi(q_2[0], j1, xi_1) );
			calka += wi_1 * wi_2 * temp;
		}
		// calka *= wi_1;
	}
	gsl_integration_glfixed_table_free(tab_1);
	gsl_integration_glfixed_table_free(tab_2);
	free(q_1);
	free(q_2);
	return calka;
}


void solve(double* A_data, double* b_data, double* x_data, int N)
{
	gsl_matrix_view A = gsl_matrix_view_array(A_data, N, N);
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

	gsl_permutation_free(p);
	gsl_vector_free(x);
}


void MES(int nx, int ny, FILE** fp, int save_to_file)
{
	int M = (nx - 1)*(ny - 1);
	double dx = PI/(nx - 1.);
	double dy = PI/(ny - 1.);

	int N_max = 4*nx*ny;


	double* x = calloc(4*M, sizeof(double));
	double* y = calloc(4*M, sizeof(double));
	double* n = calloc(4*M, sizeof(double));
	double* S = calloc(N_max*N_max, sizeof(double));
	double* F = calloc(N_max, sizeof(double));

	for(int i = 1; i < nx; ++i)
	{
		for(int j = 1; j < ny; ++j)
		{
			int m = i + (j - 1)*(nx - 1) - 1;
			x[4*m+0] = dx*(i-1);
			y[4*m+0] = dy*(j-1);
			n[4*m+0] = i + (j-1)*nx;

			x[4*m+1] = dx*(i);
			y[4*m+1] = dy*(j-1);
			n[4*m+1] = (i+1)+(j-1)*nx;

			x[4*m+2] = dx*(i);
			y[4*m+2] = dy*(j);
			n[4*m+2] = (i+1)+(j+1-1)*nx;

			x[4*m+3] = dx*(i-1);
			y[4*m+3] = dy*(j);
			n[4*m+3] = i+(j+1-1)*nx;
		}
	}

	if(nx == 3)
	{
		for(int m = 0; m < M; ++m)
		{
			printf("m = %d\tnodes: ", m+1);
			for(int temp = 1; temp <= 4; ++temp)
				printf("%g - (%.3f, %.3f),\t", n[4*m+temp-1], x[4*m+temp-1], y[4*m+temp-1]);
			printf("\n");
		}
	}

	for(int m = 1; m <= M; ++m)
	{
		for(int l1 = 1; l1 <= 4; ++l1)
		{
			for(int i1 = 0; i1 <= 1; ++i1)
			{
				for(int i2 = 0; i2 <= 1; ++i2)
				{
					double Jm = 0.25 * (x[4*(m-1)+1] - x[4*(m-1)+0]) * (y[4*(m-1)+3] - y[4*(m-1)+0]);
					
					int p = 4*(n[4*(m-1)+l1-1] - 1) + (i1 + 1) + 2*i2;
					F[p-1] += int_F(x, y, m, l1, i1, i2, Jm);

					for(int l2 = 1; l2 <= 4; l2++)
					{
						for(int j1 = 0; j1 <= 1; j1++)
						{
							for(int j2 = 0; j2 <= 1; j2++)
							{
								int q = 4*(n[4*(m-1)+l2-1] - 1) + (j1 + 1) + 2*j2;
								S[(p-1)*N_max + q-1] += int_S(x, y, m, l1, i1, i2, l2, j1, j2, Jm);
							}	
						}
					}
				}
			}
		}
	}

	for(int m = 1; m <= M; ++m)
	{
		for(int l1 = 1; l1 <= 4; ++l1)
		{
			for(int i1 = 0; i1 <= 1; ++i1)
			{
				for(int i2 = 0; i2 <= 1; ++i2)
				{
					if( (fabs(x[4*(m-1)+l1-1] - 0.) < 0.001 || fabs(x[4*(m-1)+l1-1] - PI) < 0.001 || 
						fabs(y[4*(m-1)+l1-1] - 0.) < 0.001 || fabs(y[4*(m-1)+l1-1] - PI) < 0.001) 
						&& i1*i2 == 0)
					{
						int p = 4*(n[4*(m-1)+l1-1] - 1) + (i1 + 1) + 2*i2;
						F[p-1] = 0.;
						for(int j = 0; j < N_max; ++j)
						{
							S[(p-1)*N_max + j] = S[j*N_max + (p-1)] = 0.;
						}
						S[(p-1)*N_max + (p-1)] = 1.;
					}
				}
			}
		}
	}

	/* //wyswietlanie S i F
	for(int p = 0; p < N_max; ++p)
	{
		for(int q = 0; q < N_max; ++q)
		{
			fprintf(*fp, "%.3g ", S[p*N_max + q]);
		}
		fprintf(*fp, "\n");
	}
	fprintf(*fp, "\n");
	for(int p = 0; p < N_max; ++p)
	{
		fprintf(*fp, "%.3g\n", F[p]);
	} */

	double* S_copy = calloc(N_max*N_max, sizeof(double));
	double* F_copy = calloc(N_max, sizeof(double));
	double* c = calloc(N_max, sizeof(double));
	for(int p = 0; p < N_max; ++p)
	{
		F_copy[p] = F[p];
		for(int q = 0; q < N_max; ++q)
		{
			S_copy[p*N_max + q] = S[p*N_max + q];
		}
	}
	solve(S, F, c, N_max);

	double a_num = 0.;
	for(int i = 0; i < N_max; ++i)
	{
		double temp = 0.0;
		for(int j = 0; j < N_max; ++j)
		{
			temp += 0.5*c[i]*c[j]*S_copy[i*N_max + j];
		}
		a_num += -temp - (-1)*c[i]*F_copy[i];
	}

	printf("\nnx = %d, ny = %d : a_num = %g\n", nx, ny, a_num);

	if(save_to_file == 1)
	{
		for(double xx = 0.0; xx - PI < 0.0001; xx += 0.01)
		{
			for(double yy = 0.0; yy - PI < 0.0001; yy += 0.01)
			{	
				double u_num = 0.;
				int m;
				for(m = 0; m < M; ++m)
				{
					if( xx - x[4*m+0] >= 0.0 && xx - x[4*m+1] < 0.0 &&
						yy - y[4*m+0] >= 0.0 && yy - y[4*m+3] < 0.0)
						break;
				}
				double xi_1 = (xx - ( (x[4*m+0] + x[4*m+1])/2.0 ) ) *(2.0/(x[4*m+1] - x[4*m+0]) );
				double xi_2 = (yy - ( (y[4*m+0] + y[4*m+3])/2.0 ) ) *(2.0/(y[4*m+3] - y[4*m+0]) );
				// printf("For (%.3g, %.3g) m = %d\n", xx, yy, m);
				for(int l1 = 1; l1 <= 4; ++l1)
				{
					int* q = map_l_to_pairs(l1);
					for(int i1 = 0; i1 <= 1; ++i1)
					{
						for(int i2 = 0; i2 <= 1; ++i2)
						{
							int p = 4*(n[4*m+l1-1] - 1) + (i1 + 1) + 2*i2;
							u_num += c[p-1] * phi(q[0], i1, xi_1) *phi(q[1], i2, xi_2);
						}
					}
					free(q);
				}
				fprintf(*fp, "%g %g %g\n", xx, yy, u_num);
			}
			fprintf(*fp, "\n");
		}
	}

	free(S_copy);
	free(F_copy);

	free(c);
	free(x);
	free(y);
	free(n);
	free(S);
	free(F);
}


double exact_solution(double x, double y)
{
  double deno = exp(2.0*PI)-exp(-2.0*PI);
  double first = (1.0/16.0)*(exp(2.0*x)*(exp(-2.0*PI)-1.))/deno;
  double second = -(1.0/16.0)*(exp(-2.0*x)*(exp(2.0*PI)-1.))/deno;
  double third = 1.0/8.0-(1.0/16.0)*cos(2.0*x);

  return sin(2*y)*(first + second + third);
}



int main()
{
	int* q;
	printf("\nTest of mapping l to pairs (alpha, beta):\n");
	for(int i = 1; i <= 4; ++i)
	{
		q = map_l_to_pairs(i);
		printf("%d  --  (%d, %d)\n", i, q[0], q[1]);
		free(q);
	}
	printf("\n");

	FILE* fp_3 = fopen("n3.dat", "w+");
	MES(3, 3, &fp_3, 1);
	MES(5, 5, &fp_3, 0);

	FILE* fp_10 = fopen("n10.dat", "w+");
	MES(10, 10, &fp_10, 1);
	MES(15, 15, &fp_10, 0);
	MES(20, 20, &fp_10, 0);

	FILE* fp = fopen("ndokl.dat", "w+");
	for(double xx = 0.0; xx - PI < 0.0001; xx += 0.01)
	{
		for(double yy = 0.0; yy - PI < 0.0001; yy += 0.01)
		{	
			fprintf(fp, "%g %g %g\n", xx, yy, exact_solution(xx, yy));
		}
		fprintf(fp, "\n");
	}

	fclose(fp_3);
	fclose(fp_10);
	fclose(fp);
	return 0;
}