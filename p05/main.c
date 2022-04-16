#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
// #define PI 4.*atan(1.0)


double phi(int i, double zeta, double eta)
{
	if(i == 0)
		return -0.5*(eta+zeta);
	if(i == 1)
		return 0.5*(1.+zeta);
	if(i == 2)
		return 0.5*(eta+1.);
	return 0.0;
}


double f(double x, double y)
{
	return exp(-0.5 * (x*x + y*y));
}

double local_to_physical(double* x_or_y, int m, double zeta, double eta)
{
	double result = 0.;
	for(int i = 0; i < 3; ++i)
	{
		result += x_or_y[3*m + i]*phi(i, zeta, eta);
	}
	return result;
}

double diff_zeta(double* x_or_y, int m, double zeta, double eta)
{
	double dx = 0.001;
	double result = 0.;
	for(int i = 0; i < 3; ++i)
	{
		double diff = ( phi(i, zeta+dx, eta) - phi(i, zeta-dx, eta) )/(2.0*dx);
		result += x_or_y[3*m + i]*diff;
	}
	return result;
}

double diff_eta(double* x_or_y, int m, double zeta, double eta)
{
	double dx = 0.001;
	double result = 0.;
	for(int i = 0; i < 3; ++i)
	{
		double diff = ( phi(i, zeta, eta+dx) - phi(i, zeta, eta-dx) )/(2.0*dx);
		result += x_or_y[3*m + i]*diff;
	}
	return result;
}

double Jm(double* x, double* y, int m, double zeta, double eta)
{
	double result = diff_zeta(x, m, zeta, eta) * diff_eta(y, m, zeta, eta);
	result -= diff_zeta(y, m, zeta, eta) * diff_eta(x, m, zeta, eta);
	return result;
}


// double second_derivative(double x, int i, int j)
// {
// 	double dx = 0.001;
// 	return ( phi(i, j, x+dx) - 2.0* phi(i, j, x) + phi(i, j, x-dx) )/pow(dx, 2);
// }

// double first_derivative(double x, double dx, int i)
// {
// 	double dx = 0.001;
// 	return ( phi(i, x+dx) - phi(i, x-dx) )/(2.0*dx); 
// }


double int_F(int m, int i, double* x, double* y)
{
	size_t N = 7;

	double ksi1[] = {
		-0.333333333333333,
		-0.059715871789770,
		-0.059715871789770,
		-0.880568256420460,
		-0.797426985353088,
		-0.797426985353088,
		 0.594853970706174
	};

	double ksi2[] = {
		-0.333333333333333,
		-0.059715871789770,
		-0.880568256420460,
		-0.059715871789770,
		-0.797426985353088,
		 0.594853970706174,
		-0.797426985353088
	};

	double weight[] = {
		0.450000000000000,
		0.264788305577012,
		0.264788305577012,
		0.264788305577012,
		0.251878361089654,
		0.251878361089654,
		0.251878361089654
	};

	double dx = 0.001;
	double intValue = 0.;

	for(size_t k = 0; k < N; k++)
	{
		double J_m = Jm(x, y, m, ksi1[k], ksi2[k]);
		double temp = f(
			local_to_physical(x, m, ksi1[k], ksi2[k]),
			local_to_physical(y, m, ksi1[k], ksi2[k])  
		);
		intValue += weight[k] * J_m * phi(i, ksi1[k], ksi2[k]) * temp;
		
	}

	return intValue;
}



double int_E(int m, int i, int j, double* x, double* y)
{
	size_t N = 7;

	double ksi1[] = {
		-0.333333333333333,
		-0.059715871789770,
		-0.059715871789770,
		-0.880568256420460,
		-0.797426985353088,
		-0.797426985353088,
		 0.594853970706174
	};

	double ksi2[] = {
		-0.333333333333333,
		-0.059715871789770,
		-0.880568256420460,
		-0.059715871789770,
		-0.797426985353088,
		 0.594853970706174,
		-0.797426985353088
	};

	double weight[] = {
		0.450000000000000,
		0.264788305577012,
		0.264788305577012,
		0.264788305577012,
		0.251878361089654,
		0.251878361089654,
		0.251878361089654
	};

	double dx = 0.001;
	double intValue = 0.;

	for(size_t k = 0; k < N; k++)
	{
		double J_m = Jm(x, y, m, ksi1[k], ksi2[k]);
		double* grad_i = calloc(2, sizeof(double));
		double* grad_j = calloc(2, sizeof(double));


		double dphi_dzeta = ( phi(i, ksi1[k]+dx, ksi2[k]) - phi(i, ksi1[k]-dx, ksi2[k]) )/(2.0*dx);
		double dphi_deta  = ( phi(i, ksi1[k], ksi2[k]+dx) - phi(i, ksi1[k], ksi2[k]-dx) )/(2.0*dx);

		grad_i[0] = dphi_dzeta / J_m * diff_eta(y, m, ksi1[k], ksi2[k]);
		grad_i[0] -= dphi_deta / J_m * diff_zeta(y, m, ksi1[k], ksi2[k]);
	
		grad_i[1] = -dphi_dzeta / J_m * diff_eta(x, m, ksi1[k], ksi2[k]);
		grad_i[1] += dphi_deta / J_m * diff_zeta(x, m, ksi1[k], ksi2[k]);
	
		dphi_dzeta = ( phi(j, ksi1[k]+dx, ksi2[k]) - phi(j, ksi1[k]-dx, ksi2[k]) )/(2.0*dx);
		dphi_deta  = ( phi(j, ksi1[k], ksi2[k]+dx) - phi(j, ksi1[k], ksi2[k]-dx) )/(2.0*dx);

		grad_j[0] = dphi_dzeta / J_m * diff_eta(y, m, ksi1[k], ksi2[k]);
		grad_j[0] -= dphi_deta / J_m * diff_zeta(y, m, ksi1[k], ksi2[k]);
	
		grad_j[1] = -dphi_dzeta / J_m * diff_eta(x, m, ksi1[k], ksi2[k]);
		grad_j[1] += dphi_deta / J_m * diff_zeta(x, m, ksi1[k], ksi2[k]);

		intValue += weight[k] * J_m * (grad_i[0]*grad_j[0] + grad_i[1]*grad_j[1]);
		free(grad_i);
		free(grad_j);
	}

	return intValue;
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

double calculate_triangle_area_coordinates(double x1, double x2, double x3, double y1, double y2, double y3)
{
	return 0.5*fabs( x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2) );
}

int find_m(double xx, double yy, double* x, double* y, int M)
{
	int m = 0;
	for(; m < M; ++m)
	{
		double A = calculate_triangle_area_coordinates(x[3*m+0], x[3*m+1], x[3*m+2], y[3*m+0], y[3*m+1], y[3*m+2]);
		double A0 = calculate_triangle_area_coordinates(x[3*m+0], x[3*m+1], xx, y[3*m+0], y[3*m+1], yy);
		double A1 = calculate_triangle_area_coordinates(x[3*m+1], x[3*m+2], xx, y[3*m+1], y[3*m+2], yy);
		double A2 = calculate_triangle_area_coordinates(x[3*m+2], x[3*m+0], xx, y[3*m+2], y[3*m+0], yy);
		if( fabs(A0+A1+A2 - A) < 1e-4 )
			break;
	}
	return m;
}

double* map_xy_to_zeta_eta(double xx, double yy, double* x, double* y, int m)
{
	double* temp = calloc(2, sizeof(double));
	double M = -y[3*m+0]*x[3*m+2] + y[3*m+0]*x[3*m+1] - y[3*m+1]*x[3*m+0] + y[3*m+1]*x[3*m+2] + y[3*m+2]*x[3*m+0] - y[3*m+2]*x[3*m+1];
	double A_zeta = -y[3*m+1]*x[3*m+0] - y[3*m+2]*x[3*m+0] + 2.*yy*x[3*m+0] - 2.*y[3*m+0]*xx + y[3*m+1]*x[3*m+2] + 
		y[3*m+0]*x[3*m+2] - 2.*x[3*m+2]*yy + y[3*m+0]*x[3*m+1] + 2.*xx*y[3*m+2] - y[3*m+2]*x[3*m+1];
	double A_eta = y[3*m+0]*x[3*m+1] + y[3*m+0]*x[3*m+2] - 2.*y[3*m+0]*xx - y[3*m+1]*x[3*m+0] - y[3*m+1]*x[3*m+2] +
		2.*y[3*m+1]*xx - y[3*m+2]*x[3*m+0] + y[3*m+2]*x[3*m+1] + 2.*yy*x[3*m+0] - 2.*yy*x[3*m+1];

	temp[0] = -A_zeta/M;
	temp[1] = A_eta/M;
	return temp;
}


void MES(int nx, int ny)
{
	//podzial na nx*ny podprzedzialow
	int N = (nx+1)*(ny+1);
	int M = 2*nx*ny;
	double x_min, y_min, x_max, y_max;
	x_min = y_min = -5.;
	x_max = y_max = 5.;
	double dx = (x_max - x_min)/nx;
	double dy = (y_max - y_min)/ny;

	double* x = calloc(3*M, sizeof(double));
	double* y = calloc(3*M, sizeof(double));
	int* n = calloc(3*M, sizeof(int));
	double* S = calloc(N*N, sizeof(double));
	double* F = calloc(N, sizeof(double));
	double* c = calloc(N, sizeof(double));

	//9 because of double loop for range(3)
	double* S_loc = calloc(9*M, sizeof(double));
	double* F_loc = calloc(3*M, sizeof(double));

	for(int i = 0; i < nx; ++i)
	{		
		for(int j = 0; j < ny; ++j)
		{
			int m = 2*(i*ny + j);

			x[3*m+0] = x_min + i*dx;
			y[3*m+0] = y_min + j*dy;
			n[3*m+0] = j*(ny+1) + i + 1;

			x[3*m+1] = x_min + (i+1)*dx;
			y[3*m+1] = y_min + j*dy;
			n[3*m+1] = j*(ny+1) + (i+1) + 1;

			x[3*m+2] = x_min + i*dx;
			y[3*m+2] = y_min + (j+1)*dy;
			n[3*m+2] = (j+1)*(ny+1) + i + 1;

			x[3*(m+1)+0] = x_min + (i+1)*dx;
			y[3*(m+1)+0] = y_min + (j+1)*dy;
			n[3*(m+1)+0] = (j+1)*(ny+1) + (i+1) + 1;

			x[3*(m+1)+1] = x_min + i*dx;
			y[3*(m+1)+1] = y_min + (j+1)*dy;
			n[3*(m+1)+1] = (j+1)*(ny+1) + i + 1;

			x[3*(m+1)+2] = x_min + (i+1)*dx;
			y[3*(m+1)+2] = y_min + j*dy;
			n[3*(m+1)+2] = j*(ny+1) + (i+1) + 1;
		}
	}

	//graph of nodes
	FILE* fp = fopen("nodes.dat", "w+");
	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i <= 3; ++i)
		{
			int idx = i%3;
			fprintf(fp, "%g %g %d", x[3*m+idx], y[3*m+idx], n[3*m+idx]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}
	fclose(fp);


	FILE* fp_E = fopen("E.dat", "w+");
	FILE* fp_F = fopen("F.dat", "w+");
	
	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				S_loc[m*9 + 3*i + j] = int_E(m, i, j, x, y);
				fprintf(fp_E, "%g\t", S_loc[m*9 + 3*i + j]);
			}
			F_loc[m*3 + i] = int_F(m, i, x, y);
			fprintf(fp_E, "\n");
			fprintf(fp_F, "%g\n", F_loc[m*3 + i]);
		}
		fprintf(fp_E, "\n");
		fprintf(fp_F, "\n");
	}
	fclose(fp_E);
	fclose(fp_F);

	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i < 3; ++i)
		{
			int p = n[3*m + i] - 1;
			for(int j = 0; j < 3; ++j)
			{
				int q = n[3*m + j] - 1;
				S[p*N + q] += S_loc[m*9 + 3*i + j];
			}
			F[p] += F_loc[m*3 + i];
		}
	}

	 //wyswietlanie S i F
	fp = fopen("SiF.dat", "w+");
	for(int p = 0; p < N; ++p)
	{
		for(int q = 0; q < N; ++q)
		{
			fprintf(fp, "%.3g ", S[p*N + q]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	for(int p = 0; p < N; ++p)
	{
		fprintf(fp, "%.3g\n", F[p]);
	} 
	fclose(fp);


	//warunki brzegowe
	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i < 3; ++i)
		{
			if(fabs(x[3*m+i] - x_min) < 0.001 || fabs(y[3*m+i] - y_min) < 0.001 ||
				fabs(x[3*m+i] - x_max) < 0.001 || fabs(y[3*m+i] - y_max) < 0.001
			)
			{
				int p = n[3*m + i] - 1;
				for(int q = 0; q < N; ++q)
				{
					S[p*N + q] = S[q*N + p] = 0.;
				}
				S[p*N + p] = 1.;
				F[p] = 0.;
			}
		}
	}

	solve(S, F, c, N);

	FILE* fp_num = fopen("wyniki.dat", "w+");
	for(double xx = x_min; xx - x_max < 0.0001; xx += 0.01)
	{
		for(double yy = y_min; yy - y_max < 0.0001; yy += 0.01)
		{	
			double u_num = 0.;
			int m = find_m(xx, yy, x, y, M);
			// printf("%.3g, %.3g\t%d\n", xx, yy, m);
			double* zeta_eta = map_xy_to_zeta_eta(xx, yy, x, y, m);
			for(int i = 0; i < 3; ++i)
			{
				int p = n[3*m + i] - 1;
				u_num += c[p]*phi(i, zeta_eta[0], zeta_eta[1]);
			}
			fprintf(fp_num, "%g %g %g\n", xx, yy, u_num);
			free(zeta_eta);
		}
		fprintf(fp_num, "\n");
	}
	fclose(fp_num);

	free(c);
	free(S_loc);
	free(F_loc);
	free(x);
	free(y);
	free(n);
	free(S);
	free(F);
}



int main()
{

	MES(9, 9);

	return 0;
}