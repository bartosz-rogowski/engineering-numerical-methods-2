#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#define PI 4.*atan(1.0)


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


double int_O(int m, int i, int j, double* x, double* y)
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
		intValue += weight[k] * J_m * phi(i, ksi1[k], ksi2[k]) * phi(j, ksi1[k], ksi2[k]);
		
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


void eigen_solver(gsl_matrix *E, gsl_matrix *O, gsl_matrix *c, gsl_vector *lambda)
{
 	const int N = E->size1;

	gsl_eigen_gensymmv_workspace *work = gsl_eigen_gensymmv_alloc(N);
	gsl_eigen_gensymmv(E, O, lambda, c, work);
	gsl_eigen_gensymmv_free(work);
	gsl_eigen_gensymmv_sort(lambda, c, GSL_EIGEN_SORT_VAL_ASC);
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
	double* E = calloc(N*N, sizeof(double));
	double* O = calloc(N*N, sizeof(double));

	//9 because of double loop for range(3)
	double* E_loc = calloc(9*M, sizeof(double));
	double* O_loc = calloc(9*M, sizeof(double));

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
	FILE* fp_O = fopen("O.dat", "w+");
	
	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i < 3; ++i)
		{
			for(int j = 0; j < 3; ++j)
			{
				E_loc[m*9 + 3*i + j] = int_E(m, i, j, x, y);
				fprintf(fp_E, "%g\t", E_loc[m*9 + 3*i + j]);

				O_loc[m*9 + 3*i + j] = int_O(m, i, j, x, y);
				fprintf(fp_O, "%g\t", O_loc[m*9 + 3*i + j]);
			}
		}
		fprintf(fp_E, "\n");
		fprintf(fp_O, "\n");
	}
	fclose(fp_E);
	fclose(fp_O);

	for(int m = 0; m < M; ++m)
	{
		for(int i = 0; i < 3; ++i)
		{
			int p = n[3*m + i] - 1;
			for(int j = 0; j < 3; ++j)
			{
				int q = n[3*m + j] - 1;
				E[p*N + q] += E_loc[m*9 + 3*i + j];
				O[p*N + q] += O_loc[m*9 + 3*i + j];
			}
		}
	}


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
					E[p*N + q] = E[q*N + p] = 0.;
					O[p*N + q] = O[q*N + p] = 0.;
				}
				E[p*N + p] = 2000.;
				O[p*N + p] = 1.;
			}
		}
	}

	gsl_matrix* O_gsl = gsl_matrix_calloc(N, N);
	gsl_matrix* O_gsl_copy = gsl_matrix_calloc(N, N);
	gsl_matrix* E_gsl_copy = gsl_matrix_calloc(N, N);
	
	gsl_matrix* E_gsl = gsl_matrix_calloc(N, N);
	gsl_matrix* c     = gsl_matrix_calloc(N, N);
	gsl_vector* lambda = gsl_vector_calloc(N);

	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < N; ++j)
		{
			gsl_matrix_set(E_gsl, i, j, E[i*N+j]);
			gsl_matrix_set(O_gsl, i, j, O[i*N+j]);
		}
	}

	free(E);
	free(O);

	gsl_matrix_memcpy(O_gsl_copy, O_gsl);
	gsl_matrix_memcpy(E_gsl_copy, E_gsl);

	eigen_solver(E_gsl, O_gsl, c, lambda);

	gsl_matrix_memcpy(O_gsl, O_gsl_copy);
	gsl_matrix_memcpy(E_gsl, E_gsl_copy);
	gsl_matrix_free(O_gsl_copy);
	gsl_matrix_free(E_gsl_copy);


	FILE* fp_num = fopen("wyniki.dat", "w+");
	for(int c_count = 0; c_count < 10; ++c_count)
	{
		for(double xx = x_min; xx - x_max < 0.0001; xx += 0.04)
		{
			for(double yy = y_min; yy - y_max < 0.0001; yy += 0.04)
			{	
				double u_num = 0.;
				int m = find_m(xx, yy, x, y, M);
				// printf("%.3g, %.3g\t%d\n", xx, yy, m);
				double* zeta_eta = map_xy_to_zeta_eta(xx, yy, x, y, m);
				for(int i = 0; i < 3; ++i)
				{
					int p = n[3*m + i] - 1;
					u_num += gsl_matrix_get(c, p, c_count)*phi(i, zeta_eta[0], zeta_eta[1]);
				}
				fprintf(fp_num, "%g %g %g\n", xx, yy, u_num);
				free(zeta_eta);
			}
			fprintf(fp_num, "\n");
		}
		fprintf(fp_num, "\n\n");
	}
	fclose(fp_num);

	FILE* fp_lambdy = fopen("lambdy.dat", "w+");
	for(int i = 0; i < 10; ++i)
	{
		fprintf(fp_lambdy, "%g\n", gsl_vector_get(lambda, i));
	}
	fclose(fp_lambdy);


	for(int i = 0; i < N; ++i)
	{	
		double norm;
		gsl_vector_view column = gsl_matrix_column(c, i);
		gsl_vector* temp = gsl_vector_calloc(N);
		// gsl_vector_memcpy(c, temp);
		gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &column.vector, 0.0, temp);
		gsl_blas_ddot(&column.vector, temp, &norm);
		for(int j = 0; j < N; ++j)
		{
			double value = gsl_matrix_get(c, j, i)/sqrt(norm);
			gsl_matrix_set(c, j, i, value);
		}
		gsl_vector_free(temp);
	}


	//------------------------------------- schemat Newmarka

	gsl_vector_view y_k = gsl_matrix_column(c, 1);
	gsl_vector* y_k_next = gsl_vector_calloc(N);
	gsl_vector_view v_k = gsl_matrix_column(c, 2);
	gsl_vector_scale(&v_k.vector, sqrt( gsl_vector_get(lambda, 2) ) );
	gsl_vector* v_k_next = gsl_vector_calloc(N);

	double t_max = 2.0*PI/sqrt(gsl_vector_get(lambda, 1) );
	double dt = t_max/1e4;
	double beta = 0.25;
	int iterations = ceil(t_max/dt);

	
	gsl_matrix* A = gsl_matrix_calloc(N, N);
	gsl_matrix* O_LLT = gsl_matrix_calloc(N, N);
	gsl_matrix* E_copy = gsl_matrix_calloc(N, N);
	gsl_matrix_memcpy(E_copy, E_gsl);
	gsl_matrix_scale(E_copy, beta*dt*dt);
	gsl_matrix_memcpy(A, O_gsl);
	gsl_matrix_memcpy(O_LLT, O_gsl);
	gsl_matrix_add(A, E_copy);
	gsl_linalg_cholesky_decomp1(A);
	gsl_linalg_cholesky_decomp1(O_LLT);

	FILE* fp_vals = fopen("vals.dat", "w+");
	fp_num = fopen("maps.dat", "w+");

	// gsl_vector* temp_c3 = gsl_vector_calloc(N);
	// double val_c3 = 0.;
	// gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &y_k.vector, 0.0, temp_c3);
	// gsl_blas_ddot(&y_k.vector, temp_c3, &val_c3);
	// fprintf(fp_vals, "%g %g %g %g %g\n", 0.0, 1.0, 0.0, val_c3, 0.0);
	// gsl_vector_free(temp_c3);


	for(int it = 1; it <= iterations; ++it)
	{
		gsl_vector* temp_1 = gsl_vector_calloc(N);
		gsl_vector* temp_2 = gsl_vector_calloc(N);
		gsl_vector* temp_3 = gsl_vector_calloc(N);

		gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &y_k.vector, 0.0, temp_1);
		gsl_blas_dgemv(CblasNoTrans, dt, O_gsl, &v_k.vector, 0.0, temp_2);
		gsl_blas_dgemv(CblasNoTrans, 0.5*(2.0*beta - 1.0)*dt*dt, E_gsl, &y_k.vector, 0.0, temp_3);

		gsl_vector_add(temp_2, temp_1);
		gsl_vector_add(temp_3, temp_2);

		gsl_linalg_cholesky_solve(A, temp_3, y_k_next);

		//-------------------------------------

		gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &v_k.vector, 0.0, temp_1);
		gsl_blas_dgemv(CblasNoTrans, -dt/2.0, E_gsl, &y_k.vector, 0.0, temp_2);
		gsl_blas_dgemv(CblasNoTrans, -dt/2.0, E_gsl, y_k_next, 0.0, temp_3);
		
		gsl_vector_add(temp_2, temp_1);
		gsl_vector_add(temp_3, temp_2);

		gsl_linalg_cholesky_solve(O_LLT, temp_3, v_k_next);


		gsl_vector_free(temp_1);
		gsl_vector_free(temp_2);
		gsl_vector_free(temp_3);

		gsl_vector_memcpy(&y_k.vector, y_k_next);
		gsl_vector_memcpy(&v_k.vector, v_k_next);

		if(it % 100 == 0)
		{
			double val_1, val_2, val_3, val_4;
			val_1 = val_2 = val_3 = val_4 = 0.;

			gsl_vector* temp_1 = gsl_vector_calloc(N);
			gsl_vector* temp_2 = gsl_vector_calloc(N);
			gsl_vector* temp_3 = gsl_vector_calloc(N);
			gsl_vector* temp_4 = gsl_vector_calloc(N);

			gsl_vector_view c_2 = gsl_matrix_row(c, 1);
			gsl_vector_view c_3 = gsl_matrix_row(c, 2);

			gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &c_2.vector, 0.0, temp_1);
			gsl_blas_ddot(&y_k.vector, temp_1, &val_1);
			

			gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &c_3.vector, 0.0, temp_2);
			gsl_blas_ddot(&y_k.vector, temp_2, &val_2);
			

			gsl_blas_dgemv(CblasNoTrans, 1.0, O_gsl, &y_k.vector, 0.0, temp_3);
			gsl_blas_ddot(&y_k.vector, temp_3, &val_3);
			

			gsl_blas_dgemv(CblasNoTrans, 1.0, E_gsl, &y_k.vector, 0.0, temp_4);
			gsl_blas_ddot(&y_k.vector, temp_4, &val_4);
			fprintf(fp_vals, "%g %g %g %g %g\n", dt*it, val_1, val_2, val_3, val_4);

			gsl_vector_free(temp_1);
			gsl_vector_free(temp_2);
			gsl_vector_free(temp_3);
			gsl_vector_free(temp_4);
		}

		
		if(it%1000 == 0)
		{
			for(double xx = x_min; xx - x_max < 0.0001; xx += 0.04)
			{
				for(double yy = y_min; yy - y_max < 0.0001; yy += 0.04)
				{	
					double u_num = 0.;
					int m = find_m(xx, yy, x, y, M);
					// printf("%.3g, %.3g\t%d\n", xx, yy, m);
					double* zeta_eta = map_xy_to_zeta_eta(xx, yy, x, y, m);
					for(int i = 0; i < 3; ++i)
					{
						int p = n[3*m + i] - 1;
						u_num += gsl_vector_get(&y_k.vector, p)*phi(i, zeta_eta[0], zeta_eta[1]);
					}
					fprintf(fp_num, "%g %g %g\n", xx, yy, u_num);
					free(zeta_eta);
				}
				fprintf(fp_num, "\n");
			}
			fprintf(fp_num, "\n\n");
		}
	}

	fclose(fp_num);
	fclose(fp_vals);

	//------------------------------------------------------
	// gsl_vector_free(y_k);
	gsl_vector_free(y_k_next);
	// gsl_vector_free(v_k);
	gsl_vector_free(v_k_next);

	gsl_matrix_free(E_gsl);
	gsl_matrix_free(O_gsl);
	gsl_matrix_free(O_LLT);
	gsl_matrix_free(c);
	gsl_vector_free(lambda);
	gsl_matrix_free(A);
	gsl_matrix_free(E_copy);

	free(E_loc);
	free(O_loc);
	free(x);
	free(y);
	free(n);
}



int main()
{

	MES(9, 9);

	return 0;
}