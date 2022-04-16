#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/home/NR/numerical_recipes.c/bessj0.c"
#include "/home/NR/numerical_recipes.c/bessj1.c"

double L = 1.0;
int n = 100;
double TOL = 1e-06;

void f_diff(int l, double d_r, double E, double R[])
{
	for(int i = 1; i < n; ++i)
	{
		R[i+1] = (2.0/pow(d_r, 2) + pow(l, 2)/pow(d_r*i, 2) - 2.0 * E) * R[i] 
			+ (-1.0/pow(d_r, 2) + 1.0/(d_r*i * 2.0 * d_r)) * R[i-1];
		R[i+1] /= (1.0/pow(d_r, 2) + 1.0/(d_r*i * 2.0 * d_r));
	}
}

void Numerov(int l, double d_r, double E, double R[])
{
	double g_i = 2.0*E + ( 1.0 - 4.0*pow(l, 2)/( 4.0 * pow(d_r, 2) ) );
	R[2] = 2.0*(1.0 - 5.0*pow(d_r, 2)/12.0 * g_i) * R[1];
	double g_i_plus_1 = 2.0*E + ( 1.0 - 4.0*pow(l, 2)/( 4.0 * pow(d_r*2.0, 2) ) );
	R[2] /= (1.0 + pow(d_r, 2)/12.0 * g_i_plus_1);
	for(int i = 2; i < n; ++i)
	{
		g_i = 2.0*E + ( 1.0 - 4*pow(l, 2)/( 4.0 * pow(d_r*i, 2) ) );
		double g_i_minus_1 = 2.0*E + ( 1.0 - 4.0*pow(l, 2)/( 4.0 * pow(d_r*(i-1), 2) ) );
		R[i+1] = 2.0*(1.0 - 5.0*pow(d_r, 2)/12.0 * g_i) * R[i]
			- (1.0 + pow(d_r, 2)/12.0 * g_i_minus_1) * R[i-1];

		g_i_plus_1 = 2.0*E + ( 1.0 - 4.0*pow(l, 2)/( 4.0 * pow(d_r*(i+1), 2) ) );
		R[i+1] /= (1.0 + pow(d_r, 2)/12.0 * g_i_plus_1);
	}

	for(int i = 1; i <= n; ++i)
	{
		R[i] /= sqrt(d_r*i);
	}

	//normowanie wartosci
	double norm = 0.0;
	for(int i = 0; i <= n; ++i)
	{
		norm += pow(R[i], 2) * i * pow(d_r, 2);
	}
	for(int i = 0; i <= n; ++i)
	{
		R[i] /= sqrt(norm);
	}
}

int main()
{
	double alfa[2][5] = {
		{0.0, 2.4048, 5.5200, 8.6537, 11.7915},
		{0.0, 3.8317, 7.0155, 10.1734, 13.3236}
	}; //pierwsze wyrazy sa zerami, aby zachowac zgodnosc indeksow z tymi w zadaniu
	double d_r = 0.01;
	double d_E = 0.2;
	double* R = calloc(n+1, sizeof(double));

	R[0] = R[1] = 1.0;
	FILE* fp1 = fopen("roznice.dat", "w+");
	FILE* fp2 = fopen("zad_2.txt", "w+");
	FILE* fp3 = fopen("numerov.dat", "w+");

	fprintf(fp2, "Rozwiazania numeryczne E:\n");

	// --------------------------------------------------- zad 1
	for(double E = d_E; E - 150.1 < 0; E += d_E)
	{
		R[0] = R[1] = 1.0;
		f_diff(0, d_r, E, R);
		fprintf(fp1, "%g %g\n", E, R[n]);
	}

	
	free(R);
	fprintf(fp1, "\n\n");

	// --------------------------------------------------- zad 2
	R = calloc(n+1, sizeof(double));

	//zmienne do rozwiazan numerycznych
	double* E_num = calloc(5, sizeof(double));
	int rozw = 0;

	R[0] = R[1] = 1.0;
	f_diff(0, d_r, 0.0, R);
	double R_n_prev = R[n];
	for(double E = d_E; E - 150.1 < 0; E += d_E)
	{
		f_diff(0, d_r, E, R);
		double Rn = R[n];		
		if(R_n_prev * R[n] < 0)
		{
			double E_k = E;
			double E_k_minus_1 = E - d_E;
			double E_k_plus_1;
			while(1)
			{
				f_diff(0, d_r, E_k, R);
				E_k_plus_1 = E_k - R[n]*(E_k - E_k_minus_1)/(R[n] - R_n_prev);
				if(fabs(E_k_plus_1 - E_k) - TOL < 0.0)
					break;
				E_k_minus_1 = E_k;
				E_k = E_k_plus_1;
			}
			f_diff(0, d_r, E_k_plus_1, R);
			fprintf(fp2, "%g\n", E_k_plus_1);
			E_num[rozw] = E_k_plus_1;
			++rozw;
		}

		R_n_prev = Rn;
	}

	fprintf(fp2, "\n\nRozwiazania analityczne E:\n");
	//wyliczenie dokladnych rozwiazan
	for (int p = 1; p <= 4; ++p)
	{
		double E_a = 0.5*pow(alfa[0][p]/L, 2);
		fprintf(fp2, "%g\n", E_a);
		double* R_dokl = calloc(n+1, sizeof(double));
		double* R_num = calloc(n+1, sizeof(double));
		R_num[0] = R_num[1] = 1.0;
		f_diff(0, d_r, E_num[p-1], R_num);
		fprintf(fp1, "\n\n");
		for(int i = 0; i <= n; ++i)
		{
			R_dokl[i] = bessj0(alfa[0][p]*d_r*i/L);
			fprintf(fp1, "%g %g %g\n", d_r*i, R_dokl[i] - R_num[i], R_num[i]);
		}
		free(R_dokl);
		free(R_num);
	}

	free(R);

	//=============================================================
	//------------------------- Numerov ---------------------------
	//=============================================================

	//----------------------------------------------- zad 1 Numerov
	R = calloc(n+1, sizeof(double));

	R[0] = 0.0;
	R[1] = 1.0;

	for(double E = d_E; E - 150.1 < 0; E += d_E)
	{
		Numerov(1, d_r, E, R);
		fprintf(fp3, "%g %g\n", E, R[n]);
	}

	free(R);
	fprintf(fp3, "\n\n");

	// --------------------------------------------------- zad 2 Numerov
	fprintf(fp2, "\n\nRozwiazania numeryczne E (Numerov):\n");
	R = calloc(n+1, sizeof(double));

	//zmienne do rozwiazan numerycznych
	free(E_num);
	E_num = calloc(5, sizeof(double));
	rozw = 0;

	R[0] = 0.0;
	R[1] = 1.0;
	Numerov(1, d_r, 0, R);
	R_n_prev = R[n];
	for(double E = d_E; E - 150.1 < 0; E += d_E)
	{
		Numerov(1, d_r, E, R);
		double Rn = R[n];
		if(R_n_prev * R[n] < 0)
		{
			double E_k = E;
			double E_k_minus_1 = E - d_E;
			double E_k_plus_1;
			while(1)
			{
				Numerov(1, d_r, E_k, R);
				E_k_plus_1 = E_k - R[n]*(E_k - E_k_minus_1)/(R[n] - R_n_prev);
				if(fabs(E_k_plus_1 - E_k) - TOL < 0.0)
					break;
				E_k_minus_1 = E_k;
				E_k = E_k_plus_1;
			}
			Numerov(1, d_r, E_k_plus_1, R);
			fprintf(fp2, "%g\n", E_k_plus_1);
			E_num[rozw] = E_k_plus_1;
			++rozw;
		}

		R_n_prev = Rn;
	}

	fprintf(fp2, "\n\nRozwiazania analityczne E (Numerov):\n");
	//wyliczenie dokladnych rozwiazan
	for (int p = 1; p <= 4; ++p)
	{
		double E_a = 0.5*pow(alfa[1][p]/L, 2);
		fprintf(fp2, "%g\n", E_a);
		double* R_dokl = calloc(n+1, sizeof(double));
		double* R_num = calloc(n+1, sizeof(double));
		R_num[0] = 0.0;
		R_num[1] = 1.0;
		Numerov(1, d_r, E_num[p-1], R_num);
		fprintf(fp3, "\n\n");
		
		//normowanie wartosci
		double norm = 0.0;
		for(int i = 0; i <= n; ++i)
		{
			R_dokl[i] = bessj1(alfa[1][p]*d_r*i/L);
			norm += pow(R_dokl[i], 2) * i * pow(d_r, 2);
		}

		for(int i = 0; i <= n; ++i)
		{
			fprintf(fp3, "%g %g %g %g\n", d_r*i, R_dokl[i]/sqrt(norm) - R_num[i], R_num[i], R_dokl[i]/sqrt(norm));
		}

		free(R_dokl);
		free(R_num);
	}


	//----------- zwolnienie pamieci i zamkniecie plikow
	free(R);
	free(E_num);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	return 0;
}