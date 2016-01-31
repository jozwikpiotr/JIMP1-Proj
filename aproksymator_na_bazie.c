#include "makespl.h"
#include "piv_ge_solver.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbÄ™ uĹĽywanych f. bazowych moĹĽna ustawiÄ‡ przez wartoĹ›Ä‡
          zmiennej Ĺ›rodowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi(double a, double b, int n, int i, double x)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 1 / h3 * (x - hx[0]) * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (h3 + 3 * h * h * (x - hx[1]) + 3 * h * (x - hx[1]) * (x - hx[1]) - 3 * (x - hx[1]) * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (h3 + 3 * h * h * (hx[3] - x) + 3 * h * (hx[3] - x) * (hx[3] - x) - 3 * (hx[3] - x) * (hx[3] - x) * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return 1 / h3 * (hx[4] - x) * (hx[4] - x) * (hx[4] - x);
}

/* Pierwsza pochodna fi */
double
dfi(double a, double b, int n, int i, double x)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 3 / h3 * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (3 * h * h + 6 * h * (x - hx[1]) - 9 * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (-3 * h * h - 6 * h * (hx[3] - x) + 9 * (hx[3] - x) * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return -3 / h3 * (hx[4] - x) * (hx[4] - x);
}

/* Druga pochodna fi */
double
d2fi(double a, double b, int n, int i, double x)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3 * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (6 * h - 18 * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (6 * h  -18 * (hx[3] - x));
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return 6 / h3 * (hx[4] - x);
}

/* Trzecia pochodna fi */
double
d3fi(double a, double b, int n, int i, double x)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3;
	else if (x > hx[1] && x <= hx[2])
		return -18 / h3;
	else if (x > hx[2] && x <= hx[3])
		return 18 / h3;
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return -6 / h3;
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{
	//---------------------------//
	gsl_matrix  * gsl_eqs = NULL;
	gsl_vector * gsl_v_b = NULL;
	gsl_vector * gsl_v_x = NULL;
	//
	matrix_t       *eqs= NULL;
	//--------------------------//
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );
	//------------------------//
	gsl_eqs = gsl_matrix_alloc (nb, nb);
	gsl_v_b = gsl_vector_calloc(nb);
	gsl_v_x = gsl_vector_calloc(nb);
	//
	eqs = make_matrix(nb, nb + 1);
	//--------------------------//

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	double matr_val;
	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
		{
			for (k = 0; k < pts->n; k++)
			{
				//------------------------//
				matr_val = gsl_matrix_get (gsl_eqs, j, i);
				gsl_matrix_set (gsl_eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]) + matr_val);
				//
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));
				//------------------------//
			}
		}
		for (k = 0; k < pts->n; k++)
		{
			//--------------------------//
			matr_val = gsl_vector_get (gsl_v_b, j);
			gsl_vector_set (gsl_v_b, j, y[k] * fi(a, b, nb, j, x[k]) + matr_val);
			//
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
			//--------------------------//
		}
	}
	
#ifdef DEBUG
	//--------------------//
	gsl_vector_fprintf (stdout, gsl_v_b, "%.5lf");
	putchar('\n');
	//
	write_matrix(eqs, stdout);
	//-------------------//
#endif

	//------------------------//
	gsl_linalg_cholesky_decomp(gsl_eqs);
	gsl_linalg_cholesky_solve(gsl_eqs, gsl_v_b, gsl_v_x);
	//
	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
	//-----------------------//
	
#ifdef DEBUG
	//-----------------------//
	gsl_vector_fprintf (stdout, gsl_v_x, "%.5lf");
	putchar('\n');
	//
	write_matrix(eqs, stdout);
	//-----------------------//
#endif

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				//-----------------------//
				double 		ck = gsl_vector_get (gsl_v_x, k);
				//
				/*double		ck = get_entry_matrix(eqs, k, nb)*/
				//-----------------------//
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							//------------------------//
							yi += gsl_vector_get(gsl_v_x, k) * fi(a, b, nb, k, xi);
							dyi += gsl_vector_get(gsl_v_x, k) * dfi(a, b, nb, k, xi);
							d2yi += gsl_vector_get(gsl_v_x, k) * d2fi(a, b, nb, k, xi);
							d3yi += gsl_vector_get(gsl_v_x, k) * d3fi(a, b, nb, k, xi);
							//
							/*yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);*/
							//------------------------//
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
