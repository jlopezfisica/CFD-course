/***** Description:
****** This program solves steady conduction problems.
****** Version: Central differences.

****** References: 1. Computational Fluid Dynamics, H.K. Versteeg and W.
******			    Malalasekera, Longman Group Ltd, 1995
******/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables.h"
#include "constants.h"
#include "functions.h"

/* ################################################################# */
int main(void)
/* ################################################################# */
{
	init();     /* initialize variables */
	bound();    /* apply boundary conditions */

	/* calculate coefficients T equation */
	Tcoeff1D(aE, aW, aP, b);

	for (iter = 0; iter <= OUTER_ITER; iter++) {
		/* solve T equation */
		solve1D(T, b, aE, aW, aP);

		/* print intermediate results to standard output */
		printConv();
	} /* for iteration loop */

	output(); /* write output to file */

	return 0;

} /* main */


/* ################################################################# */
void grid(void)
/* ################################################################# */
{
/***** Purpose: Defining the geometrical variables ******/
/*****          See fig. 6.2-6.4 in ref. 1 ******/
	int	I, i;
	double	Dx;

	/* Length of volume element */

	Dx = XMAX/NPI;

	/* Length variable for the scalar points in the x direction */

	x[0] = 0.;
	x[1] = 0.5*Dx;

	for (I = 2; I <= NPI; I++)
		x[I] = x[I-1] + Dx;

	x[NPI+1] = x[NPI] + 0.5*Dx;

	/* Length variable for the grid cell interfaces in the x direction */

	x_u[0] = 0.;
	x_u[1] = 0.;

	for (i = 2; i <= NPI + 1; i++)
		x_u[i] = x_u[i-1] + Dx;

} /* grid */

/* ################################################################# */
void init(void)
/* ################################################################# */
{
/***** Purpose: To initialise all parameters. ******/
	int	I;

	memalloc(); /* allocate memory for variables */
	grid();     /* define geometrical variables */

	/* Initialising all variables  */

	for (I = 0; I <= NPI + 1; I++) {
		Gamma[I] = 10.; /* thermal conductivity */
		T    [I] = 0.; /* Temperature */
	} /* for I */

	/* Definition of first internal node for the different fi variables. */
	/* See fig. 9.1. */

	Istart = 1;
	Iend   = NPI;

} /* init */


/* ################################################################# */
void bound(void)
/* ################################################################# */
{
/***** Purpose: Specify boundary conditions for a calculation ******/

      T[Istart-1] = 100.;
      T[Iend  +1] = 500.;

} /* bound */


/* ################################################################# */
void printConv(void)
/* ################################################################# */
{
/***** Purpose: print convergence to screen ******/

	fprintf(stdout, "%4d T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f T[%d] = %7.3f\n",
	iter, 1, T[1], 2, T[2], 3, T[3], 4, T[4], NPI, T[NPI]);

} /* printConv */

/* ################################################################# */
void output(void)
/* ################################################################# */
{
/***** Purpose: print out result table to file ******/
	int	I;
	FILE *fp;

	if (iter == OUTER_ITER + 1) {
		fp = fopen("output.dat", "w");

		for (I = 0; I <= NPI + 1; I++)
			fprintf(fp, "%11.4e\t%11.4e\n", x[I], T[I]);

		fclose(fp);
	}

} /* output */


/* ################################################################# */
void solve1D(double *fi, double *b, double *aE, double *aW, double *aP)
/* ################################################################# */
{
/***** Purpose: To solve the steady conduction equation ******/

	int	I;

	/* Solving from left to right */

	for (I = Istart; I <= Iend; I++)
		fi[I] = (aE[I]*fi[I+1] + aW[I]*fi[I-1] + b[I])/aP[I];

	/* Solving from right to left */

	for (I = Iend; I >= Istart; I--)
		fi[I] = (aE[I]*fi[I+1] + aW[I]*fi[I-1] + b[I])/aP[I];

} /* solve1D */


/* ################################################################# */
void Tcoeff1D(double *aE, double *aW, double *aP, double *b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/

	int	i, I;
	double	Dw, De, AREAw, AREAe;

	for (I = Istart; I <= Iend; I++) {
		i = I;

		/* Geometrical parameters */
		/* Areas of the cell faces */

		AREAw = 1.;
		AREAe = 1.;

		/* The diffusion conductance D=(Gamma/Dx)*AREA defined in eq. 5.8b  */

		Dw = ((Gamma[I-1] + Gamma[I]  )/(2*(x[I]   - x[I-1])))*AREAw;
		De = ((Gamma[I]   + Gamma[I+1])/(2*(x[I+1] - x[I]  )))*AREAe;

		/* The source terms */

		SP[I] = 0.;
		Su[I] = 0.;

		/* The coefficients (central differencing sheme) */

		aW[I] = Dw;
		aE[I] = De;

		aP[I] = aW[I] + aE[I] - SP[I];

		/* Setting the source term b = Su */

		b[I] = Su[I];

		} /* for I */

} /* Tcoeff1D */


/* ################################################################# */
int *int_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type int */
	int *a;

	a = (int *) calloc(np, sizeof(int));

	return a;

} /* int_1D_array */

/* ################################################################# */
double *double_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type double */
	double *a;

	a = (double *) calloc(np, sizeof(double));

	return a;

} /* double_1D_array */


/* ################################################################# */
void memalloc(void)
/* ################################################################# */
{
/***** Purpose: Allocating memory for dependent variables. ******/

	x    = double_1D_array(NPI + 2);
	x_u  = double_1D_array(NPI + 2);

	T      = double_1D_array(NPI + 2);
	Gamma  = double_1D_array(NPI + 2);

	aP     = double_1D_array(NPI + 2);
	aE     = double_1D_array(NPI + 2);
	aW     = double_1D_array(NPI + 2);
	b      = double_1D_array(NPI + 2);

	SP     = double_1D_array(NPI + 2);
	Su     = double_1D_array(NPI + 2);

} /* memalloc */
