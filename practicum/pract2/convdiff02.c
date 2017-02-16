/***** Solves: Steady, compressible convection-diffusion problems.

****** Description:
****** This program solves steady convection-diffusion problems	
****** using the simple algorithm described in ch. 6.4 in "Computational 
****** Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
****** variables follow exactly the notations in this reference, and all 
****** equations cited are from this reference unless mentioned otherwise.

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
int main(int argc, char *argv[])
/* ################################################################# */
{
	int	iter_T;
	
	init();
	bound(); /* set boundary conditions */

	for (iter = 0; iter <= OUTER_ITER; iter++) { /* outer iteration loop */
		Tcoeff(aE, aW, aN, aS, aP, b);
		for (iter_T = 0; iter_T <= T_ITER; iter_T++)
			solveSOR(T, b, aE, aW, aN, aS, aP);
		printConv();
	} /* for outer iteration loop */

	output();

    system("pause"); /* prevent the command window  from automatically closing */
	return 0;

} /* main */


/* ################################################################# */
void grid(void)
/* ################################################################# */
{
/***** Purpose: Defining the geometrical variables ******/
/*****          See fig. 6.2-6.4 in ref. 1 ******/
	int	I, J, i, j;
	double	Dx, Dy;

	/* Length of volume element */

	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;

	/* Length variable for the scalar points in the x direction */

	x[0] = 0.;
	x[1] = 0.5*Dx;

	for (I = 2; I <= NPI; I++)
		x[I] = x[I-1] + Dx;

	x[NPI+1] = x[NPI] + 0.5*Dx;

	/* Length variable for the scalar points fi[i][j] in the y direction */

	y[0] = 0.;
	y[1] = 0.5*Dy;

	for (J = 2; J <= NPJ; J++)
		y[J] = y[J-1] + Dy;

	y[NPJ+1] = y[NPJ] + 0.5*Dy;

	/* Length variable for the velocity components u[i][j] in the x direction */

	x_u[0] = 0.;
	x_u[1] = 0.;

	for (i = 2; i <= NPI + 1; i++)
		x_u[i] = x_u[i-1] + Dx;


	/* Length variable for the velocity components v[i][j] in the y direction */

	y_v[0] = 0.;
	y_v[1] = 0.;
	for (j = 2; j <= NPJ + 1; j++)
		y_v[j] = y_v[j-1] + Dy;

} /* grid */

/* ################################################################# */
void init(void)
/* ################################################################# */
{
/***** Purpose: To initilise all parameters. ******/
	int	I, J, i, j, NFI;

	memalloc();
	grid();

	/* Initialising all variables  */

	omega = 0.9; /* Over-relaxation factor for solver */

	for (I = 0; I <= NPI + 1; I++) {
		i = I;
		for (J = 0; J <= NPI + 1; J++) {
			j = J;

			u    [i][J] = 1.;    /* Velocity in x-direction */
			v    [I][j] = 4.;    /* Velocity in y-direction */
			T    [I][J] = 0.;    /* Temperature */
			rho  [I][J] = 1.;    /* Density */
			mu   [I][J] = 0.;    /* Turbulent viscosity ( = 10*mu_laminar) */
			Gamma[I][J] = 1.;    /* Thermal conductivity */
			Cp   [I][J] = 1.;    /* J/(K*kg) Heat capacity - assumed constant for this problem */
		} /* for J */
	} /* for I */

	/* Definition of first internal node for the different fi variables. */
	/* See fig. 9.1. */

	Istart = 1;
	Iend   = NPI;
	Jstart = 1;
	Jend   = NPJ;

} /* init */

/* ################################################################# */
void bound(void)
/* ################################################################# */
{
/***** Purpose: Specify the boundary conditions *****/
	int I, J;

	/* Fixed temperature at the upper and lower wall */

	for (I = 0; I <= NPI + 1; I++) {
		T[I][0]     = 100.; /* Temperature in Kelvin */
		T[I][NPJ+1] = 0.;   /* Temperature in Kelvin */
	} /* for I */

	/* Fixed temperature at the left and right wall */

	for (J = 0; J <= NPJ; J++) {
		T[0][J] = 100.;   /* Temperature in Kelvin */
		T[NPI+1][J] = 0.; /* Temperature in Kelvin */
	} /* for J */

} /* bound */


/* ################################################################# */
void output(void)
/* ################################################################# */
{
/***** Purpose: print out result table to file ******/
	int	I, J;
	FILE *fp;

	fp = fopen("output.dat", "w");

	for (I = 1; I <= NPI; I++) {
		for (J = 1; J <= NPJ; J++)
			fprintf(fp, "%e\t%e\t%e\n", x[I], y[J], T[I][J]);
		fprintf(fp, "\n");
	} /* for I */

	fclose(fp);

} /* output */

/* ################################################################# */
void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Successive Over-relaxation ******/
	int	I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           * omega /aP[I][J] - (omega - 1)*fi[I][J];

}

/* ################################################################# */
void conv(void)
/* ################################################################# */
{
/***** Purpose: To calculate the convective mass flux component pr. unit ******/
/*****          area defined in eq. 5.7 ******/
	int	I, J, i, j;

	for (I = 1; I <= NPI + 1; I++) {
		i = I;
		for (J = 1; J <= NPJ + 1; J++) {
			j = J;
			F_u[i][J] = (rho[I-1][J  ]*(x[I] - x_u[i]) + rho[I][J]*(x_u[i] - x[I-1]))*u[i][J]/(x[I] - x[I-1]); /* = F(i, J) */
			F_v[I][j] = (rho[I  ][J-1]*(y[J] - y_v[j]) + rho[I][J]*(y_v[j] - y[J-1]))*v[I][j]/(y[J] - y[J-1]); /* = F(I, j) */
		} /* for J */
	} /* for I */

} /* conv */

/* ################################################################# */
void printConv(void)
/* ################################################################# */
{
	fprintf(stdout, "%3d T[%d][%d]= %7.3f T[%d][%d]= %7.3f T[%d][%d]= %7.3f T[%d][%d]= %7.3f\n", 
	iter, 
    NPI/4+1, NPJ/4+1, T[NPI/4+1][NPJ/4+1],           // point A
    3*NPI/4+1, NPJ/4+1, T[3*NPI/4+1][NPJ/4+1],       // point B
    NPI/4+1, 3*NPJ/4+1, T[NPI/4+1][3*NPJ/4+1],       // point C
    3*NPI/4+1, 3*NPJ/4+1, T[3*NPI/4+1][3*NPJ/4+1]);  // point D
}

/* ################################################################# */
void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/
	int	i, j, I, J;
	double	Fw, Fe, Fs, Fn, Pw, Pe, Ps, Pn, 
		Dw, De, Ds, Dn, 
		AREAw, AREAe, AREAs, AREAn;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*Cp[I][J]*AREAw;
			Fe = F_u[i+1][J  ]*Cp[I][J]*AREAe;
			Fs = F_v[I  ][j  ]*Cp[I][J]*AREAs;
			Fn = F_v[I  ][j+1]*Cp[I][J]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = 0.5*(Gamma[I-1][J  ] + Gamma[I  ][J  ])/(x[I  ] - x[I-1])*AREAw;
			De = 0.5*(Gamma[I  ][J  ] + Gamma[I+1][J  ])/(x[I+1] - x[I  ])*AREAe;
			Ds = 0.5*(Gamma[I  ][J-1] + Gamma[I  ][J  ])/(y[J  ] - y[J-1])*AREAs;
			Dn = 0.5*(Gamma[I  ][J  ] + Gamma[I  ][J+1])/(y[J+1] - y[J  ])*AREAn;

			/* The source terms */

			SP[I][J] = -2.*AREAw*AREAs; // "-b*T" part of the source term //
			Su[I][J] = 10.*AREAw*AREAs; // "a" part of the source term //

			/* The coefficients (central differencing sheme) */

			aW[I][J] = Dw + Fw/2;
			aE[I][J] = De - Fe/2;
			aS[I][J] = Ds + Fs/2;
			aN[I][J] = Dn - Fn/2;

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J];

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J];

			/* now the Thomas algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */

} /* Tcoeff */

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
double **double_2D_matrix (int nm, int np)
/* ################################################################# */
{
/* create an 2D matrix with size [nm, np] of type double */
	int i;
	double **m;

	m = (double **) calloc(nm, sizeof(double *));
	for ( i = 0; i < nm; i++)
		m[i] = (double *) calloc(np, sizeof(double));

	return m;

} /* double_2D_matrix */

/* ################################################################# */
void memalloc(void)
/* ################################################################# */
{
	x    = double_1D_array(NPI + 2);
	x_u  = double_1D_array(NPI + 2);
	y    = double_1D_array(NPJ + 2);
	y_v  = double_1D_array(NPJ + 2);

	u      = double_2D_matrix(NPI + 2, NPJ + 2);
	v      = double_2D_matrix(NPI + 2, NPJ + 2);
	T      = double_2D_matrix(NPI + 2, NPJ + 2);
	rho    = double_2D_matrix(NPI + 2, NPJ + 2);
	mu     = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma  = double_2D_matrix(NPI + 2, NPJ + 2);
	Cp     = double_2D_matrix(NPI + 2, NPJ + 2);

	aP     = double_2D_matrix(NPI + 2, NPJ + 2);
	aE     = double_2D_matrix(NPI + 2, NPJ + 2);
	aW     = double_2D_matrix(NPI + 2, NPJ + 2);
	aN     = double_2D_matrix(NPI + 2, NPJ + 2);
	aS     = double_2D_matrix(NPI + 2, NPJ + 2);
	b      = double_2D_matrix(NPI + 2, NPJ + 2);

	SP     = double_2D_matrix(NPI + 2, NPJ + 2);
	Su     = double_2D_matrix(NPI + 2, NPJ + 2);

	F_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	F_v    = double_2D_matrix(NPI + 2, NPJ + 2);

} /* memalloc */
